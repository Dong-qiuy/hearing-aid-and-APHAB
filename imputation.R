rm(list = ls())
library(mice)
library(dplyr)
library(tidyr)
library(mmrm)
library(emmeans)
library(ggplot2)
load("data/final data.RDATA")


dfm <- dat[,c("id", "group", "aphab_t0", "aphab_t1", "aphab_t2", "aphab_t3")]
dfm$group <- relevel(dfm$group, ref = "Self-fitting")

dfm$is_trt <- as.numeric(dfm$group != "Self-fitting")

# ------------------------------------------------------------------------------
# 2. Custom MICE Imputation Function (Embedded Conditional MNAR)
# ------------------------------------------------------------------------------
# This function applies the delta shift inside the chained equations, ensuring
# that penalized early dropouts influence the predictions of their later visits.
mice.impute.mnar.custom <- function(y, ry, x, wy = NULL, delta = 0, ...) {
  if (is.null(wy)) wy <- !ry

  # 1. Perform standard MAR imputation (Bayesian linear regression)
  imp <- mice.impute.norm(y, ry, x, wy = wy, ...)

  # 2. Identify treatment dropouts
  # x is the predictor matrix. We check our 'is_trt' column for the missing rows (wy)
  trt_dropouts <- x[wy, "is_trt"] == 1

  # 3. Apply the conditional penalty (Shift)
  imp[trt_dropouts] <- imp[trt_dropouts] + delta

  # 4. Enforce logical boundaries (Squeeze) inside the loop
  imp <- pmax(0, pmin(100, imp))

  return(imp)
}

# ------------------------------------------------------------------------------
# 3. Configure MICE (Forward-Only Predictor Matrix) ================
# ------------------------------------------------------------------------------

# Assume 'data_wide' has columns: USUBJID, ARM, aphab_t0, aphab_t1, aphab_t2, aphab_t3
init <- mice(dfm, maxit = 0)
meth <- init$method
pred <- init$predictorMatrix

visit_cols <- c("aphab_t1", "aphab_t2", "aphab_t3")


# Set MNAR method for all post-baseline visits
meth[visit_cols] <- "mnar.custom"

pred[, "id"] <- 0
pred["id", ] <- 0
pred["group", ] <- 0
pred["aphab_t0", ] <- 0

# is_trt predicts the visits, but is not predicted by anything
pred[, "is_trt"] <- 0
pred[visit_cols, "is_trt"] <- 1

# Baseline and group predict everything
pred[visit_cols, "aphab_t0"] <- 1
pred[visit_cols, "group"]      <- 1

# STRICT CHRONOLOGICAL CONSTRAINTS:
# Baseline is NOT predicted by any future visits
pred["aphab_t0", visit_cols] <- 0

# 1 month is NOT predicted by 3m or 12m
pred["aphab_t1", c("aphab_t2", "aphab_t3")] <- 0

# 3 month is NOT predicted by 12m
pred["aphab_t2", "aphab_t3"] <- 0

print(pred)


# ------------------------------------------------------------------------------
# 4. Analysis Function (MMRM on Change Scores) =========================
# ------------------------------------------------------------------------------
run_mmrm_analysis <- function(wide_data) {
  long_data <- wide_data |>
    mutate(
      chg_1m = aphab_t0 - aphab_t1,
      chg_3m = aphab_t0 - aphab_t2,
      chg_12m = aphab_t0 - aphab_t3
    ) |>
    pivot_longer(
      cols = starts_with("chg_"),
      names_to = "VISIT",
      values_to = "CHG",
      names_prefix = "chg_"
    ) |>
    mutate(
      VISIT = as.factor(VISIT),
      id = as.factor(id)
    )

  # All using the original data
  baseline_lm <- lm(aphab_t0 ~ group, data = dfm)
  pooled_sd_base <- sigma(baseline_lm)
  baseline_df <- df.residual(baseline_lm) # This gives N - 3

  # Fit MMRM (aphab_t0 is used as a covariate here)
  fit <- mmrm(CHG ~ VISIT * group + aphab_t0 + us(VISIT|id), data = long_data)

  # Extract Dunnett-adjusted contrasts at 3 months
  emm <- emmeans(fit, trt.vs.ctrl ~ group | VISIT, adj = "dunnett")
  res_raw <- summary(emm, infer = c(TRUE, TRUE))$contrasts |>
    as.data.frame()

  # Standardize using the fixed baseline SD
  res_eff <- eff_size(emm$emmeans,
                      method = "trt.vs.ctrl",
                      sigma = pooled_sd_base,
                      edf = baseline_df,
                      adj = "dunnett")

  res_eff_df <- as.data.frame(res_eff) |>
    select(contrast, VISIT, effect.size, eff_SE = SE,
           eff_lower.CL = lower.CL, eff_upper.CL = upper.CL)

  # Combine Raw and Effect Size data
  final_res <- inner_join(res_raw, res_eff_df, by = c("contrast", "VISIT"))
  return(final_res)
}


# ------------------------------------------------------------------------------
# 5. Tipping Point Loop (Sensitivity Analysis) ==========================
# ------------------------------------------------------------------------------
deltas <- seq(0, 100, by = 10)
all_tp_results <- list()

all_tp_results <- sapply(deltas, function(d){

  message(sprintf("Processing Delta: %.1f", d))

  curr_blots <- list(
    aphab_t1  = list(delta = d),   # penalt
    aphab_t2  = list(delta = d),  # penalt
    aphab_t3 = list(delta = d)   # penalt
  )

  imp <- mice(
    data = dfm,
    method = meth,
    predictorMatrix = pred,
    blots = curr_blots,
    m = 20,
    seed = 20260218,
    print = FALSE
  )

  m_results <- list()
  for (i in 1:imp$m) {
    m_results[[i]] <- run_mmrm_analysis(complete(imp, i))
  }

  pooled_df <- bind_rows(m_results) %>%
    group_by(contrast, VISIT) %>%
    summarize(
      delta = d,
      # Between-imputation variance
      B = var(estimate),
      # Within-imputation variance
      W = mean(SE^2),
      est = mean(estimate),
      se  = sqrt(W + (1 + 1/imp$m) * B),
      eff_est = mean(effect.size),
      eff_se  = sqrt(mean(eff_SE^2) + (1 + 1/imp$m) * var(effect.size)),
      .groups = "drop"
    ) %>%
    mutate(
      # Rubin's Rules Degrees of Freedom (Barnard-Rubin adjustment)
      r = (1 + 1/imp$m) * B / W,
      df_mice = (imp$m - 1) * (1 + 1/r)^2,
      # MULTIPLICITY ADJUSTMENT (Bonferroni-adjusted alpha for 2 treatment arms)
      # Alpha = 0.05 / 2 = 0.025 for a two-tailed 95% CI
      alpha_adj = 0.05 / 2,
      t_crit = qt(1 - alpha_adj/2, df = df_mice),
      lower = est - t_crit * se,
      upper = est + t_crit * se,
      p_val = 2 * (1 - pt(abs(est / se), df = df_mice)),
      p_val = 2 * (1 - pnorm(abs(est / se))),
      eff_lower = eff_est - t_crit * eff_se,
      eff_upper = eff_est + t_crit * eff_se
    )

  return(pooled_df)

}, simplify = FALSE)


tipping_table <- bind_rows(all_tp_results)

tipping_table$Visit <- factor(tipping_table$VISIT,
                              levels = c("1m", "3m", "12m"),
                              labels = c("1 M", "3 M", "12 M"),
                              ordered = TRUE)

# ------------------------------------------------------------------------------
# 6. Visualization - Plot 1: Mean Differences ==================================
# ------------------------------------------------------------------------------
clin_sig_threshold <- 17

diff_plot <- ggplot(tipping_table,
                    aes(x = delta, y = est, color = Visit, fill = Visit)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  #geom_hline(yintercept = clin_sig_threshold, linetype = "dotted", color = "darkgreen") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = p_val < 0.05), size = 2) +
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 19),
                     name = "Significant (p < 0.05)") +
  facet_wrap(~contrast) +
  # scale_x_reverse() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Multiple imputation on APHAB: Mean Difference",
       subtitle = "Forward-only Imputation | Using Satterthwaite adjusted DFs from MMRM",
       x = "Delta (Penalty applied to Treatment)",
       y = "Mean Difference vs Control") +
  theme_minimal()

# ------------------------------------------------------------------------------
# 7. Visualization - Plot 2: Cohen's d (Effect Size) ===========================
# ------------------------------------------------------------------------------
eff_plot <- ggplot(tipping_table,
                   aes(x = delta, y = eff_est, color = Visit, fill = Visit)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_ribbon(aes(ymin = eff_lower, ymax = eff_upper), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = p_val < 0.05), size = 2) +
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 19),
                     name = "Significant (p < 0.05)") +
  facet_wrap(~contrast) +
  # scale_x_reverse() +
  scale_fill_manual(values = c('#00008B',"#E21A1C","#f1b656"))+
  scale_color_manual(values = c('#00008B',"#E21A1C","#f1b656"))+
  labs(title = "Multiple imputation on APHAB: Cohen's d (Effect Size)",
       subtitle = "Standardized by Unaided Baseline SD",
       x = "Delta (Penalty applied to Treatment)",
       y = "Cohen's d [95% CI]") +
  theme_minimal()

print(eff_plot)
print(tipping_table)
