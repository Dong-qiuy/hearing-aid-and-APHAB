rm(list = ls())
library(dplyr)
library(tidyr)
library(mmrm)
library(emmeans)
load("data/final data.RDATA")
#####################analysis
var_primary <- c("aphab", "aphab_ec", "aphab_bn", "aphab_rv", "aphab_av")

var_scondary <- c("hhie_situation", "anjing_shengmu","anjing_yunmu",
                  "anjing_shengdiao","anjing_danyinjie","zaosheng_juzi")

var_all <- c(var_primary, var_scondary)

var_lab <- c("APHAB","APHAB-ease of communication","APHAB-background noise",
             "APHAB-reverberation","APHAB-aversiveness",
             "HHIE","HHIE emotion","HHIE social",
             "anjing_shengmu","anjing_yunmu","anjing_shengdiao",
             "anjing_danyinjie","zaosheng_juzi")

run_mmrm <- function(data, var, pairwise = FALSE) {
  data[["baseline_val"]] <- data[[paste0(var, "_t0")]]
  data[["m1_val"]] <- data[[paste0(var, "_t1")]]
  data[["m3_val"]] <- data[[paste0(var, "_t2")]]
  data[["m12_val"]] <- data[[paste0(var, "_t3")]]

  # 2. Calculate Pooled Baseline SD AND Baseline DF
  # Using the shortcut method (fitting a simple linear model at baseline)
  # This perfectly extracts both the pooled SD and its associated degrees of freedom.
  baseline_lm <- lm(baseline_val ~ group, data = data)
  pooled_sd_base <- sigma(baseline_lm)
  baseline_df <- df.residual(baseline_lm) # This gives N - 3

  long_data <- data |>
    mutate(
      chg_1m = m1_val - baseline_val,
      chg_3m = m3_val - baseline_val,
      chg_12m = m12_val - baseline_val
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

  # Fit MMRM
  fit <- mmrm(CHG ~ VISIT * group + baseline_val + us(VISIT|id), data = long_data)

  if(!pairwise){
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
  }else{
    # Extract Dunnett-adjusted contrasts at 3 months
    emm <- emmeans(fit, ~ group | VISIT, adjust = "tukey")
    pairwise_contrasts <- pairs(emm)
    res_raw <- summary(pairwise_contrasts, infer = c(TRUE, TRUE)) |>
      as.data.frame()

    # Standardize using the fixed baseline SD
    res_eff <- eff_size(emm,
                        sigma = pooled_sd_base,
                        edf = baseline_df)

    res_eff_df <- as.data.frame(res_eff) |>
      select(contrast, VISIT, effect.size, eff_SE = SE,
             eff_lower.CL = lower.CL, eff_upper.CL = upper.CL)
  }


  # Combine Raw and Effect Size data
  final_res <- inner_join(res_raw, res_eff_df, by = c("contrast", "VISIT"))
  return(final_res)
}

# Dunnett
r_dunet <- sapply(seq_along(var_all), function(i){
  r <- run_mmrm(data = dat, var = var_all[i])
  r$variable <- var_lab[i]
  return(r)
}, simplify = FALSE)

r_dunet <- bind_rows(r_dunet)

# Pairwise
r_tukey <- sapply(seq_along(var_all), function(i){
  r <- run_mmrm(data = dat, var = var_all[i], pairwise = TRUE)
  r$variable <- var_lab[i]
  return(r)
}, simplify = FALSE)

r_tukey <- bind_rows(r_tukey)
