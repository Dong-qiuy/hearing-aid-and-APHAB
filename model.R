run_mmrm <- function(data, var, pairwise = FALSE) {
  if(var=="ioi_total"){
    data[["m1_val"]] <- data[[paste0(var, "_t1")]]
    data[["m3_val"]] <- data[[paste0(var, "_t2")]]
    data[["m12_val"]] <- data[[paste0(var, "_t3")]]
    baseline_lm <- lm(m1_val ~ group, data = data)
    
  }else{
    data[["baseline_val"]] <- data[[paste0(var, "_t0")]]
    data[["m1_val"]] <- data[[paste0(var, "_t1")]]
    data[["m3_val"]] <- data[[paste0(var, "_t2")]]
    data[["m12_val"]] <- data[[paste0(var, "_t3")]]
    baseline_lm <- lm(baseline_val ~ group, data = data)
  }
  # 2. Calculate Pooled Baseline SD AND Baseline DF
  # Using the shortcut method (fitting a simple linear model at baseline)
  # This perfectly extracts both the pooled SD and its associated degrees of freedom. # This gives N - 3
  pooled_sd_base <- sigma(baseline_lm)
  baseline_df <- df.residual(baseline_lm)
  if (var %in% c("aphab", "aphab_ec", "aphab_bn", "aphab_rv", "aphab_av")){
    long_data <- data |>
      mutate(
        chg_1m = baseline_val - m1_val,
        chg_3m = baseline_val - m3_val,
        chg_12m = baseline_val - m12_val
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
  }else if (var=="ioi_total"){
    long_data <- data |>
      mutate(
        chg_1m = m1_val,
        chg_3m = m3_val,
        chg_12m = m12_val
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
  }else{
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
  }
  if(!pairwise){
    if (var=="ioi_total"){
      fit <- mmrm(CHG ~ VISIT * group+city + us(VISIT|id), data = long_data)
    }else{
      fit <- mmrm(CHG ~ VISIT * group + baseline_val+city + us(VISIT|id), data = long_data)
    }
    # Extract Dunnett-adjusted contrasts at 3 months
    emm <- emmeans(fit, trt.vs.ctrl ~ group | VISIT, adj = "dunnett")
    res_raw <- summary(emm, infer = c(TRUE, TRUE))$contrasts |>
      as.data.frame()
    
    # Standardize using the fixed baseline SD
    res_eff <- eff_size(emm$emmeans,
                        method = "trt.vs.ctrl",
                        sigma = pooled_sd_base,
                        edf = baseline_df,
                        adj = "dunnett",
                        SE.sigma = 0)
    
    res_eff_df <- as.data.frame(res_eff) |>
      select(contrast, VISIT, effect.size, eff_SE = SE,
             eff_lower.CL = lower.CL, eff_upper.CL = upper.CL)
    res<-as.data.frame(emm$emmeans)
    # Combine Raw and Effect Size data
    final_res <- inner_join(res_raw, res_eff_df, by = c("contrast", "VISIT"))
    return(list(
      contrasts = final_res,
      emm_means = res
    ))
  }else{
    long_data$group <- factor(long_data$group,
                              levels = c("CHW-support","Remote-support", "Self-fitting"))
    if (var=="ioi_total"){
      fit <- mmrm(CHG ~ VISIT * group+city + us(VISIT|id), data = long_data)
    }else{
      fit <- mmrm(CHG ~ VISIT * group+city + baseline_val + us(VISIT|id), data = long_data)
    }
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
    # Combine Raw and Effect Size data
    final_res <- inner_join(res_raw, res_eff_df, by = c("contrast", "VISIT"))
    return(list(
      contrasts = final_res
    ))
  }
  
  
}

############sensitivity analysis

run_mmrm <- function(data, var, pairwise = FALSE) {
  if(var=="ioi_total"){
    data[["m1_val"]] <- data[[paste0(var, "_t1")]]
    data[["m3_val"]] <- data[[paste0(var, "_t2")]]
    data[["m12_val"]] <- data[[paste0(var, "_t3")]]
    baseline_lm <- lm(m1_val ~ group, data = data)
    
  }else{
    data[["baseline_val"]] <- data[[paste0(var, "_t0")]]
    data[["m1_val"]] <- data[[paste0(var, "_t1")]]
    data[["m3_val"]] <- data[[paste0(var, "_t2")]]
    data[["m12_val"]] <- data[[paste0(var, "_t3")]]
    baseline_lm <- lm(baseline_val ~ group, data = data)
  }
  # 2. Calculate Pooled Baseline SD AND Baseline DF
  # Using the shortcut method (fitting a simple linear model at baseline)
  # This perfectly extracts both the pooled SD and its associated degrees of freedom. # This gives N - 3
  pooled_sd_base <- sigma(baseline_lm)
  baseline_df <- df.residual(baseline_lm)
  if (var %in% c("aphab", "aphab_ec", "aphab_bn", "aphab_rv", "aphab_av")){
    long_data <- data |>
      mutate(
        chg_1m = baseline_val - m1_val,
        chg_3m = baseline_val - m3_val,
        chg_12m = baseline_val - m12_val
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
  }else if (var=="ioi_total"){
    long_data <- data |>
      mutate(
        chg_1m = m1_val,
        chg_3m = m3_val,
        chg_12m = m12_val
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
  }else{
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
  }
  if(!pairwise){
    if (var=="ioi_total"){
      fit <- mmrm(CHG ~ VISIT * group+city+education_c2+income_c+duju+iphone + us(VISIT|id), data = long_data)
    }else{
      fit <- mmrm(CHG ~ VISIT * group + baseline_val+city+education_c2+income_c+duju+iphone + us(VISIT|id), data = long_data)
    }
    # Extract Dunnett-adjusted contrasts at 3 months
    emm <- emmeans(fit, trt.vs.ctrl ~ group | VISIT, adj = "dunnett")
    res_raw <- summary(emm, infer = c(TRUE, TRUE))$contrasts |>
      as.data.frame()
    
    # Standardize using the fixed baseline SD
    res_eff <- eff_size(emm$emmeans,
                        method = "trt.vs.ctrl",
                        sigma = pooled_sd_base,
                        edf = baseline_df,
                        adj = "dunnett",
                        SE.sigma = 0)
    
    res_eff_df <- as.data.frame(res_eff) |>
      select(contrast, VISIT, effect.size, eff_SE = SE,
             eff_lower.CL = lower.CL, eff_upper.CL = upper.CL)
    res<-as.data.frame(emm$emmeans)
  }else{
    long_data$group <- factor(long_data$group,
                              levels = c("CHW-support","Remote-support", "Self-fitting"))
    if (var=="ioi_total"){
      fit <- mmrm(CHG ~ VISIT * group+city+education_c2+income_c+duju+iphone + us(VISIT|id), data = long_data)
    }else{
      fit <- mmrm(CHG ~ VISIT * group + baseline_val+city+education_c2+income_c+duju+iphone + us(VISIT|id), data = long_data)
    }
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
  return(list(
    contrasts = final_res,
    emm_means = res
  ))
}