library(mmrm)
library(emmeans)
library(tidyr)
library(dplyr)
rm(list = ls())
load("data/final data.RDATA")
table(dat$center)
dat$city<-ifelse(dat$center %in% c(1,2,3,4,5,6,7),1,ifelse(
  dat$center %in% c(8,9),2,ifelse(
    dat$center==10,3,NA
  )
))
var_primary <- c("aphab", "aphab_ec", "aphab_bn", "aphab_rv", "aphab_av")
var_scondary <- c("hhie","hhie_situation","hhie_emotion", "anjing_shengmu","anjing_yunmu",
                  "anjing_shengdiao","zaosheng_juzi","ioi_total")
var_all <- c(var_primary, var_scondary)

var_lab <- c("APHAB","APHAB-ease of communication","APHAB-background noise",
             "APHAB-reverberation","APHAB-aversiveness",
             "HHIE","HHIE social","HHIE emotion",
             "anjing_shengmu","anjing_yunmu","anjing_shengdiao",
             "zaosheng_juzi","ioi")
dat$center<-as.character(dat$center)



# Dunnett
r_dunet <- sapply(seq_along(var_all), function(i){
  r <- run_mmrm(data = dat, var = var_all[i])
  r$contrasts$variable <- var_lab[i]
  r$emm_means$variable <- var_lab[i]
  return(r)
}, simplify = FALSE)

all_contrasts <- bind_rows(lapply(r_dunet, `[[`, "contrasts"))
all_emm_means <- bind_rows(lapply(r_dunet, `[[`, "emm_means"))
names(all_emm_means)
all_contrasts$difference<-sprintf("%0.2f (%0.2f, %0.2f)",
                                  all_contrasts$estimate,all_contrasts$lower.CL,all_contrasts$upper.CL)
all_contrasts$effectsize<-sprintf("%0.2f (%0.2f, %0.2f)",
                                  all_contrasts$effect.size,all_contrasts$eff_lower.CL,all_contrasts$eff_upper.CL)
all_emm_means$meansd<-sprintf("%0.2f (%0.2f)",
                              all_emm_means$emmean,all_emm_means$SE)
write.csv(all_contrasts,file = "result/dunnett_result.csv",row.names = F)
write.csv(all_emm_means,file = "result/emmean_result.csv",row.names = F)

# Pairwise
r_tukey <- sapply(seq_along(var_all), function(i){
  r <- run_mmrm(data = dat, var = var_all[i], pairwise = TRUE)
  r$contrasts$variable <- var_lab[i]
  return(r)
}, simplify = FALSE)

all_contrasts <- bind_rows(lapply(r_tukey, `[[`, "contrasts"))
all_contrasts$difference<-sprintf("%0.2f (%0.2f, %0.2f)",
                                  all_contrasts$estimate,all_contrasts$lower.CL,all_contrasts$upper.CL)
all_contrasts$effectsize<-sprintf("%0.2f (%0.2f, %0.2f)",
                                  all_contrasts$effect.size,all_contrasts$eff_lower.CL,all_contrasts$eff_upper.CL)
write.csv(all_contrasts,file = "result/tukey_result.csv",row.names = F)

########sensitivity analysis
# Dunnett
r_dunet <- sapply(seq_along(var_all), function(i){
  r <- run_mmrm(data = dat, var = var_all[i])
  r$contrasts$variable <- var_lab[i]
  r$emm_means$variable <- var_lab[i]
  return(r)
}, simplify = FALSE)

all_contrasts <- bind_rows(lapply(r_dunet, `[[`, "contrasts"))
all_emm_means <- bind_rows(lapply(r_dunet, `[[`, "emm_means"))

all_contrasts$difference<-sprintf("%0.2f (%0.2f, %0.2f)",
                                  all_contrasts$estimate,all_contrasts$lower.CL,all_contrasts$upper.CL)
all_contrasts$effectsize<-sprintf("%0.2f (%0.2f, %0.2f)",
                                  all_contrasts$effect.size,all_contrasts$eff_lower.CL,all_contrasts$eff_upper.CL)
all_emm_means$meansd<-sprintf("%0.2f (%0.2f)",
                              all_emm_means$emmean,all_emm_means$SE)
write.csv(all_contrasts,file = "result/dunnett_result(sensitivity analysis).csv",row.names = F)
write.csv(all_emm_means,file = "result/emmean_result(sensitivity analysis).csv",row.names = F)







