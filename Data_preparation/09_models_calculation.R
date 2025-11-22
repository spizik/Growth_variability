# Model formula ####
form_main <- log_cv_RWI ~ 
  
  species * (median_age + 
             median_range + 
             mid_TRW + 
             mean_temp + 
             mean_cwb + 
             sox + 
             pH_L1 +
             C.N_FH
  ) +
  species * (median_age:mean_temp +
               median_age:mean_cwb +
               median_range:mean_temp +
               median_range:mean_cwb +
             pH_L1 : mid_TRW +
             pH_L1 : mean_temp +
             sox : mid_TRW +
             sox : mean_temp)

form_species <- log_cv_RWI ~ 
  
  median_age + 
  median_range + 
  mid_TRW + 
  mean_temp + 
  mean_cwb + 
  sox + 
  pH_L1 +
  C.N_FH +
  median_age:mean_temp +
  median_age:mean_cwb +
  median_range:mean_temp +
  median_range:mean_cwb +
  pH_L1 : mid_TRW +
  pH_L1 : mean_temp +
  sox : mid_TRW +
  sox : mean_temp

## Individual models calculations ####
mod_main <- lme(form_main,    random = ~ 1 | year, correlation = corARMA(p = 1, q = 0), data = main_mod_dataset_scaled, na.action = na.omit)
mod_abal <- lme(form_species, random = ~ 1 | year, correlation = corARMA(p = 1, q = 0), data = abal_mod_dataset_scaled, na.action = na.omit)
mod_pcab <- lme(form_species, random = ~ 1 | year, correlation = corARMA(p = 1, q = 0), data = pcab_mod_dataset_scaled, na.action = na.omit)
mod_pisy <- lme(form_species, random = ~ 1 | year, correlation = corARMA(p = 1, q = 0), data = pisy_mod_dataset_scaled, na.action = na.omit)
mod_fasy <- lme(form_species, random = ~ 1 | year, correlation = corARMA(p = 1, q = 0), data = fasy_mod_dataset_scaled, na.action = na.omit)
mod_qusp <- lme(form_species, random = ~ 1 | year, correlation = corARMA(p = 1, q = 0), data = qusp_mod_dataset_scaled, na.action = na.omit)

## Saving models ####
saveRDS(mod_main, "Calculated_datasets/Calculated_models/model_MAIN.rds")
saveRDS(mod_abal, "Calculated_datasets/Calculated_models/model_ABAL.rds")
saveRDS(mod_fasy, "Calculated_datasets/Calculated_models/model_FASY.rds")
saveRDS(mod_pcab, "Calculated_datasets/Calculated_models/model_PCAB.rds")
saveRDS(mod_pisy, "Calculated_datasets/Calculated_models/model_PISY.rds")
saveRDS(mod_qusp, "Calculated_datasets/Calculated_models/model_QUSP.rds")
