library(caret)
library(lme4)
library(DHARMa)
library(lmtest)
library(influence.ME)

## test overfittingu ####
test.overfitting <- function(in_dataset, form, nrep=100){
  # in_dataset = main_mod_dataset_scaled
  # form = form_main
  # nrep = 100
  
  years <- unique(in_dataset$year)
  results <- data.frame(replication = c(1:nrep), RMSE = NA)
  
  for (i in results$replication) {
    # print(i)
    
    set.seed(123)  # pro reprodukovatelnost
    
    # vytvoří logický vektor TRUE (trénink) / FALSE (validace)
    # print("vyber indexu")
    if(length(unique(in_dataset$species == 1))){
      train_idx <- sample(seq_len(nrow(in_dataset)), size = 0.5 * nrow(in_dataset))
    } else {
      train_idx <- in_dataset %>%
        group_by(species) %>%
        sample_frac(0.5) %>%
        ungroup() %>%
        pull(row_number())
    }
    
    # print("inicializace datasetu")
    train_data <- in_dataset[train_idx, ]
    test_data  <- in_dataset[-train_idx, ]
    
    # Fit model na trénovacích datech
    # print("model calc")
    # environment(form) <- environment()
    tryCatch({
      model_cv <- lme(
        form,
        random = ~ 1 | year,
        correlation = corARMA(p = 1, q = 0),
        data = train_data
      )
      
      # Predikce na testovacích datech
      # print("predikce")
      preds <- predict(model_cv, newdata = test_data, level = 0)  # only fixed effects
      
      # Výpočet RMSE
      # print("RSME")
      rmse <- sqrt(mean((test_data$log_cv_RWI - preds)^2, na.rm = TRUE))
      results$RMSE[i] <- rmse
    }, error = function(e) {
      # message(sprintf("Chyba v iteraci %d: %s", i, e$message))
      results$RMSE[i] <- NA  # Zapsat NA pro selhanou iteraci
    })
  }
  
  # Průměrná chyba
  output <- data.frame(RSME_validace = mean(results$RMSE, na.rm=T),
                       SD_RSME = sd(results$RMSE, na.rm = T),
                       SD_dastasetu = sd(in_dataset$log_cv_RWI, na.rm = TRUE))
  
  return(output)
}

print("---------------------------------------------")
print("Overfitting test - all data")
test.overfitting(main_mod_dataset_scaled, form_main)

print("---------------------------------------------")
print("Overfitting test - ABAL")
test.overfitting(abal_mod_dataset_scaled, form_species)

print("---------------------------------------------")
print("Overfitting test - PCAB")
test.overfitting(pcab_mod_dataset_scaled, form_species)

print("---------------------------------------------")
print("Overfitting test - PISY")
test.overfitting(pisy_mod_dataset_scaled, form_species)

print("---------------------------------------------")
print("Overfitting test - FASY")
test.overfitting(fasy_mod_dataset_scaled, form_species)

print("---------------------------------------------")
print("Overfitting test - QUSP")
test.overfitting(qusp_mod_dataset_scaled, form_species)

print("---------------------------------------------")
print("---------------------------------------------")
print("---------------------------------------------")

## test qqplots ####
par(mfrow=c(1,2))

print("---------------------------------------------")
print("QQplot a Disperze - all data")
in_model = mod_main

qqnorm(resid(in_model)) ; qqline(resid(in_model))
plot(fitted(in_model), resid(in_model)) ; abline(h = 0, col = "red")

print("---------------------------------------------")
print("QQplot a Disperze - ABAL")
in_model = mod_abal

qqnorm(resid(in_model)) ; qqline(resid(in_model))
plot(fitted(in_model), resid(in_model)) ; abline(h = 0, col = "red")

print("---------------------------------------------")
print("QQplot a Disperze - PCAB")
in_model = mod_pcab

qqnorm(resid(in_model)) ; qqline(resid(in_model))
plot(fitted(in_model), resid(in_model)) ; abline(h = 0, col = "red")

print("---------------------------------------------")
print("QQplot a Disperze - PISY")
in_model = mod_pisy

qqnorm(resid(in_model)) ; qqline(resid(in_model))
plot(fitted(in_model), resid(in_model)) ; abline(h = 0, col = "red")

print("---------------------------------------------")
print("QQplot a Disperze - FASY")
in_model = mod_fasy

qqnorm(resid(in_model)) ; qqline(resid(in_model))
plot(fitted(in_model), resid(in_model)) ; abline(h = 0, col = "red")

print("---------------------------------------------")
print("QQplot a Disperze - QUSP")
in_model = mod_qusp

qqnorm(resid(in_model)) ; qqline(resid(in_model))
plot(fitted(in_model), resid(in_model)) ; abline(h = 0, col = "red")