library(dplyr)
library(tibble)
library(ggplot2)
library(purrr)

## vypocetni funkce
prepare.cormat.data <- function(species_data, random_effects){
  
  # ## Testing arguments
  # species_data = fasy_mod_dataset
  # random_effects = ranef(mod_fasy)
  
  annual_means <- species_data %>%
    group_by(year) %>%
    summarise(
      sox = mean(sox, na.rm = TRUE),
      mean_temp = mean(mean_temp, na.rm = TRUE),
      mean_cwb = mean(mean_cwb, na.rm = TRUE),
      .groups = "drop"
    )
  
  names(random_effects) <- "RandomEffect"
  random_effects$year <- as.numeric(as.character(rownames(random_effects)))
  
  data_combined <- annual_means %>%
    left_join(random_effects, by = "year")
  
  output <- data_combined
  
  spp <- unique(species_data$species)
  
  if(length(spp) > 1){
    spp <- "All"
  }
  
  output$species <- spp
  
  return(output)
  
}

## dataset making
dataset <- rbind(prepare.cormat.data(main_mod_dataset, ranef(mod_main)),
                 prepare.cormat.data(abal_mod_dataset, ranef(mod_abal)),
                 prepare.cormat.data(pcab_mod_dataset, ranef(mod_pcab)),
                 prepare.cormat.data(pisy_mod_dataset, ranef(mod_pisy)),
                 prepare.cormat.data(fasy_mod_dataset, ranef(mod_fasy)),
                 prepare.cormat.data(qusp_mod_dataset, ranef(mod_qusp)))


## Graf
corr_data <- dataset %>%
  pivot_longer(
    cols = c(sox, mean_temp, mean_cwb),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(species, variable) %>%
  summarise(
    cor = cor(value, RandomEffect, use = "complete.obs", method = "pearson"),
    p_value = cor.test(value, RandomEffect, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    variable = recode(
      variable,
      sox = "SOx",
      mean_temp = "Temperature",
      mean_cwb = "CWB"
    ),
    label = ifelse(p_value < 0.05, sprintf("%.2f", cor), "")
  )


corr_data <- corr_data %>%
  mutate(
    species = factor(
      species,
      levels = c("All", "ABAL", "PCAB", "PISY", "FASY", "QUSP")
    ),
    variable = factor(
      variable,
      levels = c("SOx", "Temperature", "CWB")
    )
  )


ggplot(corr_data, aes(x = species, y = variable, fill = cor)) +
  geom_tile(color = "black", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson's r"
  ) +
  scale_y_discrete(
    limits = c("CWB", "Temperature", "SOx")
  ) +
  labs(x = "Species", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank()
  )
