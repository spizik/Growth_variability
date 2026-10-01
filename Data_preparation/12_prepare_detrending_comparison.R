output_dir <- "Calculated_datasets/Methods_compare_data/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c(
  # mean = "mean", 
  # GAM = "GAM",
  qGAM = "qGAM", 
  spline = "Spline",
  spline50 = "Spline50"
)

for(method in names(methods)){
  
  print(method)
  
  ## RWI
  rwi.data <- list()
  
  for(i in names(prepared.data.tree.means)){
    print(i)
    rwi.data[[i]] <- create.rwi(
      prepared.data.tree.means[[i]],
      methods[[method]]
    )
  }
  
  ## Chronologie
  crn.data <- create.chronologies(rwi.data)
  
  ## Between-site: postup z 06_cv_chronologies.R
  dataset <- NULL
  
  for(i in names(crn.data)){
    temp <- crn.data[[i]]
    temp$site_code <- i
    temp$species <- substr(i, 8, 11)
    temp$year <- as.numeric(rownames(temp))
    temp$site_age <- seq_len(nrow(temp))
    dataset <- rbind(dataset, temp)
  }
  
  selected_sites <- subset(
    dataset, year == 1961 & samp.depth >= 5
  )$site_code
  
  dataset <- subset(
    dataset,
    site_code %in% selected_sites &
      site_age >= 20 & year >= 1961 & year <= 2020
  )
  dataset <- unify.categories(dataset)
  
  between <- NULL
  set.seed(1234)
  
  for(sp in c("ABAL", "PCAB", "PISY", "FASY", "QUSP")){
    between <- rbind(
      between,
      bootstrap.chronologies(dataset, sp)
    )
  }
  
  write.table(
    between,
    paste0(output_dir, "chronologies_data_", method, ".txt"),
    dec = ".", sep = ";"
  )
  
  ## Within-site: postup z 07_site_data_calculations.R
  raw_all <- calculate.all.sites.data(
    prepared.data.tree.means, site.list
  )
  
  within <- prepare.data(
    raw_all,
    min_age = 20,
    min_trees = 5,
    min_final_year = 2015,
    min_starting_year = 2015
  )
  within <- unify.categories(within)
  
  write.table(
    within,
    paste0(output_dir, "dataset_", method, ".txt"),
    dec = ".", sep = ";"
  )
}
