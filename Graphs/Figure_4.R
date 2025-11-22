## Funkce ####
## Data preparation and calculations
# Prepares a site-level chronology dataset by selecting single-species sites and restricting years to 1961–2017.
# Returns a long-format data frame with site, species, year, standardized and residual indices (std, res), and placeholder columns for coordinates.
prepare.crn.data<-function(input){
  
  ## Testing arguments
  # input=df.crn.cutted
  
  output.df <- data.frame(site = character(),
                          species = character(),
                          year = numeric(),
                          std = numeric(),
                          res = numeric(),
                          x = numeric(),
                          y = numeric())
  
  for(i in unique(input$site_code)){
    
    # i="V002460PCAB"
    # print(i)
    
    if(length(site.list$species[which(site.list$site_code == i)]) == 1){
      sub<-subset(input, site_code == i)
      sub<-sub[which(sub$year %in% c(1961:2017)),]
      
      temp <- data.frame(site = i,
                         species = site.list$species[which(site.list$site_code == i)],
                         year = sub$year,
                         std = sub$std,
                         res = sub$res
      )
      
      output.df <- rbind(output.df, temp)
    }
  }
  
  return(output.df)
}

# Prepares intra-site variability data by subsetting to the target period and selecting key variability and location variables.
# Returns a data frame with site code, species, year, RWI variability metrics (sd_RWI, cv_RWI), and site coordinates.
prepare.variability.data<-function(input){
  
  ## Testing arguments
  # input = clim.dataset
  
  input <- subset(input, year > 1960 & year < 2024)
  output.df <- input[,c("site_code", "species", "year", "sd_RWI", "cv_RWI", "x", "y")]
  
  return(output.df)
}

# Computes spatial variograms of standardized chronology values (std) in moving time windows across sites.
# Merges site metadata, transforms coordinates to a metric projection, and returns a combined data frame of semivariance versus distance for each window year.
calculate.variogram.std <- function(input, window_length=10, dist=5, cutoff=100){
  
  ## Testing arguments
  # input = subset(variogram_crn_data, species=="PCAB")
  # dist = 5
  # cutoff = 100
  # window_length = 10
  
  input<-merge(input, site.list[,c("site_code", "elevation", "lon", "lat")], by.x="site", by.y="site_code")
  
  # Initialize an empty list to store results
  spatial_variograms <- list()
  
  # Get unique years in dataset
  unique_years <- seq(min(input$year)+window_length/2-1, max(input$year)+window_length/2+1, window_length/2)
  
  # Loop through each year and compute the spatial variogram
  for (yr in unique_years) {
    
    # print(yr)
    
    sub_data <- subset(input, year>=(yr-window_length/2) & year<=(yr+window_length/2))
    
    if(length(unique(sub_data$site)) > 2){
      
      dataset <- aggregate(std ~ site, data = sub_data, FUN = mean)
      dataset$sd <- aggregate(std ~ site, data = sub_data, FUN = sd)$std
      dataset$elevation <- aggregate(elevation ~ site, data = sub_data, FUN = unique)$elevation
      dataset$lon <- aggregate(lon ~ site, data = sub_data, FUN = unique)$lon
      dataset$lat <- aggregate(lat ~ site, data = sub_data, FUN = unique)$lat
      
      dataset<-na.omit(dataset)
      
      ## Define WGS84 coordinate system (longitude, latitude)
      coords <- SpatialPoints(dataset[, c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      ## Transform to a metric projection (e.g., UTM Zone 33N for Central Europe)
      coords_utm <- spTransform(coords, CRS("+proj=utm +zone=33 +datum=WGS84 +units=km"))
      
      ## Add transformed coordinates back to the dataframe
      dataset$x_m <- coordinates(coords_utm)[,1]
      dataset$y_m <- coordinates(coords_utm)[,2]
      
      ## Ensure we have enough data points
      if(nrow(dataset) > 5) {
        
        ## Compute spatial variogram for this year
        variog_year <- variogram(std ~ elevation, locations = ~ x_m + y_m, data = dataset, width = dist, cutoff = cutoff)
        
        ## Convert to data frame and add year column
        variog_year_df <- as.data.frame(variog_year)
        variog_year_df$year <- yr  # Store the year
        
        ## Store in list
        spatial_variograms[[as.character(yr)]] <- variog_year_df
      }
    }
    
  }
  
  # Combine all variograms into a single data frame
  spatial_variograms_df <- bind_rows(spatial_variograms)
  
  return(spatial_variograms_df)
}

# Computes spatial variograms of intra-site variability (cv_RWI) in moving time windows across sites.
# Merges site metadata, transforms coordinates, and returns a combined data frame of semivariance versus distance for each window year.
calculate.variogram.cv <- function(input, window_length=10, dist=5, cutoff=100){
  
  ## Testing arguments
  # input = subset(variogram_site_data, species=="PCAB")
  # dist = 5
  # cutoff = 100
  # window_length = 10
  
  input<-merge(input, site.list[,c("site_code", "elevation", "lon", "lat")], by.x="site_code", by.y="site_code")
  
  #£ Initialize an empty list to store results
  spatial_variograms <- list()
  
  #£ Get unique years in dataset
  unique_years <- seq(min(input$year)+window_length/2-1, max(input$year)+window_length/2+1, window_length/2)
  
  # Loop through each year and compute the spatial variogram
  for (yr in unique_years) {
    
    # print(yr)
    
    sub_data <- subset(input, year>=(yr-window_length/2) & year<=(yr+window_length/2))
    
    if(length(unique(sub_data$site_code)) > 2){
      dataset <- aggregate(cv_RWI ~ site_code, data = sub_data, FUN = mean)
      dataset$sd_cv <- aggregate(cv_RWI ~ site_code, data = sub_data, FUN = sd)$cv_RWI
      dataset$elevation <- aggregate(elevation ~ site_code, data = sub_data, FUN = unique)$elevation
      dataset$lon <- aggregate(lon ~ site_code, data = sub_data, FUN = unique)$lon
      dataset$lat <- aggregate(lat ~ site_code, data = sub_data, FUN = unique)$lat
      
      dataset<-na.omit(dataset)
      
      ## Define WGS84 coordinate system (longitude, latitude)
      coords <- SpatialPoints(dataset[, c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      ## Transform to a metric projection (e.g., UTM Zone 33N for Central Europe)
      coords_utm <- spTransform(coords, CRS("+proj=utm +zone=33 +datum=WGS84 +units=km"))
      
      ## Add transformed coordinates back to the dataframe
      dataset$x_m <- coordinates(coords_utm)[,1]
      dataset$y_m <- coordinates(coords_utm)[,2]
      
      ## Ensure we have enough data points
      if(nrow(dataset) > 5) {
        
        ## Compute spatial variogram for this year
        variog_year <- variogram(cv_RWI ~ elevation, locations = ~ x_m + y_m, data = dataset, width = dist, cutoff = cutoff)
        
        ## Convert to data frame and add year column
        variog_year_df <- as.data.frame(variog_year)
        variog_year_df$year <- yr  # Store the year
        
        ## Store in list
        spatial_variograms[[as.character(yr)]] <- variog_year_df
      }
    }
  }
  
  # Combine all variograms into a single data frame
  spatial_variograms_df <- bind_rows(spatial_variograms)
  
  return(spatial_variograms_df)
}

## Tests

# Quantifies the effect of elevation on spatial variograms of standardized chronology values by comparing models with and without elevation as a trend.
# Returns a combined data frame with semivariances from both models and their difference (gamma_diff) for each distance class and window year.
calculate.variogram.std.eleveff <- function(input, window_length=10, dist=5, cutoff=100){
  
  ## Testing arguments
  # input = subset(variogram_crn_data, species=="PCAB")
  # dist = 5
  # cutoff = 100
  # window_length = 10
  
  input<-merge(input, site.list[,c("site_code", "elevation", "lon", "lat")], by.x="site", by.y="site_code")
  
  # Initialize an empty list to store results
  spatial_variograms <- list()
  
  # Get unique years in dataset
  # unique_years <- c(1961:2014)
  unique_years <- seq(min(input$year)+window_length/2-1, max(input$year)+window_length/2+1, window_length/2)
  
  # Loop through each year and compute the spatial variogram
  for (yr in unique_years) {
    
    # print(yr)
    
    sub_data <- subset(input, year>=(yr-window_length/2) & year<=(yr+window_length/2))
    
    if(length(unique(sub_data$site)) > 2){
      
      dataset <- aggregate(std ~ site, data = sub_data, FUN = mean)
      dataset$sd <- aggregate(std ~ site, data = sub_data, FUN = sd)$std
      dataset$elevation <- aggregate(elevation ~ site, data = sub_data, FUN = unique)$elevation
      dataset$lon <- aggregate(lon ~ site, data = sub_data, FUN = unique)$lon
      dataset$lat <- aggregate(lat ~ site, data = sub_data, FUN = unique)$lat
      
      dataset<-na.omit(dataset)
      
      # Define WGS84 coordinate system (longitude, latitude)
      coords <- SpatialPoints(dataset[, c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      # Transform to a metric projection (e.g., UTM Zone 33N for Central Europe)
      coords_utm <- spTransform(coords, CRS("+proj=utm +zone=33 +datum=WGS84 +units=km"))
      
      # Add transformed coordinates back to the dataframe
      dataset$x_m <- coordinates(coords_utm)[,1]
      dataset$y_m <- coordinates(coords_utm)[,2]
      
      # Ensure we have enough data points
      if(nrow(dataset) > 5) {
        
        # Compute spatial variogram for this year
        variog_elev <- variogram(std ~ elevation, locations = ~ x_m + y_m, data = dataset, width = dist, cutoff = cutoff)
        variog_noelev <- variogram(std ~ 1, locations = ~ x_m + y_m, data = dataset, width = dist, cutoff = cutoff)
        
        # Convert to data frame and add year column
        variog_elev <- as.data.frame(variog_elev)
        variog_noelev <- as.data.frame(variog_noelev)
        
        variog_year_df <- variog_elev[,c("np","dist","gamma")]
        names(variog_year_df) <- c("np","dist","gamma_elev")
        variog_year_df$gamma_noelev <- variog_noelev$gamma
        variog_year_df$gamma_diff <- variog_year_df$gamma_noelev - variog_year_df$gamma_elev
        variog_year_df$year <- yr  # Store the year
        
        # Store in list
        spatial_variograms[[as.character(yr)]] <- variog_year_df
      }
    }
    
  }
  
  # Combine all variograms into a single data frame
  spatial_variograms_df <- bind_rows(spatial_variograms)
  
  return(spatial_variograms_df)
}

# Quantifies the effect of elevation on spatial variograms of intra-site variability (cv_RWI) by comparing models with and without elevation as a trend.
# Returns a combined data frame with semivariances from both models and their difference (gamma_diff) for each distance class and window year.
calculate.variogram.cv.eleveff <- function(input, window_length=10, dist=5, cutoff=100){
  
  ## Testing arguments
  # input = subset(variogram_site_data, species=="PCAB")
  # dist = 5
  # cutoff = 100
  # window_length = 10
  
  input<-merge(input, site.list[,c("site_code", "elevation", "lon", "lat")], by.x="site_code", by.y="site_code")
  
  # Initialize an empty list to store results
  spatial_variograms <- list()
  
  # Get unique years in dataset
  unique_years <- seq(min(input$year)+window_length/2-1, max(input$year)+window_length/2+1, window_length/2)
  
  # Loop through each year and compute the spatial variogram
  for (yr in unique_years) {
    
    # print(yr)
    
    sub_data <- subset(input, year>=(yr-window_length/2) & year<=(yr+window_length/2))
    
    if(length(unique(sub_data$site_code)) > 2){
      dataset <- aggregate(cv_RWI ~ site_code, data = sub_data, FUN = mean)
      dataset$sd_cv <- aggregate(cv_RWI ~ site_code, data = sub_data, FUN = sd)$cv_RWI
      dataset$elevation <- aggregate(elevation ~ site_code, data = sub_data, FUN = unique)$elevation
      dataset$lon <- aggregate(lon ~ site_code, data = sub_data, FUN = unique)$lon
      dataset$lat <- aggregate(lat ~ site_code, data = sub_data, FUN = unique)$lat
      
      dataset<-na.omit(dataset)
      
      # Define WGS84 coordinate system (longitude, latitude)
      coords <- SpatialPoints(dataset[, c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      # Transform to a metric projection (e.g., UTM Zone 33N for Central Europe)
      coords_utm <- spTransform(coords, CRS("+proj=utm +zone=33 +datum=WGS84 +units=km"))
      
      # Add transformed coordinates back to the dataframe
      dataset$x_m <- coordinates(coords_utm)[,1]
      dataset$y_m <- coordinates(coords_utm)[,2]
      
      # Ensure we have enough data points
      if(nrow(dataset) > 5) {
        
        # Compute spatial variogram for this year
        variog_elev <- variogram(cv_RWI ~ elevation, locations = ~ x_m + y_m, data = dataset, width = dist, cutoff = cutoff)
        variog_noelev <- variogram(cv_RWI ~ 1, locations = ~ x_m + y_m, data = dataset, width = dist, cutoff = cutoff)
       
        # Convert to data frame and add year column
        variog_elev <- as.data.frame(variog_elev)
        variog_noelev <- as.data.frame(variog_noelev)
        
        variog_year_df <- variog_elev[,c("np","dist","gamma")]
        names(variog_year_df) <- c("np","dist","gamma_elev")
        variog_year_df$gamma_noelev <- variog_noelev$gamma
        variog_year_df$gamma_diff <- variog_year_df$gamma_noelev - variog_year_df$gamma_elev
        variog_year_df$year <- yr  # Store the year
        
        # Store in list
        spatial_variograms[[as.character(yr)]] <- variog_year_df
      }
    }
  }
  
  # Combine all variograms into a single data frame
  spatial_variograms_df <- bind_rows(spatial_variograms)
  
  return(spatial_variograms_df)
}

# Evaluates temporal trends in semivariance by fitting a GAMM of semivariance (gamma) against year with distance as a random effect.
# Prints the p-value for the year effect and R² values summarizing variance explained by year and by distance.
eval_variogram <- function(input){
  
  ## Testing arguments
  # input = calculate.variogram.std(subset(variogram_crn_data, species=="PCAB"))
  # input = calculate.variogram.cv(subset(variogram_site_data, species=="PCAB"))
  
  variogram_data <- input
  variogram_data$dist <- floor(variogram_data$dist)
  
  model <- gamm(gamma ~ s(year), random = list(dist = ~1), data = variogram_data)
  model_gam <- summary(model$gam)
  
  pval_year <- model_gam$s.table[4]
  rsqr <- r.squaredGLMM(model$lme)
  
  
  print(cat("pval pro year je: ", pval_year, "\n", 
            "rsqr pro year je: ", round(model_gam$r.sq,3), "\n", 
            "rsqr pro dist je: ", round(rsqr[2]-rsqr[1],3), "."))
}

# Assesses how often including elevation improves or worsens the variogram fit by counting positive and negative semivariance differences (gamma_diff).
# Prints the ratio of positive to negative elevation effects and the proportions of cases favoring the model with or without elevation.
eval_variogram_elev_effect <- function(input){
  
  ## Testing arguments
  # input = calculate.variogram.std.eleveff(subset(variogram_crn_data, species=="PCAB"))
  # input = calculate.variogram.cv.eleveff(subset(variogram_site_data, species=="PCAB"))
  
  n_pos <- sum(input$gamma_diff > 0)
  n_neg <- sum(input$gamma_diff < 0)
  
  perc_neg <- n_neg / (n_pos + n_neg)
  perc_pos <- n_pos / (n_pos + n_neg)
  ratio <- n_pos / n_neg
  
  
  print(cat("Ratio (pos_elev/neg_elev): ", round(ratio,2), "\n",
            "Percentage pos_elev: ", round(perc_pos,2), "\n",
            "Percentage neg_elev: ", round(perc_neg,2)))
}

## Graph functions
my.theme<-function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_line(colour="black"),
                     axis.text.x = element_text(colour="black"),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
prepare.tile.data <- function(input){
  matrix_data <- data.frame(dist = rep(seq(0, 100, 5), times = length(unique(input$year))),
                            years=rep(unique(input$year), each = length(seq(0, 100, 5))),
                            val=NA)
  
  for(i in unique(matrix_data$years)){
    new_df <- subset(matrix_data, years == i)
    sub_df <- subset(input, year == i)
    
    matrix_data$val[which(matrix_data$years == i)] <- as.numeric(predict(loess(gamma ~ dist, data=sub_df, span=0.75),newdata=new_df))
    
  }
  
  return(matrix_data)
}
plot.betweensite.variance.tile <- function(input){
  
  ## Testing arguments
  # input=prepare.tile.data(calculate.variogram.std(subset(variogram_crn_data, species=="PCAB"),10))
  
  paleta <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(10))
  
  g <- ggplot(input)
  g <- g + geom_tile(aes(x = years, y = dist, fill= val))
  g <- g + scale_fill_gradientn(colors = paleta,
                                limits = c(0, 0.2),
                                breaks = seq(0, 0.2, length.out = 11))
  g <- g + scale_y_continuous("Distance (km)", limits = c(0,100, breaks = seq(0, 200, 10)))
  g <- g + scale_x_continuous("Calendar year", limits = c(1960,2020, breaks = seq(0, 3000, 10)))
  g <- my.theme(g)
  g
}
plot.withinsite.variance.tile <- function(input){
  
  ## Testing arguments
  # input=prepare.tile.data(calculate.variogram.cv(subset(variogram_site_data, species=="PCAB"),10))
  
  paleta <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(10))
  
  g <- ggplot(input)
  g <- g + geom_tile(aes(x = years, y = dist, fill= val))
  g <- g + scale_fill_gradientn(colors = paleta,
                                limits = c(0, 0.1),
                                breaks = seq(0, 0.1, length.out = 11))
  g <- g + scale_y_continuous("Distance (km)", limits = c(0,100, breaks = seq(0, 200, 10)))
  g <- g + scale_x_continuous("Calendar year", limits = c(1960,2020, breaks = seq(0, 3000, 10)))
  g <- my.theme(g)
  g
}

## Calculations ####
df.crn.cutted <- df.crn.all[which(df.crn.all$site_code %in% clim.dataset$site_code),]
variogram_crn_data<-prepare.crn.data(df.crn.cutted)
variogram_site_data<-prepare.variability.data(clim.dataset)

variogram_crn_data<-subset(variogram_crn_data, year>1960 & year<2017)
variogram_site_data<-subset(variogram_site_data, year>1960 & year<2017)
variogram_crn_data <- variogram_crn_data[which(variogram_crn_data$site %in% variogram_site_data$site_code),]

## Testing ####
print("-----------------------------------------------------------------------")
print("-----------------------------------------------------------------------")
print("----------------------- Variogram evaluation --------------------------")
print("-----------------------------------------------------------------------")
print("--------------- std variance ABAL ---------------")
eval_variogram(calculate.variogram.std(subset(variogram_crn_data, species=="ABAL")))
print("--------------- std variance PCAB ---------------")
eval_variogram(calculate.variogram.std(subset(variogram_crn_data, species=="PCAB")))
print("--------------- std variance PISY ---------------")
eval_variogram(calculate.variogram.std(subset(variogram_crn_data, species=="PISY")))
print("--------------- std variance FASY ---------------")
eval_variogram(calculate.variogram.std(subset(variogram_crn_data, species=="FASY")))
print("--------------- std variance QUSP ---------------")
eval_variogram(calculate.variogram.std(subset(variogram_crn_data, species=="QUSP")))

print("--------------- withinsite variance ABAL ---------------")
eval_variogram(calculate.variogram.cv(subset(variogram_site_data, species=="ABAL")))
print("--------------- withinsite variance PCAB ---------------")
eval_variogram(calculate.variogram.cv(subset(variogram_site_data, species=="PCAB")))
print("--------------- withinsite variance PISY ---------------")
eval_variogram(calculate.variogram.cv(subset(variogram_site_data, species=="PISY")))
print("--------------- withinsite variance FASY ---------------")
eval_variogram(calculate.variogram.cv(subset(variogram_site_data, species=="FASY")))
print("--------------- withinsite variance QUSP ---------------")
eval_variogram(calculate.variogram.cv(subset(variogram_site_data, species=="QUSP")))

print("-----------------------------------------------------------------------")
print("-----------------------------------------------------------------------")
print("----------------- Evaluatng the effect of elevation -------------------")
print("-----------------------------------------------------------------------")
print("--------------- std elev effect ABAL ---------------")
eval_variogram_elev_effect(calculate.variogram.std.eleveff(subset(variogram_crn_data, species=="ABAL")))
print("--------------- std elev effect PCAB ---------------")
eval_variogram_elev_effect(calculate.variogram.std.eleveff(subset(variogram_crn_data, species=="PCAB")))
print("--------------- std elev effect PISY ---------------")
eval_variogram_elev_effect(calculate.variogram.std.eleveff(subset(variogram_crn_data, species=="PISY")))
print("--------------- std elev effect FASY ---------------")
eval_variogram_elev_effect(calculate.variogram.std.eleveff(subset(variogram_crn_data, species=="FASY")))
print("--------------- std elev effect QUSP ---------------")
eval_variogram_elev_effect(calculate.variogram.std.eleveff(subset(variogram_crn_data, species=="QUSP")))


print("--------------- withinsite elev effect ABAL ---------------")
eval_variogram_elev_effect(calculate.variogram.cv.eleveff(subset(variogram_site_data, species=="ABAL")))
print("--------------- withinsite elev effect PCAB ---------------")
eval_variogram_elev_effect(calculate.variogram.cv.eleveff(subset(variogram_site_data, species=="PCAB")))
print("--------------- withinsite elev effect PISY ---------------")
eval_variogram_elev_effect(calculate.variogram.cv.eleveff(subset(variogram_site_data, species=="PISY")))
print("--------------- withinsite elev effect FASY ---------------")
eval_variogram_elev_effect(calculate.variogram.cv.eleveff(subset(variogram_site_data, species=="FASY")))
print("--------------- withinsite elev effect QUSP ---------------")
eval_variogram_elev_effect(calculate.variogram.cv.eleveff(subset(variogram_site_data, species=="QUSP")))


## Figure ####
figure<-ggarrange(ggarrange(plot.betweensite.variance.tile(prepare.tile.data(calculate.variogram.std(subset(variogram_crn_data, species=="ABAL"),10))),
                            plot.betweensite.variance.tile(prepare.tile.data(calculate.variogram.std(subset(variogram_crn_data, species=="PCAB"),10))),
                            plot.betweensite.variance.tile(prepare.tile.data(calculate.variogram.std(subset(variogram_crn_data, species=="PISY"),10))),
                            plot.betweensite.variance.tile(prepare.tile.data(calculate.variogram.std(subset(variogram_crn_data, species=="FASY"),10))),
                            plot.betweensite.variance.tile(prepare.tile.data(calculate.variogram.std(subset(variogram_crn_data, species=="QUSP"),10))),
                            nrow = 5, ncol = 1, labels = c("A","C","E","G","I"), align = "hv", common.legend = T, legend = "bottom"),
                  
                  ggarrange(plot.withinsite.variance.tile(prepare.tile.data(calculate.variogram.cv(subset(variogram_site_data, species=="ABAL"),10))),
                            plot.withinsite.variance.tile(prepare.tile.data(calculate.variogram.cv(subset(variogram_site_data, species=="PCAB"),10))),
                            plot.withinsite.variance.tile(prepare.tile.data(calculate.variogram.cv(subset(variogram_site_data, species=="PISY"),10))),
                            plot.withinsite.variance.tile(prepare.tile.data(calculate.variogram.cv(subset(variogram_site_data, species=="FASY"),10))),
                            plot.withinsite.variance.tile(prepare.tile.data(calculate.variogram.cv(subset(variogram_site_data, species=="QUSP"),10))),
                            nrow = 5, ncol = 1, labels = c("B","D","F","H","J"), align = "hv", common.legend = T, legend = "bottom"),
                  
                  nrow = 1, ncol = 2, align = "hv")

