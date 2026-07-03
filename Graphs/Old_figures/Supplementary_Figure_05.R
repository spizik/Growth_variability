## Functions ####
## Calculations
# Performs bootstrapped estimation of inter-site variability (CV of standardized chronologies) across elevation zones for a given species.
# Returns a combined data frame of bootstrapped confidence intervals for each zone with associated species and zone labels.
boot.crns<-function(sp, n_boot=1000){
  
  ## Testing arguments
  # sp="PISY"
  # i=1
  
  zones<-unique(working_df$elev_zone)
  
  out<-data.frame(year=numeric(),
                  species=character(),
                  min=numeric(),mid=numeric(),max=numeric(),
                  zone=character())
  
  for(i in zones){
    # print(i)
    sub<-subset(working_df, species==sp & elev_zone == as.character(i))
    
    if(nrow(sub)>=5){
      booted.data<-boot.crn(sub, n_boot)
      booted.data$zone<-i
      out<-rbind(out,booted.data)
    }
  }
  
  return(out)
}

# Computes bootstrap confidence intervals for yearly CV of standardized chronologies within one elevation zone.
# Returns min, median, and max bootstrapped values for each year for the given species subset.
boot.crn<-function(input, n_boot=1000){
  
  ## Testing arguments
  # input=sub
  # n_boot=1000
  
  yrs<-unique(input$year)
  yrs<-yrs[order(yrs)]
  species<-unique(input$species)
  species<-species[order(species)]
  
  out<-data.frame(year=rep(yrs,times=length(species)),
                  species=rep(species,each=length(yrs)),
                  min=NA,mid=NA,max=NA)
  
  for(i in 1:nrow(out)){
    
    sub<-subset(input, year == out$year[i])
    
    if(nrow(sub)>=5){
      sam <- numeric(n_boot) 
      
      for (j in 1:n_boot) {
        sample_indices <- sample(1:nrow(sub), size = nrow(sub), replace = TRUE)
        sample_data <- sub[sample_indices, "std_chron"]
        
        cv <- sd(sample_data)/mean(sample_data)
        sam[j] <- cv
      }
      
      out[i,c("min","mid","max")]<-quantile(sam, c(0.025, 0.500, 0.975))
    }
  }
  
  return(out)
}

# Calculates synchrony (mean inter-series correlation r̄) between standardized chronologies for each elevation zone and species.
# Returns a data frame with r̄ values computed from pairwise correlations of series within each zone.
calc.synchro<-function(sp){
  
  ## Testing arguments
  # sp = "PCAB"
  # i=1
  
  zones<-unique(working_df$elev_zone)
  out<-data.frame(zone=zones, rbar=NA)
  
  for(i in zones){
    # print(i)
    input<-subset(working_df, species==sp & elev_zone == as.character(i))
    
    if(nrow(input)>=5){
      input <- as.data.frame(input)
      input <- input[, c("site_code", "year", "std_chron")]
      rownames(input) <- c(1:nrow(input))
      
      # Convert to wide format (year × series)
      rwi_matrix <- input %>%
        tidyr::pivot_wider(names_from = site_code, values_from = std_chron)
      
      rwi_matrix <- as.data.frame(rwi_matrix)
      rownames(rwi_matrix) <- rwi_matrix$year
      rwi_matrix <- rwi_matrix[, -1]
      
      cor_matrix <- cor(rwi_matrix, use = "pairwise.complete.obs")
      
      out$rbar[which(out$zone==i)] <- mean(cor_matrix[lower.tri(cor_matrix)], na.rm = TRUE)
    }
  }
  return(out)
}

# Performs bootstrapped estimation of mean TRW values across elevation zones for a given species.
# Returns min, median, and max bootstrapped envelopes for each year and elevation class.
boot.trws<-function(sp, n_boot=1000){
  
  ## Testing arguments
  # sp="PISY"
  # i=1
  
  zones<-unique(working_df$elev_zone)
  
  out<-data.frame(year=numeric(),
                  species=character(),
                  min=numeric(),mid=numeric(),max=numeric(),
                  zone=character())
  
  for(i in zones){
    # print(i)
    sub<-subset(working_df, species==sp & elev_zone == as.character(i))
    
    if(nrow(sub)>=5){
      booted.data<-boot.trw(sub, n_boot)
      booted.data$zone<-i
      out<-rbind(out,booted.data)
    }
  }
  
  return(out)
}

# Computes bootstrap envelopes for mean tree-ring width (TRW) across years for a given subset.
# Returns bootstrapped minimum, median, and maximum TRW values for each year.
boot.trw<-function(input, n_boot=1000){
  
  ## Testing arguments
  # input=sub
  # n_boot=1000
  
  yrs<-unique(input$year)
  yrs<-yrs[order(yrs)]
  species<-unique(input$species)
  species<-species[order(species)]
  
  out<-data.frame(year=rep(yrs,times=length(species)),
                  species=rep(species,each=length(yrs)),
                  min=NA,mid=NA,max=NA)
  
  for(i in 1:nrow(out)){
    
    sub<-subset(input, year == out$year[i])
    
    if(nrow(sub)>=5){
      sam <- numeric(n_boot) 
      
      for (j in 1:n_boot) {
        sample_indices <- sample(1:nrow(sub), size = nrow(sub), replace = TRUE)
        sample_data <- sub[sample_indices, "mid_TRW"]
        
        cv <- mean(sample_data)
        sam[j] <- cv
      }
      
      out[i,c("min","mid","max")]<-quantile(sam, c(0.025, 0.500, 0.975))
    }
  }
  
  return(out)
}

# Extracts unique soil chemistry and pH values per site and assigns elevation zones to each location.
# Returns a simplified soil dataset classified by elevation, ready for linking with TRW or variability metrics.
prepare.soil<-function(input){
  
  ## Testing arguments
  # input=pcab_mod_dataset
  
  unique_combinations <- input %>%
    distinct(site_code,elevation,  pH_L1, C.N_FH)
  
  out <- unique_combinations
  out$zone <- cut(x = out$elevation, 
                  breaks = elevation_breaks, 
                  labels = elevation_labels)
  
  return(out)
}


## Figures functions
my.theme<-function(graph,legend.pos="none"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_line(colour="black"),
                     axis.text.x = element_text(colour="black"),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
plot.crns.variability <- function(sp){
  input <- boot.crns(sp)
  
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  input$y <- elev_y[as.character(input$zone)]
  
  g <- ggplot(input)
  g <- g + geom_boxplot(aes(x = mid, y = y, fill = zone), alpha = 0.75)
  g <- g + scale_x_continuous("CRN variability", limits = c(0, 0.8), breaks=seq(0, 1, 0.2), labels=formatC(seq(0, 1, 0.2),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.trws.means <- function(sp){
  input <- boot.trws(sp)
  
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  input$y <- elev_y[as.character(input$zone)]
  
  g<-ggplot(input)
  g <- g + geom_boxplot(aes(x = mid, y = y, fill = zone), alpha = 0.75)
  g <- g + scale_x_continuous("mean TRW", limits = c(0, 2.5), breaks=seq(0, 5, 0.5), labels=formatC(seq(0, 5, 0.5),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.rbar <- function(sp){
  input <- calc.synchro(sp)
  
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  input$y <- elev_y[as.character(input$zone)]
  
  g<-ggplot(input)
  g <- g + geom_point(aes(x=rbar,y=y,colour=zone), size=3)
  g <- g + scale_x_continuous("rBar", limits = c(.15, 0.8), breaks=seq(0, 1, 0.1), labels=formatC(seq(0, 1, 0.1),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_colour_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.pH <- function(input){
  
  ## Testing arguments
  # input <- soil.dataset
  
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  input$y <- elev_y[as.character(input$zone)]
  
  g<-ggplot(input)
  g <- g + geom_boxplot(aes(x = pH, y = y, fill = zone), alpha = 0.75)
  g <- g + scale_x_continuous("pH", limits = c(3.75, 5.75), breaks=seq(0, 10, 0.5), labels=formatC(seq(0, 10, 0.5),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.c.n <- function(input){
  
  ## Testing arguments
  # input <- soil.dataset
  
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  input$y <- elev_y[as.character(input$zone)]
  
  g<-ggplot(input)
  g <- g + geom_boxplot(aes(x = C.N, y = y, fill = zone), alpha = 0.75)
  g <- g + scale_x_continuous("C:N", limits = c(10, 40), breaks=seq(0, 100, 10), labels=formatC(seq(0, 100, 10),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.data <- function(sp, soil.dataset){
  
  ## Testing arguments
  # sp="PISY"
  # soil.dataset=pisy_mod_dataset
  
  soil.dataset <- prepare.soil(soil.dataset)
  
  g<-ggarrange(plot.rbar(sp),
               plot.crns.variability(sp),
               plot.trws.means(sp),
               # plot.pH(soil.dataset),
               # plot.c.n(soil.dataset),
               nrow=1,ncol=3,widths=c(0.2, 0.4, 0.4))
  g
}

## Calcuations ####
working_df <- clim.dataset
working_df$elev_zone <- cut(x = working_df$elevation, 
                            breaks = elevation_breaks, 
                            labels = elevation_labels)

## Figure making ####
figure<-ggarrange(plot.data("ABAL", abal_mod_dataset),
                  plot.data("PCAB", pcab_mod_dataset),
                  plot.data("PISY", pisy_mod_dataset),
                  plot.data("FASY", fasy_mod_dataset),
                  plot.data("QUSP", qusp_mod_dataset),
                  nrow=5,ncol=1,labels=LETTERS[1:5])



  







