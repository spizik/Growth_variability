## Functions ####
## Data preparation
# Fits GLS models with AR(1) errors for CV_RWI, temperature, and CWB trends at each site and classifies each trend by direction and significance.
# Returns a site-level summary table containing trend categories for variability, temperature, and water balance for the selected species.
prepare_data <- function(in.data){
  
  ## Testing arguments
  # input = clim.dataset
  
  output <- data.frame(site = character(),
                       species = character(),
                       cv_trend = character(),
                       temp_trend = character(),
                       cwb_trend = character())
  
  for(i in unique(in.data$site_code)){
    
    # i="P810114PCAB"
    # print(i)
    
    input <- subset(in.data, site_code == i)
    
    gls_rwi <- gls(cv_RWI ~ year, correlation = corAR1(form = ~ year), data = input)
    gls_temp <- gls(mean_temp ~ year, correlation = corAR1(form = ~ year), data = input)
    gls_cwb <- gls(cwb ~ year, correlation = corAR1(form = ~ year), data = input)
    
    gls_rwi <- summary(gls_rwi)[["tTable"]]
    gls_temp <- summary(gls_temp)[["tTable"]]
    gls_cwb <- summary(gls_cwb)[["tTable"]]
    
    out_gls_rwi <- NA
    out_temp_rwi <- NA
    out_cwb_rwi <- NA
    
    if(gls_rwi["year","p-value"]<0.05){
      if(gls_rwi["year","Value"]>0) out_gls_rwi<-"increasing_significant"
      if(gls_rwi["year","Value"]<0) out_gls_rwi<-"decreasing_significant"
    } else{
      if(gls_rwi["year","Value"]>0) out_gls_rwi<-"increasing_nonsignificant"
      if(gls_rwi["year","Value"]<0) out_gls_rwi<-"decreasing_nonsignificant"
    }
    if(gls_temp["year","p-value"]<0.05){
      if(gls_temp["year","Value"]>0) out_temp_rwi<-"increasing_significant"
      if(gls_temp["year","Value"]<0) out_temp_rwi<-"decreasing_significant"
    } else{
      if(gls_temp["year","Value"]>0) out_temp_rwi<-"increasing_nonsignificant"
      if(gls_temp["year","Value"]<0) out_temp_rwi<-"decreasing_nonsignificant"
    }
    if(gls_cwb["year","p-value"]<0.05){
      if(gls_cwb["year","Value"]>0) out_cwb_rwi<-"increasing_significant"
      if(gls_cwb["year","Value"]<0) out_cwb_rwi<-"decreasing_significant"
    } else{
      if(gls_cwb["year","Value"]>0) out_cwb_rwi<-"increasing_nonsignificant"
      if(gls_cwb["year","Value"]<0) out_cwb_rwi<-"decreasing_nonsignificant"
    }
    
    temp_output <- data.frame(site = i,
                              species = unique(input$species),
                              cv_trend = out_gls_rwi,
                              temp_trend = out_temp_rwi,
                              cwb_trend = out_cwb_rwi)
    
    output <- rbind(output, temp_output)
  }
  
  return(output)
}

# Summarizes the relative proportions of CV_RWI trend categories within each CWB trend class.
# Returns a normalized contingency table expressing how variability trends distribute across drought trend types.
sort_data <- function(input){
  
  ## Testing arguments
  # input = prepare_data(clim.dataset)
  
  output <- data.frame(matrix(NA,nrow=4,ncol=5))
  names(output) <- c("cwb_trend", "decreasing_nonsignificant", "decreasing_significant", "increasing_nonsignificant", "increasing_significant")
  output$cwb_trend <- c("decreasing_nonsignificant", "decreasing_significant", "increasing_nonsignificant", "increasing_significant")
  
  for(i in 1:nrow(output)){
    
    sub <- subset(input, cwb_trend == output$cwb_trend[i])
    tab <- table(sub$cv_trend)
    
    for (name in names(tab)) {
      output[i, name] <- tab[name]
    }
  }
  
  output[is.na(output)] <- 0
  num_cols <- setdiff(names(output), "cwb_trend")
  row_sums <- rowSums(output[, num_cols])
  output[, num_cols] <- output[, num_cols] / row_sums
  
  return(output)
}

# Converts the summarized contingency table into long format and assigns ordered x-axis positions for plotting.
# Returns a melted data frame where each bar segment is associated with a numeric x-axis coordinate.
prepare_graph_data <- function(input){
  
  ## Testing arguments
  # input <- sort_data(prepare_data(input))
  
  melted <- melt(input)
  melted$x_axis <- NA
  
  melted$x_axis[which(melted$cwb_trend == "decreasing_significant")] <- 1
  melted$x_axis[which(melted$cwb_trend == "decreasing_nonsignificant")] <- 2
  melted$x_axis[which(melted$cwb_trend == "increasing_nonsignificant")] <- 3
  melted$x_axis[which(melted$cwb_trend == "increasing_significant")] <- 4
  
  return(melted)
}

## Figure functions
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
plot_data <- function(input){
  
  ## Testing arguments
  # input = clim.dataset
  
  prep_dataset <- prepare_data(input)
  sorted_dta <- sort_data(prep_dataset)
  graph_data <- prepare_graph_data(sorted_dta)
  
  no_samples <- table(prep_dataset$cwb_trend)
  no_samples_df <- data.frame(
    cwb_trend = names(no_samples),
    no_samples = as.numeric(no_samples)
  )
  
  no_samples_df$x_axis <- NA
  no_samples_df$x_axis[which(no_samples_df$cwb_trend == "decreasing_significant")] <- 1
  no_samples_df$x_axis[which(no_samples_df$cwb_trend == "decreasing_nonsignificant")] <- 2
  no_samples_df$x_axis[which(no_samples_df$cwb_trend == "increasing_nonsignificant")] <- 3
  no_samples_df$x_axis[which(no_samples_df$cwb_trend == "increasing_significant")] <- 4
  
  graph_data$variable <- as.character(graph_data$variable)
  graph_data$variable[which(graph_data$variable == "decreasing_significant")] <- "1_decreasing_significant"
  graph_data$variable[which(graph_data$variable == "decreasing_nonsignificant")] <- "2_decreasing_nonsignificant"
  graph_data$variable[which(graph_data$variable == "increasing_nonsignificant")] <- "3_increasing_nonsignificant"
  graph_data$variable[which(graph_data$variable == "increasing_significant")] <- "4_increasing_significant"
  
  cols <- c("#3498DB", "#3498DB", "#E74C3C", "#E74C3C")
  names(cols) <- c("1_decreasing_significant", "2_decreasing_nonsignificant", "3_increasing_nonsignificant", "4_increasing_significant")
  
  alphas <- c(1.00, 0.50, 0.50, 1.00)
  names(alphas) <- c("1_decreasing_significant", "2_decreasing_nonsignificant", "3_increasing_nonsignificant", "4_increasing_significant")
  
  g <- ggplot(graph_data)
  g <- g + geom_bar(aes(x = x_axis, y = value, fill= variable, alpha=variable), stat = "identity", position = "stack")
  g <- g + geom_text(data = no_samples_df, aes(x = x_axis, y = 1.05, label = paste0("n=", no_samples)), inherit.aes = FALSE)
  
  g <- g + scale_fill_manual(values = cols, breaks = names(cols))
  g <- g + scale_alpha_manual(values = alphas, breaks = names(alphas))
  
  g <- g + scale_y_continuous("Within-site variance trend", limits = c(0, 1.1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20))
  g <- g + scale_x_continuous("Trend in CWB", limits = c(0.5, 4.5), breaks = c(1:4), labels = c("(--)", "(-)", "(+)", "(++)S"))
  g <- my.theme(g)
  g
}

##Figure making ####

figure <- ggarrange(plot_data(clim.dataset),
                    ggarrange(plot_data(subset(clim.dataset, species == "ABAL")),
                              plot_data(subset(clim.dataset, species == "PCAB")),
                              plot_data(subset(clim.dataset, species == "PISY")),
                              plot_data(subset(clim.dataset, species == "FASY")),
                              plot_data(subset(clim.dataset, species == "QUSP")),
                              nrow=3,ncol=2,labels=LETTERS[2:6],align="hv"),
                    nrow=2, ncol=1, labels=c("A",""),heights=c(0.25,0.75))





















 