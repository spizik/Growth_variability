## Functions ####
## Calculations and data adjusting
# Standardizes heterogeneous sampling-design labels into a unified categorical variable for subsequent statistical comparison.
# Returns the input dataset with simplified sampling categories assigned based on original metadata codes.
assign_sampling_categories<-function(input){
  
  ## Testing arguments
  # input = working_aggregated_mean
  
  cat <- rep(NA, nrow(input))
  cat[which(input$data_collection_design %in% c("ALL", "all", "transect"))] <- "All"
  cat[which(input$data_collection_design %in% c("dom", "DOM"))] <- "Dom"
  cat[which(input$data_collection_design %in% c("ran", "RAN", "subjective_SEL/ALL", "grid"))] <- "Ran"
  cat[which(input$data_collection_design %in% c("windthrows"))] <- "Wind"
  cat[which(input$data_collection_design %in% c("death"))] <- "Dea"
  
  input$data_collection_design <- cat
  
  return(input)
}

# Tests whether intra-site variability differs across standardized sampling categories using Dunn’s post-hoc test.
# Prints adjusted significance levels to identify which sampling strategies produce significantly different variability outcomes.
evaluate_sampling_category <- function(input){
  
  ## Testing arguments
  # input = working_aggregated_mean
  
  input <- assign_sampling_categories(input)
  
  dunn.results<-dunnTest(cv_RWI~data_collection_design, data=input)$res
  pvals<-rep("non-significant", nrow(dunn.results))
  pvals[which(dunn.results$P.adj<0.05)]<-"significant"
  dunn.results$P.adj<-pvals
  
  print(dunn.results)
}

# Evaluates whether management types differ in intra-site variability by applying Dunn’s post-hoc test to management categories.
# Prints adjusted significance values highlighting significant pairwise contrasts among management groups.
evaluate_management_category <- function(input){
  
  ## Testing arguments
  # input = working_aggregated_mean
  
  dunn.results<-dunnTest(cv_RWI~management_type, data=input)$res
  pvals<-rep("non-significant", nrow(dunn.results))
  pvals[which(dunn.results$P.adj<0.05)]<-"significant"
  dunn.results$P.adj<-pvals
  
  print(dunn.results)
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
plot_sampling_strategy <- function(input){
  # input=working_aggregated_mean
  
  input <- assign_sampling_categories(input)
  
  sites_per_zone <- input %>% 
    distinct(site_code, data_collection_design) %>% 
    dplyr::count(data_collection_design) %>% 
    dplyr::rename(n_sites = n)
  sites_per_zone$y <- 0.95
  
  input$x <- NA ; sites_per_zone$x <- NA
  input$x[which(input$data_collection_design == "All")] <- 1 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "All")] <- 1
  input$x[which(input$data_collection_design == "Dea")] <- 2 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Dea")] <- 2
  input$x[which(input$data_collection_design == "Dom")] <- 3 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Dom")] <- 3
  input$x[which(input$data_collection_design == "Ran")] <- 4 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Ran")] <- 4
  # input$x[which(input$data_collection_design == "Sel")] <- 5 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Sel")] <- 5
  input$x[which(input$data_collection_design == "Wind")] <- 5 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Wind")] <- 5
  input$x[is.na(input$data_collection_design)] <- 6 ; sites_per_zone$x[is.na(sites_per_zone$data_collection_design)] <- 6
  
  
  
  g <- ggplot(input)
  g <- g + geom_boxplot(aes(x = x, y = cv_RWI, fill = data_collection_design), alpha = 0.75)
  g <- g + geom_label(data = sites_per_zone, aes(x = x, y = y, label = n_sites, colour=data_collection_design), fontface = "bold", size = 4, fill="#FFFFFF", alpha=0.40, label.size = 0 )
  g <- g + xlab("") + ylab("Within-site variability")
  g <- g + scale_fill_manual(values = sampling_colors, breaks = names(sampling_colors)) 
  g <- g + scale_colour_manual(values = sampling_colors, breaks = names(sampling_colors)) 
  g <- g + scale_x_continuous(limits = c(0.5, 6.5), breaks = c(1:6), labels = c("All", "Dea", "Dom", "Ran", "Wind", "NA"))
  g <- g + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = formatC(seq(0, 1, 0.25), format = "f", digits = 2))
  g <- my.theme(g)
  g
}
plot_sampling_strategy_sd <- function(input){
  # input=working_aggregated_mean
  
  input <- assign_sampling_categories(input)
  
  sites_per_zone <- input %>% 
    distinct(site_code, data_collection_design) %>% 
    dplyr::count(data_collection_design) %>% 
    dplyr::rename(n_sites = n)
  sites_per_zone$y <- 0.25
  
  input$x <- NA ; sites_per_zone$x <- NA
  input$x[which(input$data_collection_design == "All")] <- 1 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "All")] <- 1
  input$x[which(input$data_collection_design == "Dea")] <- 2 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Dea")] <- 2
  input$x[which(input$data_collection_design == "Dom")] <- 3 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Dom")] <- 3
  input$x[which(input$data_collection_design == "Ran")] <- 4 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Ran")] <- 4
  # input$x[which(input$data_collection_design == "Sel")] <- 5 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Sel")] <- 5
  input$x[which(input$data_collection_design == "Wind")] <- 5 ; sites_per_zone$x[which(sites_per_zone$data_collection_design == "Wind")] <- 5
  input$x[is.na(input$data_collection_design)] <- 6 ; sites_per_zone$x[is.na(sites_per_zone$data_collection_design)] <- 6
  
  
  
  g <- ggplot(input)
  g <- g + geom_boxplot(aes(x = x, y = cv_RWI, fill = data_collection_design), alpha = 0.75)
  g <- g + geom_label(data = sites_per_zone, aes(x = x, y = y, label = n_sites, colour=data_collection_design), fontface = "bold", size = 4, fill="#FFFFFF", alpha=0.40, label.size = 0 )
  g <- g + xlab("") + ylab("Within-site variability")
  g <- g + scale_fill_manual(values = sampling_colors, breaks = names(sampling_colors)) 
  g <- g + scale_colour_manual(values = sampling_colors, breaks = names(sampling_colors)) 
  g <- g + scale_x_continuous(limits = c(0.5, 6.5), breaks = c(1:6), labels = c("All", "Dea", "Dom", "Ran", "Wind", "NA"))
  g <- g + scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 1, 0.1), labels = formatC(seq(0, 1, 0.1), format = "f", digits = 2))
  g <- my.theme(g)
  g
}
plot_management <- function(input){
  # input=working_aggregated_mean
  
  sites_per_zone <- input %>% 
    distinct(site_code, management_type) %>% 
    dplyr::count(management_type) %>% 
    dplyr::rename(n_sites = n)
  sites_per_zone$y <- 0.95
  
  input$x <- NA ; sites_per_zone$x <- NA
  input$x[which(input$management_type == "Managed")] <- 1 ; sites_per_zone$x[which(sites_per_zone$management_type == "Managed")] <- 1
  input$x[which(input$management_type == "ManPa")] <- 2 ; sites_per_zone$x[which(sites_per_zone$management_type == "ManPa")] <- 2
  input$x[which(input$management_type == "Unmanaged")] <- 3 ; sites_per_zone$x[which(sites_per_zone$management_type == "Unmanaged")] <- 3
  input$x[is.na(input$management_type)] <- 4 ; sites_per_zone$x[is.na(sites_per_zone$management_type)] <- 4
  
  g <- ggplot(input)
  g <- g + geom_boxplot(aes(x = x, y = cv_RWI, fill = management_type), alpha = 0.75)
  g <- g + geom_label(data = sites_per_zone, aes(x = x, y = y, label = n_sites, colour=management_type), fontface = "bold", size = 4, fill="#FFFFFF", alpha=0.40, label.size = 0 )
  g <- g + xlab("") + ylab("Within-site variability")
  g <- g + scale_fill_manual(values = management_colors, breaks = names(management_colors)) 
  g <- g + scale_colour_manual(values = management_colors, breaks = names(management_colors)) 
  g <- g + scale_x_continuous(limits = c(0.5, 3.5), breaks = c(1:4), labels = c("Man", "ManPa", "UnMan", "NA"))
  g <- g + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = formatC(seq(0, 1, 0.25), format = "f", digits = 2))
  g <- my.theme(g)
  g
}
plot_management_sd <- function(input){
  # input=working_aggregated_sd
  
  sites_per_zone <- input %>% 
    distinct(site_code, management_type) %>% 
    dplyr::count(management_type) %>% 
    dplyr::rename(n_sites = n)
  sites_per_zone$y <- 0.25
  
  input$x <- NA ; sites_per_zone$x <- NA
  input$x[which(input$management_type == "Managed")] <- 1 ; sites_per_zone$x[which(sites_per_zone$management_type == "Managed")] <- 1
  input$x[which(input$management_type == "ManPa")] <- 2 ; sites_per_zone$x[which(sites_per_zone$management_type == "ManPa")] <- 2
  input$x[which(input$management_type == "Unmanaged")] <- 3 ; sites_per_zone$x[which(sites_per_zone$management_type == "Unmanaged")] <- 3
  input$x[is.na(input$management_type)] <- 4 ; sites_per_zone$x[is.na(sites_per_zone$management_type)] <- 4
  
  g <- ggplot(input)
  g <- g + geom_boxplot(aes(x = x, y = cv_RWI, fill = management_type), alpha = 0.75)
  g <- g + geom_label(data = sites_per_zone, aes(x = x, y = y, label = n_sites, colour=management_type), fontface = "bold", size = 4, fill="#FFFFFF", alpha=0.40, label.size = 0 )
  g <- g + xlab("") + ylab("Within-site variability")
  g <- g + scale_fill_manual(values = management_colors, breaks = names(management_colors)) 
  g <- g + scale_colour_manual(values = management_colors, breaks = names(management_colors)) 
  g <- g + scale_x_continuous(limits = c(0.5, 3.5), breaks = c(1:4), labels = c("Man", "ManPa", "UnMan", "NA"))
  g <- g + scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 1, 0.1), labels = formatC(seq(0, 1, 0.1), format = "f", digits = 2))
  g <- my.theme(g)
  g
}


## Calculations ####
working_df <- clim.dataset

working <- merge(working_df[,c("site_code","cv_RWI","species")],
                 site.list[, c("site_code", "data_collection_design", "management_type")],
                 by = "site_code",
                 all.x = TRUE)

working_aggregated_mean <- aggregate(cv_RWI ~ site_code + species + data_collection_design + management_type, data = working, FUN = mean)
working_aggregated_sd <- aggregate(cv_RWI ~ site_code + species + data_collection_design + management_type, data = working, FUN = sd)

sampling_colors <- c(
  All  = "#1B1B1B",   
  Dom  = "#2C7BB6",  
  Ran  = "#E6550D",  
  Wind = "#66C2A5", 
  Dea  = "#969696",  
  "NA" = "#E377C2"  
)


management_colors <- c(
  "Managed"    = "#006D2C",
  "ManPa"      = "#FCAE17", 
  "Unmanaged"  = "#984EA3", 
  "NA"         = "#CCCCCC"  
)

## Results ####
print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("--------------------- AOV a mean by sampling strategy ------------------")
print("------------------------------------------------------------------------")
evaluate_sampling_category(working_aggregated_mean)
print("--------------------- ÁBAL------------------")
evaluate_sampling_category(subset(working_aggregated_mean, species == "ABAL"))
print("--------------------- PCAB------------------")
evaluate_sampling_category(subset(working_aggregated_mean, species == "PCAB"))
evaluate_sampling_category(subset(working_aggregated_mean, species == "PISY"))
evaluate_sampling_category(subset(working_aggregated_mean, species == "FASY"))
evaluate_sampling_category(subset(working_aggregated_mean, species == "QUSP"))


print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("------------------------- AOV a mean by management ---------------------")
print("------------------------------------------------------------------------")
evaluate_management_category(working_aggregated_mean)
print("--------------------- ÁBAL------------------")
evaluate_management_category(subset(working_aggregated_mean, species == "ABAL"))
print("--------------------- PCAB------------------")
evaluate_management_category(subset(working_aggregated_mean, species == "PCAB"))
print("--------------------- PISY------------------")
evaluate_management_category(subset(working_aggregated_mean, species == "PISY"))
print("--------------------- FASY------------------")
evaluate_management_category(subset(working_aggregated_mean, species == "FASY"))
print("--------------------- QUSP------------------")
evaluate_management_category(subset(working_aggregated_mean, species == "QUSP"))

## Figure making ####
figure <- ggarrange(plot_sampling_strategy(working_aggregated_mean),
                    plot_management(working_aggregated_mean),
                    
                    plot_sampling_strategy(subset(working_aggregated_mean, species == "ABAL")),
                    plot_management(subset(working_aggregated_mean, species == "ABAL")),
                    
                    plot_sampling_strategy(subset(working_aggregated_mean, species == "PCAB")),
                    plot_management(subset(working_aggregated_mean, species == "PCAB")),
                    
                    plot_sampling_strategy(subset(working_aggregated_mean, species == "PISY")),
                    plot_management(subset(working_aggregated_mean, species == "PISY")),
                    
                    plot_sampling_strategy(subset(working_aggregated_mean, species == "FASY")),
                    plot_management(subset(working_aggregated_mean, species == "FASY")),
                    
                    plot_sampling_strategy(subset(working_aggregated_mean, species == "QUSP")),
                    plot_management(subset(working_aggregated_mean, species == "QUSP")),
                    
                    nrow=6, ncol=2, labels=LETTERS[1:12], 
                    widths = c(0.6,0.4), align = "hv", common.legend = T, legend = "bottom")