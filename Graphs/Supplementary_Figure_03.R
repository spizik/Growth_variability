## Functions ####
## Data preparation
# Extracts site-level random effects from an LME model of intra-site variability while accounting for sample depth and AR1 temporal correlation.
# Returns a data frame with random effects, species identity, and geographic attributes merged from site metadata.
calc.re<-function(mod_dataset){
  
  ## Testing arguments
  # mod_dataset = subset(clim.dataset, species=="PCAB")
  
  in_sites_lme<-subset(mod_dataset,year>1960 & year<2018)[,c("year","cv_RWI","sample_depth","site_code")]
  names(in_sites_lme)<-c("year","variance","sample_depth","site_code")
  in_sites_lme$type<-"intra"
  
  mod <- lme(fixed = variance ~ year + sample_depth, random = ~1 | site_code, correlation = corAR1(form = ~ year), data = in_sites_lme)
  
  output <- data.frame(site_code = rownames(ranef(mod)),
                       species=unique(mod_dataset$species),
                       re = ranef(mod)[,1],
                       elevation = NA,
                       x = NA,
                       y = NA)
  for(i in 1:nrow(output)){
    # print(i)
    ind <- which(output$site_code[i] == site.list$site_code)
    if(length(ind)>0){
      output[i,"elevation"] <- site.list[ind, "elevation"]
      output[i,"x"] <- site.list[ind, "lon"]
      output[i,"y"] <- site.list[ind, "lat"]
    }
  }
  return(na.omit(output))
}

## Data analysis
# Performs a Dunn post-hoc test to evaluate whether random effects differ among elevation zones for a given species.
# Prints adjusted significance levels for all pairwise comparisons between elevation categories.
testing.boxplots <- function(sp){
  
  ## Testing arguments
  # sp="PCAB"
  
  input <- subset(clim.dataset, species==sp)
  input <- calc.re(input)
  
  input$elev_zone <- cut(x = input$elevation, 
                         breaks = elevation_breaks, 
                         labels = elevation_labels)
  
  dunn.results<-dunnTest(re ~ elev_zone, data = input)$res
  pvals<-rep("non-significant", nrow(dunn.results))
  pvals[which(dunn.results$P.adj<0.05)]<-"significant"
  dunn.results$P.adj<-pvals
  
  print(dunn.results)
}

# Summarizes the proportion of sites with positive or negative random effects within each elevation zone for a species.
# Returns a table reporting total sites, counts of positive/negative effects, and their relative frequencies by elevation class.
summary.barplots <- function(sp){
  
  ## Testing arguments
  # input=subset(clim.dataset, species=="PCAB")
  
  input <- subset(clim.dataset, species==sp)
  input <- calc.re(input)
  
  input$elev_zone <- cut(x = input$elevation, 
                         breaks = elevation_breaks, 
                         labels = elevation_labels)
  
  input$value <- "negative"
  input$value[which(input$re>0)] <- "positive"
  
  sites_per_zone <- aggregate(site_code ~ elev_zone,
                              data = unique(input[c("site_code", "elev_zone")]),
                              FUN = length)
  names(sites_per_zone)[2] <- "n_sites"
  
  input <- aggregate(site_code ~ elev_zone + value, input, length)
  
  sites_per_zone <- sites_per_zone[which(sites_per_zone$n_sites>=5),]
  
  input$perc <- NA
  for(i in 1:nrow(sites_per_zone)){
    
    inds <- which(input$elev_zone==sites_per_zone$elev_zone[i])
    input$perc[inds] <- input$site_code[inds] / sites_per_zone$n_sites[i]
  }
  
  names(input) <- c("elev_zone", "direction", "NO_sites", "rel_sites")
  
  print(input)
}

# Tests whether site-level random effects vary systematically with elevation for a given species using a simple linear model.
# Prints the regression coefficient table showing the magnitude and significance of the elevation effect.
summary.slope <- function(sp){
  
  ## Testing arguments
  # input=subset(clim.dataset, species=="PCAB")
  
  input <- subset(clim.dataset, species==sp)
  input <- calc.re(input)
  
  reg <- lm(re  ~ elevation, input)
  reg <- summary(reg)
  reg <- reg[["coefficients"]]
  
  print(reg)
}


## Graphs functions
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
my.theme.nox<-function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_line(colour="black"),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x=element_blank(),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks.y = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
my.theme.nox.b<-function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_line(colour="black"),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_line(color = "black"),
                     axis.title.x=element_blank(),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks.y = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
plot.re<-function(input){
  
  ## Testing arguments
  # dataset=abal_mod_dataset
  
  reg <- lm(re ~ elevation, input)
  reg <- summary(reg)
  reg <- reg[["coefficients"]]
  
  sig <- ""
  if(reg[2, 4] < 0.05) sig <- "*"
  txt <- paste0(format(round(reg[2, 1],6), scientific = FALSE, digits = 6), " ", sig)

  
  label_data <- data.frame(
    elevation = 0,
    re = 0.4, 
    label = txt,
    species = unique(input$species)
  )
  
  g <- ggplot(input)
  g <- g + geom_point(aes(x=elevation, y=re, colour=species))
  g <- g + geom_hline(yintercept=0, linetype="dotted", linewidth=0.75, colour="#000000")
  g <- g + geom_smooth(aes(x=elevation, y=re, colour=species, fill=species), method=lm, linewidth=0.75, linetype="solid", alpha=0.25)
  g <- g + geom_text(data=label_data, aes(x=elevation, y=re, label=label, colour = species), inherit.aes=FALSE)
  g <- g + scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g <- g + scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g <- g + scale_y_continuous("Random effect", limits=c(-0.5, 0.5), breaks=seq(-1,1,0.2), labels=formatC(seq(-1,1,0.2), format="f", digits=1))
  g <- g + scale_x_continuous("Elevation (m a.s.l.)" ,limits=c(-200, 1500), breaks=seq(0,2000,250))
  g <- my.theme(g, "none")
  g
}
plot.elev.boxplot <- function(input){
  
  ## Testing arguments
  # input=subset(clim.dataset, species=="PCAB")
  
  show <- input
  show$elev_zone <- cut(x = show$elevation, 
                        breaks = elevation_breaks, 
                        labels = elevation_labels)
  
  show$x_graph <- NA
  for(i in 1:length(elevation_labels)){
    inds <- which(show$elev_zone == elevation_labels[i])
    if(length(inds) > 0) show$x_graph[inds] <- i
  }
  
  sites_per_zone <- aggregate(site_code ~ elev_zone,
                              data = unique(show[c("site_code", "elev_zone")]),
                              FUN = length)
  names(sites_per_zone)[2] <- "n_sites"
  
  sites_per_zone$x<-NA
  for(i in 1:length(elevation_labels)){
    inds <- which(sites_per_zone$elev_zone == elevation_labels[i])
    if(length(inds) > 0) sites_per_zone$x[i] <- i
  }
  sites_per_zone$y <- 0.25
  
  sites_per_zone <- sites_per_zone[which(sites_per_zone$n_sites>=5),]
  
  g<-ggplot(show)
  g<-g+geom_boxplot(aes(x=x_graph, y=re, fill=elev_zone, group=elev_zone), alpha=0.75, outliers =F)
  
  g<-g+geom_label(data = sites_per_zone, aes(x = x, y = y, label = n_sites, colour=elev_zone), fontface = "bold", size = 4, fill="#FFFFFF", alpha=0.40, label.size = 0 )
  
  g<-g+scale_colour_manual(values=elev_colors, breaks=names(elev_colors))
  g<-g+scale_fill_manual(values=elev_colors, breaks=names(elev_colors))
  g<-g+scale_y_continuous("RE values",limits=c(-0.25,0.25),breaks=seq(-1,1,0.1),labels=formatC(seq(-1,1,0.1),format="f",digits=1))
  g<-g+scale_x_continuous("",limits=c(0.5,length(elevation_labels)+0.5),breaks=seq(1,length(elevation_labels),1),labels=elevation_labels)
  g<-my.theme.nox(g,"none")
  g
  
}
plot.elev.barplot <- function(input){
  
  ## Testing arguments
  # input=subset(clim.dataset, species=="PCAB")
  
  show <- input
  show$elev_zone <- cut(x = show$elevation, 
                        breaks = elevation_breaks, 
                        labels = elevation_labels)
  
  show$value <- "negative"
  show$value[which(show$re>0)] <- "positive"
  
  sites_per_zone <- aggregate(site_code ~ elev_zone,
                              data = unique(show[c("site_code", "elev_zone")]),
                              FUN = length)
  names(sites_per_zone)[2] <- "n_sites"
  
  show <- aggregate(site_code ~ elev_zone + value, show, length)
  
  sites_per_zone <- sites_per_zone[which(sites_per_zone$n_sites>=5),]
  
  show$perc <- NA
  for(i in 1:nrow(sites_per_zone)){
    
    inds <- which(show$elev_zone==sites_per_zone$elev_zone[i])
    show$perc[inds] <- show$site_code[inds] / sites_per_zone$n_sites[i]
    
  }
  show <- na.omit(show)
  names(show) <- c("elev_zone", "val", "no_sites", "perc" )
  
  show$x_graph <- NA
  for(i in 1:length(elevation_labels)){
    inds <- which(show$elev_zone == elevation_labels[i])
    if(length(inds) > 0) show$x_graph[inds] <- i
  }
  
  
  g<-ggplot(show)
  g<-g+geom_bar(aes(x=x_graph, y=perc, group=elev_zone, fill=val),stat="identity", alpha=0.5)
  g<-g+scale_colour_manual(values=val_cols, breaks=names(val_cols))
  g<-g+scale_fill_manual(values=val_cols, breaks=names(val_cols))
  g<-g+scale_y_continuous("RE presence",limits=c(0,1),breaks=seq(0,1,0.2),labels=seq(0,100,20))
  g<-g+scale_x_continuous("",limits=c(0.5,length(elevation_labels)+0.5),breaks=seq(1,length(elevation_labels),1),labels=elevation_labels)
  g<-my.theme.nox(g,"none")
  g
  
}
plot.panel <- function(input){
  
  ## Testing arguments
  # input=subset(clim.dataset, species=="PCAB")
  
  input <- calc.re(input)
  
  figure <- ggarrange(ggarrange(plot.elev.boxplot(input),
                                plot.elev.barplot(input),
                                nrow=2,ncol=1, align="hv"),
                      plot.re(input),
                      nrow=1, ncol=2)
  figure
}

## Tests ####
print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("---------------------------------- ABAL --------------------------------")
sp <- "ABAL"
testing.boxplots(sp)
summary.barplots(sp)
summary.slope(sp)

print("------------------------------------------------------------------------")
print("---------------------------------- PCAB --------------------------------")
sp <- "PCAB"
testing.boxplots(sp)
summary.barplots(sp)
summary.slope(sp)

print("------------------------------------------------------------------------")
print("---------------------------------- PISY --------------------------------")
sp <- "PISY"
testing.boxplots(sp)
summary.barplots(sp)
summary.slope(sp)

print("------------------------------------------------------------------------")
print("---------------------------------- FASY --------------------------------")
sp <- "FASY"
testing.boxplots(sp)
summary.barplots(sp)
summary.slope(sp)

print("------------------------------------------------------------------------")
print("---------------------------------- QUSP --------------------------------")
sp <- "QUSP"
testing.boxplots(sp)
summary.barplots(sp)
summary.slope(sp)

## Figure making ####
figure<-ggarrange(plot.panel(subset(clim.dataset, species=="ABAL")),
                  plot.panel(subset(clim.dataset, species=="PCAB")),
                  plot.panel(subset(clim.dataset, species=="PISY")),
                  plot.panel(subset(clim.dataset, species=="FASY")),
                  plot.panel(subset(clim.dataset, species=="QUSP")),
                  nrow=5, ncol=1, labels=LETTERS[1:5],align="hv")






