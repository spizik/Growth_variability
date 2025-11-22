## Functions ####
## Calculation functions
# Computes GLS-based temporal slopes (AR1) of standardized growth for each site of a given species.
# Returns a table with slope estimates and significance flags (“pos” for significant, “neg” for non-significant).
calculate.curve.slopes<-function(input){
  # input=subset(df.crn.all, species=="FASY")
  
  
  output<-data.frame(site_code=unique(input$site_code),
                     species=unique(input$species),
                     slope=NA,
                     pval=NA)
  
  
  for(i in 1:nrow(output)){
    sub_dataset <- subset(input, site_code == output$site_code[i])
    
    gls_model_trend <- gls(std ~ year, correlation = corAR1(form = ~ year), data = sub_dataset)
    
    output$slope[i] <- summary(gls_model_trend)$tTable[2,1]
    output$pval[i] <- summary(gls_model_trend)$tTable[1,4]
  }
  
  output<-na.omit(output)  
  
  pval <- rep("neg", nrow(output))
  pval[which(output$pval<0.05)] <- "pos"
  
  output$pval <- pval
  
  return(output)
}

## Testing functions
# Performs a Dunn post-hoc test to evaluate whether trend slopes differ among elevation zones within a species.
# Prints adjusted significance results for all pairwise comparisons between elevation categories.
testing.boxplots <- function(sp){
  
  # sp = "PCAB"
  
  input <- subset(curve.differences, species==sp)
  input$elev_zone <- cut(x = input$elevation, 
                         breaks = elevation_breaks, 
                         labels = elevation_labels)
  
  dunn.results<-dunnTest(slope ~ elev_zone, data = input)$res
  pvals<-rep("non-significant", nrow(dunn.results))
  pvals[which(dunn.results$P.adj<0.05)]<-"significant"
  dunn.results$P.adj<-pvals
  
  print(dunn.results)
}

# Summarizes the proportion of sites with positive, negative, or neutral trends within each elevation zone for a species.
# Returns a table reporting counts and relative frequencies of trend directions per elevation band.
summary.barplots <- function(sp){
  
  # sp = "ABAL"
  
  input <- subset(curve.differences, species==sp)
  input$elev_zone <- cut(x = input$elevation, 
                         breaks = elevation_breaks, 
                         labels = elevation_labels)
  
  input$value <- "negative"
  input$value[which(input$pval == "neg")] <- "none"
  input$value[which(input$slope >0 & input$pval =="pos")] <- "positive"
  
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

# Fits a linear regression to test whether site-level slopes systematically change with elevation for a given species.
# Prints the regression coefficient table showing trend–elevation relationships and their significance.
summary.slope <- function(sp){
  
  # sp = "PCAB"
  
  input <- subset(curve.differences, species==sp)
  
  reg <- lm(slope  ~ elevation, input)
  reg <- summary(reg)
  reg <- reg[["coefficients"]]
  
  print(reg)
}

## Plotting functions
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
  graph<-graph+theme(axis.line.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x=element_blank(),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks.y = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
plot.slopes<-function(input){
  
  # input=subset(curve.differences, species == "PCAB")
  
  # input <- calc.re(fasy_mod_dataset)
  
  reg <- lm(slope  ~ elevation, input)
  reg <- summary(reg)
  reg <- reg[["coefficients"]]
  
  sig <- ""
  if(reg[2, 4] < 0.05) sig <- "*"
  # txt<- paste0(round(reg[2, 1], 5), " ", sig)
  txt <- paste0(format(round(reg[2, 1],6), scientific = FALSE, digits = 6), " ", sig)
  
  
  label_data <- data.frame(
    elevation = 0,  # nebo jiná vhodná hodnota
    slope = 0.007,          # vertikální pozice labelu
    label = txt,
    species = unique(input$species)
  )
  
  g <- ggplot(input)
  g <- g + geom_point(aes(x=elevation, y=slope , colour=species, alpha=pval))
  g <- g + geom_hline(yintercept=0, linetype="dotted", linewidth=0.75, colour="#000000")
  g <- g + geom_smooth(aes(x=elevation, y=slope, colour=species, fill=species), method=lm, linewidth=0.75, linetype="solid", alpha=0.25)
  g <- g + geom_text(data=label_data, aes(x=elevation, y=slope, label=label, colour = species), inherit.aes=FALSE)
  g <- g + scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g <- g + scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g <- g + scale_y_continuous("Slope", limits=c(-0.016, 0.012), breaks=seq(-1,1,0.004), labels=formatC(seq(-1,1,0.004), format="f", digits=3))
  g <- g + scale_x_continuous("Elevation (m a.s.l.)" ,limits=c(-200, 1500), breaks=seq(0,2000,250))
  g <- my.theme(g, "none")
  g
}
plot.elev.boxplot <- function(input){
  
  show <- input
  show$elev_zone <- cut(x = show$elevation, 
                        breaks = elevation_breaks, 
                        labels = elevation_labels)
  
  show$x_graph <- NA
  for(i in 1:length(elevation_labels)){
    inds <- which(show$elev_zone == elevation_labels[i])
    if(length(inds) > 0) show$x_graph[inds] <- i
  }
  
  # sites_per_zone <- show %>%
  #   distinct(site_code, elev_zone) %>%  # vezmeme jen unikátní kombinace site a zóna
  #   count(elev_zone, name = "n_sites")  # spočítáme je podle elev_zone
  sites_per_zone <- aggregate(site_code ~ elev_zone, show, length)
  
  sites_per_zone$x<-NA
  for(i in 1:length(elevation_labels)){
    inds <- which(sites_per_zone$elev_zone == elevation_labels[i])
    if(length(inds) > 0) sites_per_zone$x[i] <- i
  }
  sites_per_zone$y <- 0.007
  
  sites_per_zone <- sites_per_zone[which(sites_per_zone$site_code>=5),]
  
  show<-na.omit(show)
  g<-ggplot(show)
  g<-g+geom_boxplot(aes(x=x_graph, y=slope, fill=elev_zone, group=elev_zone), alpha=0.75, outliers =F)
  
  g<-g+geom_label(data = sites_per_zone, aes(x = x, y = y, label = site_code, colour=elev_zone), fontface = "bold", size = 4, fill="#FFFFFF", alpha=0.40, label.size = 0 )
  
  g<-g+scale_colour_manual(values=elev_colors, breaks=names(elev_colors))
  g<-g+scale_fill_manual(values=elev_colors, breaks=names(elev_colors))
  g<-g+scale_y_continuous("STD slope",limits=c(-0.016,0.014),breaks=seq(-1,1,0.004), labels=formatC(seq(-1,1,0.004),format="f",digits=3))
  g<-g+scale_x_continuous("",limits=c(0.5,length(elevation_labels)+0.5),breaks=seq(1,length(elevation_labels),1),labels=elevation_labels)
  g<-my.theme.nox(g,"none")
  g
  
}
plot.elev.barplot <- function(input){
  
  show <- input
  show$elev_zone <- cut(x = show$elevation, 
                        breaks = elevation_breaks, 
                        labels = elevation_labels)
  
  show$value <- "negative"
  show$value[which(show$pval == "neg")] <- "none"
  show$value[which(show$slope >0 & show$pval =="pos")] <- "positive"

  # sites_per_zone <- input %>%
  #   distinct(site_code, elev_zone) %>%  # vezmeme jen unikátní kombinace site a zóna
  #   rename(elev_zone, name = "n_sites")  # spočítáme je podle elev_zone
  
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
  # g<-g+geom_barplot(aes(x=x_graph, y=perc, fill=val, group=elev_zone), alpha=0.75, outliers =F)
  g<-g+geom_bar(aes(x=x_graph, y=perc, group=elev_zone, fill=val),stat="identity", alpha=0.5)
  g<-g+scale_colour_manual(values=val_cols, breaks=names(val_cols))
  g<-g+scale_fill_manual(values=val_cols, breaks=names(val_cols))
  g<-g+scale_y_continuous("Value presence",limits=c(0,1),breaks=seq(0,1,0.2),labels=seq(0,100,20))
  g<-g+scale_x_continuous("",limits=c(0.5,length(elevation_labels)+0.5),breaks=seq(1,length(elevation_labels),1),labels=elevation_labels)
  g<-my.theme.nox(g,"none")
  g
  
}
plot.panel <- function(input){
  
  # input=subset(curve.differences, species == "QUSP")
  
  figure <- ggarrange(ggarrange(plot.elev.boxplot(input),
                                plot.elev.barplot(input),
                                nrow=2,ncol=1, align="hv"),
                      plot.slopes(input),
                      nrow=1, ncol=2)
  figure
}

## Calculations ####
df.crn.cutted <- df.crn.all[which(df.crn.all$site_code %in% clim.dataset$site_code),]

df.crn.cutted <- unify.categories(df.crn.cutted)

curve.differences<-rbind(calculate.curve.slopes(subset(df.crn.cutted, species=="ABAL")),
                         calculate.curve.slopes(subset(df.crn.cutted, species=="PCAB")),
                         calculate.curve.slopes(subset(df.crn.cutted, species=="PISY")),
                         calculate.curve.slopes(subset(df.crn.cutted, species=="FASY")),
                         calculate.curve.slopes(subset(df.crn.cutted, species=="QUSP")))

curve.differences <- merge(curve.differences, site.list[, c("site_code", "elevation")], by = "site_code")

curve.differences <- unify.categories(curve.differences)

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
figure<-ggarrange(plot.panel(subset(curve.differences, species == "ABAL")),
                  plot.panel(subset(curve.differences, species == "PCAB")),
                  plot.panel(subset(curve.differences, species == "PISY")),
                  plot.panel(subset(curve.differences, species == "FASY")),
                  plot.panel(subset(curve.differences, species == "QUSP")),
                  nrow=5, ncol=1, labels=LETTERS[1:5],align="hv")

