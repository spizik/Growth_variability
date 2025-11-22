## Functions ####
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
plot.mean.cv<-function(sp){
  
  ## Testing agruments
  # sp="PCAB"
  
  input <- subset(clim.dataset, species == sp)
  input <- aggregate(cv_RWI ~ site_code,data = input, FUN = mean)
  
  out <- merge(input, site.list[,c("site_code", "elevation")], by = "site_code")
  
  out$zone <- cut(x = out$elevation, 
                  breaks = elevation_breaks, 
                  labels = elevation_labels)
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  out$y <- elev_y[as.character(out$zone)]
  
  g<-ggplot(out)
  g <- g + geom_boxplot(aes(x=cv_RWI,y=y,fill=zone))
  g <- g + scale_x_continuous("CV mean", limits = c(0, 1), breaks=seq(0, 1, 0.2), labels=formatC(seq(0, 1, 0.2),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.range.cv<-function(sp){

  ## Testing agruments
  # sp="PCAB"
  
  input <- subset(clim.dataset, species == sp)
  input <- aggregate(cv_RWI ~ site_code,data = input, FUN = function(x){quantile(x, probs=0.975)-quantile(x, probs=0.025)})
  
  out <- merge(input, site.list[,c("site_code", "elevation")], by = "site_code")
  
  out$zone <- cut(x = out$elevation, 
                  breaks = elevation_breaks, 
                  labels = elevation_labels)
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  out$y <- elev_y[as.character(out$zone)]
  
  g<-ggplot(out)
  g <- g + geom_boxplot(aes(x=cv_RWI,y=y,fill=zone))
  g <- g + scale_x_continuous("CV range", limits = c(0, 1), breaks=seq(0, 1, 0.2), labels=formatC(seq(0, 1, 0.2),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.rbars<-function(sp){
  
  ## Testing agruments
  # sp="PCAB"
  
  sites <- unique(subset(clim.dataset, species == sp)$site_code)
  
  out <- data.frame(site_code = sites,
                    rBar=NA)
  
  for(i in 1:nrow(out)){
    rwis <- read.table(paste0("Calculated_datasets/rwi_data/",out$site_code[i],".txt"), sep=";")
    rwis <- rwis[which(as.numeric(rownames(rwis)) %in% c(1960:2017)),]
    rwis_t <- t(rwis)
    cor_mat <- cor(rwis_t, use = "pairwise.complete.obs")
    rBar <- mean(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
    out$rBar[i] <- rBar
  }
  
  out <- merge(out, site.list[,c("site_code", "elevation")], by = "site_code")
  
  out$zone <- cut(x = out$elevation, 
                  breaks = elevation_breaks, 
                  labels = elevation_labels)
  elev_y <- c(1:length(elevation_labels))
  names(elev_y) <- elevation_labels
  out$y <- elev_y[as.character(out$zone)]
  
  g<-ggplot(out)
  g <- g + geom_boxplot(aes(x=rBar,y=y,fill=zone))
  g <- g + scale_x_continuous("rBar", limits = c(0, 1), breaks=seq(0, 1, 0.2), labels=formatC(seq(0, 1, 0.2),format = "f", digits = 1))
  g <- g + scale_y_continuous("Elevation", limits = c(0.5, 6.5), breaks = elev_y, labels = names(elev_y))
  g <- g + scale_fill_manual(breaks = names(elev_colors), values = elev_colors)
  g <- my.theme(g)
  g
}
plot.panel<-function(sp){
  
  ## Testing arguments
  # sp="PCAB"
  
  out<-ggarrange(plot.rbars(sp),
                 plot.mean.cv(sp),
                 plot.range.cv(sp),
                 nrow=1,ncol=3,align="hv")
}

## Figure ####
figure <- ggarrange(plot.panel("ABAL"),
                    plot.panel("PCAB"),
                    plot.panel("PISY"),
                    plot.panel("FASY"),
                    plot.panel("QUSP"),
                    nrow=5,ncol=1,align="hv", labels=LETTERS[1:5])