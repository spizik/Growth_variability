## Functions ####
## Figures functions
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
plot.sd.relation<-function(input, barva="black"){
  # input=mod_dataset
  # barva="black"
  
  input$log_sample_depth <- base::log(input$sample_depth)
  
  g<-ggplot(input)
  g<-g+geom_point(aes(x=cv_RWI, y=log_sample_depth,colour=species), alpha=0.10, shape=16)
  g<-g+geom_hline(yintercept = base::log(20), linetype="dotted", linewidth=0.8, colour="red")
  g<-g+geom_smooth(aes(x=cv_RWI, y=log_sample_depth),colour="black", fill="black", alpha=0.40)
  
  g<-g+geom_boxplot(aes(x = cv_RWI, y = .75), width = 1.5, fill = barva, color = "black", outliers = F, alpha=0.20)
  g<-g+scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("log Sample depth",
                          limits=c(0,base::log(300)),
                          breaks=c(0,base::log(10), base::log(25), base::log(50), base::log(100), base::log(150), base::log(200), base::log(300)),
                          labels=c(0, 10, 25, 50, 100, 150, 200, 300))
  g<-g+scale_x_continuous("within-site variability", limits=c(0, 1.5), breaks=seq(0, 1.5, 0.25), labels=formatC(seq(0, 1.5, 0.25),format="f",digits = 2))
  g<-my.theme(g,"none")
  g
}

## Figure making ####
figure <- ggarrange(plot.sd.relation(clim.dataset),
                    ggarrange(plot.sd.relation(subset(clim.dataset, species== "ABAL"), cols.species["ABAL"]),
                              plot.sd.relation(subset(clim.dataset, species== "PCAB"), cols.species["PCAB"]),
                              plot.sd.relation(subset(clim.dataset, species== "PISY"), cols.species["PISY"]),
                              plot.sd.relation(subset(clim.dataset, species== "FASY"), cols.species["FASY"]),
                              plot.sd.relation(subset(clim.dataset, species== "QUSP"), cols.species["QUSP"]),
                              ncol=2,nrow=3, align="hv",labels=LETTERS[2:6]),
                    nrow=2,ncol=1,heights=c(0.25,0.75),labels=c("A",""))
