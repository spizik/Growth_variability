
## Figure plotting ####
# plots maximum stand age of each analyzed site, with dunn test showing differences in stand ages among species
plot.stand.age <- function(input){
  ## testing arguments
  #input = main_mod_dataset
  
  input<- aggregate(median_age ~ site_code + species, data = input, FUN=max)
  input<-input[which(input$species %in% c("ABAL","PCAB","PISY","FASY","QUSP")),]
  input$x<-NA
  input$x[which(input$species=="ABAL")]<-1
  input$x[which(input$species=="PCAB")]<-2
  input$x[which(input$species=="PISY")]<-3
  input$x[which(input$species=="FASY")]<-4
  input$x[which(input$species=="QUSP")]<-5
  input<-na.omit(input)
  input$grp<-paste0(input$x,"_",input$species)
  
  aov <- dunnTest(median_age ~ grp, data = input, method = "bonferroni")
  aov$res$sigg <- ""
  aov$res$sigg[which(aov$res$P.adj<0.05)] <- "*"
  print(aov$res)
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=x,y=median_age,fill=species, group=grp), alpha=0.75, outliers = F)
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("Site age",limits=c(0,250),breaks=seq(0,500,50))
  g<-g+scale_x_continuous("",limits=c(0.5,5.5),breaks=seq(1,5,1),labels=c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
  g<-my.theme(g,"none")
  g
}

# plots maximum stand age heterogenity of each analyzed site, with dunn test showing differences in stand ages among species
plot.stand.age.heterogenity <- function(input){
  ## testing arguments
  #input = main_mod_dataset
  
  input<- aggregate(median_range ~ site_code + species, data = input, FUN=mean)
  input<-input[which(input$species %in% c("ABAL","PCAB","PISY","FASY","QUSP")),]
  input$x<-NA
  input$x[which(input$species=="ABAL")]<-1
  input$x[which(input$species=="PCAB")]<-2
  input$x[which(input$species=="PISY")]<-3
  input$x[which(input$species=="FASY")]<-4
  input$x[which(input$species=="QUSP")]<-5
  input<-na.omit(input)
  input$grp<-paste0(input$x,"_",input$species)
  
  aov <- dunnTest(median_range ~ grp, data = input, method = "bonferroni")
  aov$res$sigg <- ""
  aov$res$sigg[which(aov$res$P.adj<0.05)] <- "*"
  print(aov$res)
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=x,y=median_range,fill=species, group=grp), alpha=0.75, outliers = F)
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("Site age range",limits=c(0,250),breaks=seq(0,500,50))
  g<-g+scale_x_continuous("",limits=c(0.5,5.5),breaks=seq(1,5,1),labels=c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
  g<-my.theme(g,"none")
  g
}

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

## -----------------------------------------------------------------------------
## Figure output ####
figure<-ggarrange(plot.stand.age(main_mod_dataset),
                  plot.stand.age.heterogenity(main_mod_dataset),
                  nrow=2,ncol=1,labels=LETTERS[1:2])

