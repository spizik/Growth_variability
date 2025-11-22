## Functions ####
boot.data.general<-function(input){
  
  ## Testing arguments
  # input=show_data_withinsite
  
  yrs<-unique(input$year)
  yrs<-yrs[order(yrs)]
  method<-unique(input$method)
  method<-method[order(method)]
  
  out<-data.frame(year=rep(yrs,times=length(method)),
                  method=rep(method,each=length(yrs)),
                  cv_min=NA,cv=NA,cv_max=NA)
  
  for(i in 1:nrow(out)){
    
    sub<-subset(input,year==out$year[i] & method==out$method[i])
    if(nrow(sub)>=5){
      out[i,c("cv_min","cv","cv_max")]<-quantile(apply(replicate(100,sample(sub$cv,nrow(sub),T)),2,mean),probs=c(0.025,0.50,0.975))
    }
  }
  
  return(out)
}
corr_cv_by_method <- function(input) {

  ## Testing arguments
  # input = show_data_betweensite
  
  df_wide <- input %>%
    select(year, method, cv) %>%
    pivot_wider(names_from = method, values_from = cv)
  
  cor_mat <- df_wide %>%
    select(-year) %>%
    cor(use = "pairwise.complete.obs")
  
  return(cor_mat)
}
plot.boxplot.variance.betweensite<-function(input, graph_name){
  
  ## Testing arguments
  # input=show_data_withinsite
  
  input$x<-NA
  input$x[which(input$method=="mean")]<-1
  input$x[which(input$method=="GAM")]<-2
  input$x[which(input$method=="qGAM")]<-3
  input$x[which(input$method=="spline")]<-4
  
  input<-na.omit(input)
  input$grp<-paste0(input$x,"_",input$species)
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=x,y=cv,fill=method, group=method), alpha=0.75)
  g<-g+scale_fill_manual(values=cols.methods, breaks=names(cols.methods))
  g<-g+scale_y_continuous(graph_name,limits=c(0,0.9),breaks=seq(0,1,0.2),labels=formatC(seq(0,1,0.2),format="f",digits=1))
  g<-g+scale_x_continuous("",limits=c(0.5,4.5),breaks=seq(1,4,1),labels=c("Mean", "GAM", "qGAM", "Spline"))
  g<-my.theme(g,"none")
  g
}
plot.data.betweensite<-function(input, graph_name){
  
  ## Testing arguments
  # input=show_data_betweensite
  
  g<-ggplot(input)
  g<-g+annotate("rect",xmin=1971,xmax=1992,ymin=0,ymax=0.89,fill="#BBBBBB",alpha=0.25)
  g<-g+geom_rect_pattern(
    xmin = 1992, xmax = 1998, ymin = 0, ymax = 0.89,
    inherit.aes = FALSE,
    pattern       = "stripe",   
    pattern_colour= "#BBBBBB",
    pattern_fill  = "#BBBBBB",
    pattern_angle = -45,
    pattern_density = .1,
    alpha=0.25,
    fill = NA,                 
    colour = NA
  )
  g<-g+geom_vline(xintercept=c(1992, 2003),linetype="dotted",colour="#D73027",linewidth=0.75,alpha=0.25)
  g<-g+geom_vline(xintercept=c(1996, 2010),linetype="dotted",colour="#26466D",linewidth=0.75,alpha=0.25)
  g<-g+geom_ribbon(aes(x=year,ymin=cv_min, ymax=cv_max,fill=method, group=method), alpha=0.2)
  g<-g+geom_point(aes(x=year,y=cv,fill=method,colour=method, group=method), alpha=0.75)
  g<-g+scale_colour_manual(values=cols.methods, breaks=names(cols.methods))
  g<-g+scale_fill_manual(values=cols.methods, breaks=names(cols.methods))
  g<-g+scale_y_continuous("Between-site variability",limits=c(0,0.9),breaks=seq(0,1,0.1),labels=formatC(seq(0,1,0.1),format="f",digits=1))
  # g<-g+scale_x_continuous("",limits=c(0.5,4.5),breaks=seq(1,4,1),labels=c("Mean", "GAM", "qGAM", "Spline"))
  g<-my.theme(g,"none")
  g
}
plot.data.withinsite<-function(input, graph_name){
  
  ## Testing arguments
  # input=show_data_withinsite
  
  input <- boot.data.general(input)
  
  g<-ggplot(input)
  g<-g+annotate("rect",xmin=1971,xmax=1992,ymin=0,ymax=0.89,fill="#BBBBBB",alpha=0.25)
  g<-g+geom_vline(xintercept=c(1992, 2003),linetype="dotted",colour="#D73027",linewidth=0.75,alpha=0.25)
  g<-g+geom_vline(xintercept=c(1996, 2010),linetype="dotted",colour="#26466D",linewidth=0.75,alpha=0.25)
  g<-g+geom_ribbon(aes(x=year,ymin=cv_min, ymax=cv_max,fill=method, group=method), alpha=0.2)
  g<-g+geom_point(aes(x=year,y=cv,fill=method,colour=method, group=method), alpha=0.75)
  g<-g+scale_colour_manual(values=cols.methods, breaks=names(cols.methods))
  g<-g+scale_fill_manual(values=cols.methods, breaks=names(cols.methods))
  g<-g+scale_y_continuous("Within-site variability",limits=c(0,0.9),breaks=seq(0,1,0.1),labels=formatC(seq(0,1,0.1),format="f",digits=1))
  # g<-g+scale_x_continuous("",limits=c(0.5,4.5),breaks=seq(1,4,1),labels=c("Mean", "GAM", "qGAM", "Spline"))
  g<-my.theme(g,"none")
  g
}

## Data loading ####
sp <- "PISY"

## Between-site
dta_mean <- subset(read.table("Calculated_datasets/Methods_compare_data/chronologies_data_mean.txt",sep=";"), species == sp)
dta_GAM <- subset(read.table("Calculated_datasets/Methods_compare_data/chronologies_data_GAM.txt",sep=";"), species == sp)
dta_qGAM <- subset(read.table("Calculated_datasets/Methods_compare_data/chronologies_data_qGAM.txt",sep=";"), species == sp)
dta_spline <- subset(read.table("Calculated_datasets/Methods_compare_data/chronologies_data_spline.txt",sep=";"), species == sp)

dta_mean$method <- "mean"
dta_GAM$method <- "GAM"
dta_qGAM$method <- "qGAM"
dta_spline$method <- "spline"

show_data_betweensite <- data.frame(year = c(dta_mean$year,
                                             dta_GAM$year,
                                             dta_qGAM$year,
                                             dta_spline$year),
                                    cv = c(dta_mean$cv.std.mid,
                                           dta_GAM$cv.std.mid,
                                           dta_qGAM$cv.std.mid,
                                           dta_spline$cv.std.mid),
                                    cv_min = c(dta_mean$cv.std.min,
                                               dta_GAM$cv.std.min,
                                               dta_qGAM$cv.std.min,
                                               dta_spline$cv.std.min),
                                    cv_max = c(dta_mean$cv.std.max,
                                               dta_GAM$cv.std.max,
                                               dta_qGAM$cv.std.max,
                                               dta_spline$cv.std.max),
                                    method = c(dta_mean$method,
                                               dta_GAM$method,
                                               dta_qGAM$method,
                                               dta_spline$method))

## Within-site
dta_mean <- subset(read.table("Calculated_datasets/Methods_compare_data/dataset_mean.txt",sep=";"), species == sp)
dta_GAM <- subset(read.table("Calculated_datasets/Methods_compare_data/dataset_GAM.txt",sep=";"), species == sp)
dta_qGAM <- subset(read.table("Calculated_datasets/Methods_compare_data/dataset_qGAM.txt",sep=";"), species == sp)
dta_spline <- subset(read.table("Calculated_datasets/Methods_compare_data/dataset_spline.txt",sep=";"), species == sp)

dta_mean$method <- "mean"
dta_GAM$method <- "GAM"
dta_qGAM$method <- "qGAM"
dta_spline$method <- "spline"

show_data_withinsite <- data.frame(site_code = c(dta_mean$site_code,
                                                 dta_GAM$site_code,
                                                 dta_qGAM$site_code,
                                                 dta_spline$site_code),
                                   year = c(dta_mean$year,
                                            dta_GAM$year,
                                            dta_qGAM$year,
                                            dta_spline$year),
                                   cv = c(dta_mean$cv_RWI,
                                          dta_GAM$cv_RWI,
                                          dta_qGAM$cv_RWI,
                                          dta_spline$cv_RWI),
                                   method = c(dta_mean$method,
                                              dta_GAM$method,
                                              dta_qGAM$method,
                                              dta_spline$method))




## Tests ####
print("Betweensite-korelace")
print(corr_cv_by_method(show_data_betweensite))


print("Withinsite-korelace")
print(corr_cv_by_method(boot.data.general(show_data_withinsite)))

## Figure ####
figure<-ggarrange(plot.boxplot.variance.betweensite(show_data_betweensite, "Between-site variability"),
                  plot.data.betweensite(show_data_betweensite),
                  plot.boxplot.variance.betweensite(show_data_withinsite, "Within-site variability"),
                  plot.data.withinsite(show_data_withinsite),
                  nrow=2,ncol=2,labels=LETTERS[1:4],widths=c(0.35,0.65))
