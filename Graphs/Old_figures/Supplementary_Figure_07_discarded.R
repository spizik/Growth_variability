## Functions ####
## Calculations
# Selects the most age-heterogeneous sites for a species and splits their trees into young and old cohorts based on percentile thresholds.
# Calculates per-year CV of RWI for young and old trees and tests cohort differences using Kruskal–Wallis, paired t-tests, and GLS with AR(1) errors.
calculate_dataset<-function(sp, perc_sites=0.25, perc.trees=0.2, min_trees_per_site=20){
  
  ## Testing arguments
  # sp="PCAB"
  # perc_sites=0.25
  # perc.trees=0.2
  # min_trees_per_site=20
  
  sites<-subset(clim.dataset, species==sp)
  
  ranges_biggest<-aggregate(cambial_age_range~site_code, data=sites, mean)
  notrees_sites <-aggregate(sample_depth~site_code, data=sites, mean)
  
  ranges_biggest <- merge(ranges_biggest, notrees_sites, by = "site_code")
  ranges_biggest <- subset(ranges_biggest, sample_depth>=min_trees_per_site)
  
  ranges_biggest<-ranges_biggest[order(ranges_biggest$cambial_age_range, decreasing = T),]
  ranges_biggest<-ranges_biggest[c(1:(round(nrow(ranges_biggest)*perc_sites))),]
  
  sub_dataset<-prepared.data.tree.means[which(names(prepared.data.tree.means) %in% ranges_biggest$site_code)]
  sub_dataset_young<-list()
  sub_dataset_old<-list()
  
  for(i in names(sub_dataset)){
    sub<-sub_dataset[[i]]
    tree_ages<-aggregate(cambial.age~TreeID, data=sub, max)
    tresholds<-quantile(tree_ages$cambial.age,probs=c(perc.trees,(1-perc.trees)))
    
    young_trees<-subset(tree_ages,cambial.age<=tresholds[1])
    old_trees<-subset(tree_ages,cambial.age>=tresholds[2])
    
    sub_dataset_young[[i]]<-sub[which(sub$TreeID %in% young_trees$TreeID),]
    sub_dataset_old[[i]]<-sub[which(sub$TreeID %in% old_trees$TreeID),]
  }
  
  output.crn.young<-data.frame(cv_value=numeric(), site_code=numeric(), species=numeric(), year=character())
  output.crn.old<-data.frame(cv_value=numeric(), site_code=numeric(), species=numeric(), year=character())
  
  for(i in names(sub_dataset_young)){
    rwi<-create.rwi(sub_dataset_young[[i]])
    temp<-data.frame(cv_value=apply(rwi,1,sd,na.rm=T)/apply(rwi,1,mean,na.rm=T),
                     site_code=i,
                     species=substr(i,8,11),
                     year=as.numeric(as.character(rownames(rwi))))
    
    output.crn.young<-rbind(output.crn.young,temp)
  }
  
  for(i in names(sub_dataset_old)){
    rwi<-create.rwi(sub_dataset_old[[i]])
    temp<-data.frame(cv_value=apply(rwi,1,sd)/apply(rwi,1,mean),
                     site_code=i,
                     species=substr(i,8,11),
                     year=as.numeric(as.character(rownames(rwi))))
    
    output.crn.old<-rbind(output.crn.old,temp)
  }
  
  output.crn.young<-subset(output.crn.young, year>=1961 & year<=2017)
  output.crn.old<-subset(output.crn.old, year>=1961 & year<=2017)
  output.crn.young$class<-"young"
  output.crn.old$class<-"old"
  
  output<-data.frame(site_code=unique(output.crn.young$site_code),
                     aov=NA, ttest=NA, gls_res=NA)
  
  for(i in output$site_code){
    # print(i)
    sub_dta<-rbind(subset(output.crn.young, site_code==i),
                   subset(output.crn.old, site_code==i))
    
    if(nrow(na.omit(sub_dta))>0){
      
      av<-kruskal.test(cv_value~class, sub_dta)
      output[which(output$site_code==i),"aov"]<-av$p.value
      
      df_ttest<-data.frame(year=c(1961:2017),
                           young=NA, old=NA)
      df_ttest$young[which(df_ttest$year %in% subset(sub_dta,class=="young")$year)]<-subset(sub_dta,class=="young")$cv_value
      df_ttest$old[which(df_ttest$year %in% subset(sub_dta,class=="old")$year)]<-subset(sub_dta,class=="old")$cv_value
      df_ttest<-na.omit(df_ttest)
      
      tt<-t.test(df_ttest$young, df_ttest$old, paired = TRUE)
      output[which(output$site_code==i),"ttest"]<-round(tt$p.value,3)
      
      age<-"old"
      if(median(df_ttest$young)>median(df_ttest$old)) age<-"young"
      
      output[which(output$site_code==i),"abs"]<-age
      
      diff.gls<-NA
      sub_dta<-na.omit(sub_dta)
      common.years<-table(sub_dta$year)
      common.years<-common.years[which(common.years==2)]
      common.years<-as.numeric(as.character(names(common.years)))
      sub_dta<-sub_dta[which(sub_dta$year %in% common.years),]
      
      tryCatch({
        gls_model <- gls(cv_value ~ class, correlation = corAR1(form = ~ year | class), data = sub_dta)
        gls_model <- summary(gls_model)
        if(gls_model$tTable["classyoung","p-value"]<0.05){
          if(gls_model$coefficients["classyoung"]>0) diff.gls<-"young"
          if(gls_model$coefficients["classyoung"]<0) diff.gls<-"old"
        } else{
          diff.gls<-"none"
        } }, error = function(e) {
          diff.gls<-NA
        })
      
      output[which(output$site_code==i),"gls_res"]<-diff.gls
    }
  }
  
  output<-na.omit(output)
  
  output$abs[which(output$aov>0.05)]<-"none"
  
  output$abs_val<-0
  output$abs_val[which(output$abs=="old")]<-1
  output$abs_val[which(output$abs=="young")]<--1
  
  output$abs_val_gls<-0
  output$abs_val_gls[which(output$gls_res=="old")]<-1
  output$abs_val_gls[which(output$gls_res=="young")]<--1
  
  ranges_biggest<-ranges_biggest[which(ranges_biggest$site_code %in% output$site_code),]
  ranges_biggest<-ranges_biggest[order(ranges_biggest$site_code),]
  output<-output[order(output$site_code),]
  output$range<-ranges_biggest$cambial_age_range
  
  sub_sites<-site.list[which(site.list$site_code %in% output$site_code),]
  sub_sites<-sub_sites[order(sub_sites$site_code),]
  output$elevation<-sub_sites$elevation
  
  return(output)
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
boxplot.aov.elev<-function(input){
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=abs_val, y=elevation, fill=abs), alpha=0.75)
  g<-g+scale_fill_manual(values=cols.class, breaks=names(cols.class))
  g<-g+scale_x_continuous("",limits=c(-1.5,1.5), breaks=c(-1,0,1), labels=c("Young", "No diff.", "Old"))
  g<-g+scale_y_continuous("Elevation",limits=c(150,1603), breaks=seq(0,2000,250))
  g<-my.theme(g)
  g
  
}

boxplot.trend.elev<-function(input){
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=abs_val_gls, y=elevation, fill=gls_res), alpha=0.75)
  g<-g+scale_fill_manual(values=cols.class, breaks=names(cols.class))
  g<-g+scale_x_continuous("",limits=c(-1.5,1.5), breaks=c(-1,0,1), labels=c("Young", "No diff.", "Old"))
  g<-g+scale_y_continuous("Elevation",limits=c(150,1603), breaks=seq(0,2000,250))
  g<-my.theme(g)
  g
  
}
barplot.differences<-function(input){
  
  show.abs<-data.frame(x=1,
                       y=as.numeric(table(input$abs)/nrow(input)),
                       fill=names(table(input$abs)))
  
  show.trend<-data.frame(x=2,
                         y=as.numeric(table(input$gls_res)/nrow(input)),
                         fill=names(table(input$gls_res)))
  
  
  show<-rbind(show.abs,show.trend)
  
  show$barorder<-0
  show$barorder[which(show$fill=="young")]<--1
  show$barorder[which(show$fill=="old")]<-1
  
  g<-ggplot(show)
  g<-g+geom_bar(aes(x=x, y=y, group=barorder, fill=fill), stat = "identity", position = "stack", alpha=0.75)
  g<-g+scale_fill_manual(values=cols.class, breaks=names(cols.class))
  g<-g+scale_x_continuous("",limits=c(0.5,2.5), breaks=c(1,2), labels=c("Value", "Trend"))
  g<-g+scale_y_continuous("",limits=c(0,1), breaks=seq(0,1,0.2), labels=seq(0,100,20))
  g<-my.theme(g)
  g
  
}
plot.panel<-function(input){
  
  ## Testing arguments
  # input = calculate_dataset("PCAB")
  
  g<-ggarrange(barplot.differences(input),
               boxplot.aov.elev(input),
               boxplot.trend.elev(input),
               ncol=3, nrow=1, align="hv", common.legend=T, legend="bottom")
  g
}

## Figure making ####
figure<-ggarrange(plot.panel(calculate_dataset("ABAL")),
                  plot.panel(calculate_dataset("PCAB")),
                  plot.panel(calculate_dataset("PISY")),
                  plot.panel(calculate_dataset("FASY")),
                  plot.panel(calculate_dataset("QUSP")),
                  ncol=1, nrow=5, align="hv", common.legend=T, legend="none",labels=LETTERS[1:5])
