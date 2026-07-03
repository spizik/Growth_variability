## Functions ####
## Data calculations

# Computes bootstrap confidence intervals (2.5%, 50%, 97.5%) of TRW and RWI coefficients of variation for each species–year combination. 
# Returns a data frame with these intervals for all species–year pairs with at least five observations. 
boot.data.general<-function(input){
  
  ## Testing arguments
  # input=clim.dataset
  
  yrs<-unique(input$year)
  yrs<-yrs[order(yrs)]
  species<-unique(input$species)
  species<-species[order(species)]
  
  out<-data.frame(year=rep(yrs,times=length(species)),
                  species=rep(species,each=length(yrs)),
                  trw_min=NA,trw_mid=NA,trw_max=NA,
                  rwi_min=NA,rwi_mid=NA,rwi_max=NA)
  
  for(i in 1:nrow(out)){
    
    sub<-subset(input,year==out$year[i] & species==out$species[i])
    if(nrow(sub)>=5){
      out[i,c("trw_min","trw_mid","trw_max")]<-quantile(apply(replicate(100,sample(sub$cv_TRW,nrow(sub),T)),2,mean),probs=c(0.025,0.50,0.975))
      out[i,c("rwi_min","rwi_mid","rwi_max")]<-quantile(apply(replicate(100,sample(sub$cv_RWI,nrow(sub),T)),2,mean),probs=c(0.025,0.50,0.975))
    }
  }
  
  return(out)
}

## Data analysis
# For each site, compares intra-site variability (cv_RWI) with inter-site variability (cv.std.mid) and computes their correlation using GLS with AR(1). 
# Returns a data frame with site-specific within-site trends, differences relative to between-site trends, and the correlation category between the two curves. 
calculate.curve.differences<-function(in.crn, in.main){
  
  ## Testing arguments
  # in.crn=subset(df.crn,species=="PCAB")
  # in.main=pcab_mod_dataset
  
  
  output<-data.frame(site_code=character(),
                     trend_within=character(),
                     diff_toBetween=character(),
                     corr_withBetween=character())

  
  for(i in unique(in.main$site_code)){

    # i="C004015PCAB"
    # print(i)
    
    sub.site<-subset(in.main, site_code==i)
    sub.site<-sub.site[order(sub.site$year),]
    sub.crn<-subset(in.crn,year>=min(sub.site$year) & year<=max(sub.site$year))
    
    if(nrow(sub.site)>30){
      eval.df<-data.frame(year=sub.site$year,
                          cv.site=sub.site$cv_RWI,
                          cv.crn=NA)
      
      eval.df$cv.crn<-in.crn$cv.std.mid[which(eval.df$year %in% in.crn$year)]
      eval.df.2<-melt(eval.df, id=c("year"))
      
      # Corelations
      corr.pval<-cor.test(eval.df$cv.site,eval.df$cv.crn)$p.value
      cor.est<-cor.test(eval.df$cv.site, eval.df$cv.crn)[["estimate"]][["cor"]]
      
      if(corr.pval<0.05){
        if(cor.est>0) corr<-"positive"
        if(cor.est<0) corr<-"negative"
      } else{
        corr<-"none"
      }
      
      # Differences in the mean values
      data <- data.frame(
        value = c(eval.df$cv.site, 
                  eval.df$cv.crn),
        time = c(eval.df$year, 
                 eval.df$year),
        group = c(rep("intra", nrow(eval.df)),
                  rep("inter", nrow(eval.df)))
      )
      
      # Intra-site variability trends
      trend.intra<-NA
      # gls_model_trend <- gls(value ~ time * group, correlation = corAR1(form = ~ time | group), data = data)
      gls_model_trend <- gls(value ~ time, correlation = corAR1(form = ~ time), data = subset(data,group=="intra"))
      gls_model <- summary(gls_model_trend)
      if(gls_model$tTable["time","p-value"]<0.05){
        if(gls_model$coefficients["time"]>0) trend.intra<-"increasing"
        if(gls_model$coefficients["time"]<0) trend.intra<-"decreasing"
      } else{
        trend.intra<-"none"
      }
      
      # Differences between intra and inter-site variability
      diff.trend<-NA
      gls_model_trend <- gls(value ~ time * group, correlation = corAR1(form = ~ time | group), data = data)
      gls_model <- summary(gls_model_trend)
      if(gls_model$tTable["time:groupintra","p-value"]<0.05){
        if(gls_model$coefficients["time:groupintra"]>0) diff.trend<-"higher_thanBetween"
        if(gls_model$coefficients["time:groupintra"]<0) diff.trend<-"lower_thanBetween"
      } else{
        diff.trend<-"same"
      }
      
      # output
      temp.output<-data.frame(site_code=i,
                              trend_within=trend.intra,
                              diff_toBetween=diff.trend,
                              corr_withBetween=corr)
      
      output<-rbind(output,temp.output)
    }
  }
  
  output<-na.omit(output)  
  output$species<-unique(in.crn$species)
  
  return(output)
}

# For each site, fits a GLS model with AR(1) of standardized growth (std) against year to estimate the temporal slope. 
# Returns a table of site codes, species, slope estimates, and a simple significance flag (“pos” for p < 0.05, “neg” otherwise). 
calculate.curve.slopes<-function(input){
  
  ## Testing arguments
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

# For a given species, fits GLS/LME models to intra-site (bootstrapped and LME-based) and inter-site variability time series to assess temporal trends. 
# Prints detailed model summaries to the console to compare intra-site, inter-site, and combined variability trends. 
sumarize.trends.results<-function(sp){
  
  ## Testing arguments
  # sp="ABAL"
  # mod_dataset = pcab_mod_dataset
  
  mod_dataset<-subset(clim.dataset, species==sp)
  
  in_sites<-subset(booted.data, species==sp & year>1960 & year<2018)[,c("year","rwi_mid")]
  names(in_sites)<-c("year","variance")
  in_sites$type<-"intra"
  
  in_sites_lme<-subset(mod_dataset,year>1960 & year<2018)[,c("year","cv_RWI","sample_depth","site_code")]
  names(in_sites_lme)<-c("year","variance","sample_depth","site_code")
  in_sites_lme$type<-"intra"
  
  in_crn<-subset(df.crn, species==sp & year>1960 & year<2018)[,c("year","cv.std.mid")]
  names(in_crn)<-c("year","variance")
  in_crn$type<-"inter"
  
  combination<-rbind(in_sites, in_crn)
  
  mod_intra<-gls(variance ~ year, correlation = corAR1(form = ~ year), data = in_sites)
  mod_intra_lme <- lme(fixed = variance ~ year + sample_depth, random = ~1 | site_code, correlation = corAR1(form = ~ year), data = in_sites_lme)
  mod_inter<-gls(variance ~ year, correlation = corAR1(form = ~ year), data = in_crn)
  mod_both <- gls(variance ~ year * type, correlation = corAR1(form = ~ year | type), data = combination)
  print("----------------------------------------------------")
  print(paste0("-------------------------",sp,"-----------------------"))
  print("----------------------------------------------------")
  print("----------- Intra-site variability - LME -----------")
  print(summary(mod_intra_lme)[["tTable"]])
  print("----------------------------------------------------")
  print("------- Intra-site variability - booted data -------")
  print(summary(mod_intra)[["tTable"]])
  print("-------------- Inter-site variability --------------")
  print(summary(mod_inter)[["tTable"]])
  print("-------------- Between variability --------------")
  print(summary(mod_both)[["tTable"]])
}

# For a given species, summarizes how many sites show increasing, decreasing, or no within-site variability trend based on curve.differences. 
# Prints the proportion of sites in each within-site trend category. 
sumarize.trends.directionss<-function(sp){
  
  ## Testing arguments
  # sp="ABAL"
  
  input <- subset(curve.differences,species==sp)
  
  dta<-aggregate(site_code~trend_within, data=input, FUN=length)
  dta$site_code<-dta$site_code/nrow(subset(curve.differences,species==sp))
  print(dta)
}

# For a given species, summarizes how intra-site trends differ from between-site trends (higher, lower, or similar) using curve.differences. 
# Prints the proportion of sites in each category describing agreement between intra-site and inter-site variability trends. 
sumarize.trend.agreement<-function(sp){
  
  ## Testing arguments
  # sp="ABAL"
  
  dta<-aggregate(site_code~species+diff_toBetween,data=subset(curve.differences,species==sp), FUN=length)
  dta$site_code<-dta$site_code/nrow(subset(curve.differences,species==sp))
  print(dta)
}

# For a given species, summarizes the correlation direction (positive, negative, none) between intra-site and inter-site variability curves from curve.differences. 
# Prints the proportion of sites in each correlation category. 
sumarize.correlationss<-function(sp){
  
  ## Testing arguments
  # sp="ABAL"
  
  dta<-aggregate(site_code~species+corr_withBetween,data=subset(curve.differences,species==sp), FUN=length)
  dta$site_code<-dta$site_code/nrow(subset(curve.differences,species==sp))
  print(dta)
}

## Figure plotting
create.boxplot.dataset<-function(input, input.crn){
  
  ## Testing arguments
  # input=subset(clim.dataset,species=="FASY")
  # input.crn=subset(df.crn,species=="FASY")
  
  grp.periods<-input[,c("year","species","cv_RWI")]
  grp.periods$group<-"one"
  grp.periods$group[which(grp.periods$year>1990)]<-"two"
  grp.periods$x<-1-0.2
  grp.periods$x[which(grp.periods$group=="two")]<-2-0.2
  grp.periods$group<-paste0("sites_",grp.periods$group)
  names(grp.periods)<-c("year", "species", "cv", "group", "x")
  
  std.grp.periods<-input.crn[,c("year","species","cv.std.mid")]
  std.grp.periods$group<-"one"
  std.grp.periods$group[which(std.grp.periods$year>1990)]<-"two"
  std.grp.periods$x<-1+0.2
  std.grp.periods$x[which(std.grp.periods$group=="two")]<-2+0.2
  std.grp.periods$group<-paste0("crn_",std.grp.periods$group)
  std.grp.periods$species<-"none"
  names(std.grp.periods)<-c("year", "species", "cv", "group", "x")
  
  out<-rbind(grp.periods, std.grp.periods)
  
  return(out)
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
plot.boxplot.variance.withinsite<-function(input,input.crn){
  
  ## Testing arguments
  # input=clim.dataset
  
  input<-input[which(input$species %in% c("ABAL","PCAB","PISY","FASY","QUSP")),]
  
  input$x<-NA
  input$x[which(input$species=="ABAL")]<-1
  input$x[which(input$species=="PCAB")]<-2
  input$x[which(input$species=="PISY")]<-3
  input$x[which(input$species=="FASY")]<-4
  input$x[which(input$species=="QUSP")]<-5
  
  input<-input[,c("x","cv_RWI","species")]
  names(input)<-c("x","cv","species")
  
  input<-na.omit(input)
  input$grp<-paste0(input$x,"_",input$species)
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=x,y=cv,fill=species, group=grp), alpha=0.75)
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("Within-site variability",limits=c(0,1),breaks=seq(0,1,0.2),labels=formatC(seq(0,1,0.2),format="f",digits=1))
  g<-g+scale_x_continuous("",limits=c(0.5,5.5),breaks=seq(1,5,1),labels=c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
  g<-my.theme(g,"none")
  g
}
plot.boxplot.variance.betweensite<-function(input,input.crn){
  
  ## Testing arguments
  # input=df.crn
  
  input<-input[which(input$species %in% c("ABAL","PCAB","PISY","FASY","QUSP")),]
  
  input$x<-NA
  input$x[which(input$species=="ABAL")]<-1
  input$x[which(input$species=="PCAB")]<-2
  input$x[which(input$species=="PISY")]<-3
  input$x[which(input$species=="FASY")]<-4
  input$x[which(input$species=="QUSP")]<-5
  
  input<-input[,c("x","cv.std.mid","species")]
  names(input)<-c("x","cv","species")
  
  input<-na.omit(input)
  input$grp<-paste0(input$x,"_",input$species)
  
  g<-ggplot(input)
  g<-g+geom_boxplot(aes(x=x,y=cv,fill=species, group=grp), alpha=0.75)
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("Between-site variability",limits=c(0,0.5),breaks=seq(0,1,0.1),labels=formatC(seq(0,1,0.1),format="f",digits=1))
  g<-g+scale_x_continuous("",limits=c(0.5,5.5),breaks=seq(1,5,1),labels=c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
  g<-my.theme(g,"none")
  g
}
plot.curve.differences<-function(input){
  
  # input=subset(curve.differences, species=="ABAL")
  
  sub_crn_slopes <- subset(curve.slopes, species == unique(input$species))
  
  show<-data.frame(x=c(1,1,1, # 
                       2,2,2, # trend v intra-site variablilite +,-, žádný
                       3,3,3 # korelace True-False
                       ), 
                   y=c(
                     # Rrendy v chronologiích
                     length(which(sub_crn_slopes$pval=="pos"))/nrow(sub_crn_slopes),
                     length(which(sub_crn_slopes$pval=="none"))/nrow(sub_crn_slopes),
                     length(which(sub_crn_slopes$pval=="neg"))/nrow(sub_crn_slopes),
                     
                     # Trend v inter-site
                     length(which(input$trend_within=="increasing"))/nrow(input),
                     length(which(input$trend_within=="none"))/nrow(input),
                     length(which(input$trend_within=="decreasing"))/nrow(input),
                     
                     # korelace
                     length(which(input$corr_withBetween =="positive"))/nrow(input),
                     length(which(input$corr_withBetween =="none"))/nrow(input),
                     length(which(input$corr_withBetween =="negative"))/nrow(input) #,
                     
                     
                     # rozdíl v trendu
                     # length(which(input$diff_toBetween=="higher_thanBetween"))/nrow(input),
                     # length(which(input$diff_toBetween=="same"))/nrow(input),
                     # length(which(input$diff_toBetween=="lower_thanBetween"))/nrow(input),
                     ),
                   
                   grp=c("a.pos.cor","b.non.cor","c.neg.cor",
                         "a.incr.trend","b.stable.trend","c.decr.trend",
                         "a.pos.cor","b.non.cor","c.neg.cor" #,
                         # "a.higher.trend","b.same.trend","c.lower.trend"
                         ),
                   col=c("high","same","low",
                         "high","same","low",
                         "high","same","low" #,
                         # "high","same","low"
                         ))

  spec.cols<-c(cols.species[unique(input$species)], "#D73027", "#BBBBBB","#26466D")
  names(spec.cols)<-c(unique(input$species),"high","same","low")
  
  g<-ggplot(show)
  g<-g+geom_bar(aes(x=x,y=y, group=grp, fill=col),stat="identity", alpha=0.5)
  g<-g+scale_fill_manual(breaks=names(spec.cols), values=spec.cols)
  g<-g+scale_y_continuous("Proportion of sites (%)",limits=c(0,1),breaks=seq(0,1,0.2),labels=formatC(seq(0,1,0.2),format="f",digits=1))
  # g<-g+scale_x_continuous("",limits=c(0.5,3.5),breaks=c(1:3),labels=c("WithinT","DifToBtw", "Cor"))
  g<-g+scale_x_continuous("",limits=c(0.5,3.5),breaks=c(1:3),labels=c("BtwT","WithinT", "Cor"))
  g<-my.theme(g,"none")
  g
}
plot.booted.cv<-function(input, input.crn, input_alldata){
  
  ## Testing arguments
  # input=subset(booted.data, species=="PCAB")
  # input.crn=subset(df.crn, species=="PCAB")
  # input_alldata=subset(clim.dataset, species=="PCAB")

  input.crn<-subset(input.crn,year>=1960 & year<2018)

  l <- lme(fixed = cv_RWI ~ year + sample_depth, random = ~1 | site_code, correlation = corAR1(form = ~ year), data = input_alldata)
  
  pred.l<-data.frame(year=c(1961:2017),sample_depth=20)
  predicted<-predict(l, newdata=pred.l, level = 0, se = T)
  pred.l$fit<-predicted$fit
  pred.l$low<-pred.l$fit-predicted$se.fit
  pred.l$upper<-pred.l$fit+predicted$se.fit
  
  sd.l<-gls(cv.std.mid~year, correlation = corAR1(form = ~ year),data=subset(input.crn,year>=1960 & year<2018))
  sd.pred.l<-data.frame(year=c(1961:2017),temp=predict(sd.l))
  
  input$rwi_max[which(input$rwi_max>0.89)]<-0.89
  
  g<-ggplot(input)
  
  ## Eventy
  g<-g+annotate("rect",xmin=1971,xmax=1992,ymin=0,ymax=0.89,fill="#BBBBBB",alpha=0.25)
  g<-g+geom_vline(xintercept=c(1992, 2003),linetype="dotted",colour="#D73027",linewidth=0.75,alpha=0.25)
  g<-g+geom_vline(xintercept=c(1996, 2010),linetype="dotted",colour="#26466D",linewidth=0.75,alpha=0.25)
  
  ## Hlavní data
  g<-g+geom_ribbon(aes(x=year,ymin=rwi_min,ymax=rwi_max),linetype="dotted",alpha=0.2, fill=cols.species[unique(input$species)])
  g<-g+geom_point(aes(x=year,y=rwi_mid), colour=cols.species[unique(input$species)])
  g<-g+geom_ribbon(data=input.crn,mapping=aes(x=year,ymin=cv.std.min,ymax=cv.std.max),linetype="dotted",alpha=0.2, fill="#777777",inherit.aes=F)
  g<-g+geom_point(data=input.crn,mapping=aes(x=year,y=cv.std.mid), colour="#777777",inherit.aes=F)
  
  ## Lineární regrese
  g<-g+geom_ribbon(data=pred.l, aes(x=year,ymin=low,ymax=upper), fill=cols.species[unique(input$species)], inherit.aes=F, alpha=0.25)
  g<-g+geom_line(data=pred.l, aes(x=year,y=fit), colour=cols.species[unique(input$species)], linetype="solid", inherit.aes=F,linewidth=0.5)
  
  g<-g+geom_line(data=sd.pred.l, aes(x=year,y=temp), colour="#777777", linetype="solid", inherit.aes=F,linewidth=0.5)
  
  ## Zbytek grafu
  g<-g+scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("Mean variance",limits=c(0.0,0.90),breaks=seq(0,10,0.2),labels=formatC(seq(0,10,0.2),format="f",digits=1))
  g<-g+scale_x_continuous("Calendar year",limits=c(1960,2020),breaks=seq(1900,2020,10))
  g<-my.theme(g)
  g
}
plot.booted.cv.2<-function(input, input.crn, input_alldata){
  
  ## Testing arguments
  # input=subset(booted.data, species=="PCAB")
  # input.crn=subset(df.crn, species=="PCAB")
  # input_alldata=subset(clim.dataset, species=="PCAB")
  
  input.crn<-subset(input.crn,year>=1960 & year<2018)
  
  l <- lme(fixed = cv_RWI ~ year + sample_depth, random = ~1 | site_code, correlation = corAR1(form = ~ year), data = input_alldata)
  
  pred.l<-data.frame(year=c(1961:2017),sample_depth=20)
  predicted<-predict(l, newdata=pred.l, level = 0, se = T)
  pred.l$fit<-predicted$fit
  pred.l$low<-pred.l$fit-predicted$se.fit
  pred.l$upper<-pred.l$fit+predicted$se.fit
  
  sd.l<-gls(cv.std.mid~year, correlation = corAR1(form = ~ year),data=subset(input.crn,year>=1960 & year<2018))
  sd.pred.l<-data.frame(year=c(1961:2017),temp=predict(sd.l))
  
  input$rwi_max[which(input$rwi_max>0.89)]<-0.89
  
  g<-ggplot(input)
  
  ## Eventy
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
  
  ## Hlavní data
  g<-g+geom_ribbon(aes(x=year,ymin=rwi_min,ymax=rwi_max),linetype="dotted",alpha=0.2, fill=cols.species[unique(input$species)])
  g<-g+geom_point(aes(x=year,y=rwi_mid), colour=cols.species[unique(input$species)])
  g<-g+geom_ribbon(data=input.crn,mapping=aes(x=year,ymin=cv.std.min,ymax=cv.std.max),linetype="dotted",alpha=0.2, fill="#777777",inherit.aes=F)
  g<-g+geom_point(data=input.crn,mapping=aes(x=year,y=cv.std.mid), colour="#777777",inherit.aes=F)
  
  ## Lineární regrese
  g<-g+geom_ribbon(data=pred.l, aes(x=year,ymin=low,ymax=upper), fill=cols.species[unique(input$species)], inherit.aes=F, alpha=0.25)
  g<-g+geom_line(data=pred.l, aes(x=year,y=fit), colour=cols.species[unique(input$species)], linetype="solid", inherit.aes=F,linewidth=0.5)
  
  g<-g+geom_line(data=sd.pred.l, aes(x=year,y=temp), colour="#777777", linetype="solid", inherit.aes=F,linewidth=0.5)
  
  ## Zbytek grafu
  g<-g+scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_y_continuous("Mean variance",limits=c(0.0,0.90),breaks=seq(0,10,0.2),labels=formatC(seq(0,10,0.2),format="f",digits=1))
  g<-g+scale_x_continuous("Calendar year",limits=c(1960,2020),breaks=seq(1900,2020,10))
  g<-my.theme(g)
  g
}
plot.panel<-function(sp, mod_dataset){
  
  ## Testing arguments
  # sp = "ABAL"
  # mod_dataset = clim.dataset
  
  ggarrange(plot.curve.differences(subset(curve.differences, species==sp)),
            plot.booted.cv(subset(booted.data,species==sp),
                           subset(df.crn, species==sp),
                           subset(mod_dataset, species==sp)),
            nrow=1,ncol=2,align="hv",widths=c(0.20,0.80))
}
plot.panel.2<-function(sp, mod_dataset){
  
  ## Testing arguments
  # sp = "ABAL"
  # mod_dataset = clim.dataset
  
  ggarrange(plot.curve.differences(subset(curve.differences, species==sp)),
            plot.booted.cv.2(subset(booted.data,species==sp),
                           subset(df.crn, species==sp),
                           subset(mod_dataset, species==sp)),
            nrow=1,ncol=2,align="hv",widths=c(0.20,0.80))
}

## -----------------------------------------------------------------------------
## Data calculations ####
booted.data<-boot.data.general(clim.dataset)

curve.differences<-rbind(calculate.curve.differences(subset(df.crn,species=="ABAL"), subset(clim.dataset, species=="ABAL")),
                         calculate.curve.differences(subset(df.crn,species=="PCAB"), subset(clim.dataset, species=="PCAB")),
                         calculate.curve.differences(subset(df.crn,species=="PISY"), subset(clim.dataset, species=="PISY")),
                         calculate.curve.differences(subset(df.crn,species=="FASY"), subset(clim.dataset, species=="FASY")),
                         calculate.curve.differences(subset(df.crn,species=="QUSP"), subset(clim.dataset, species=="QUSP")))

df.crn.cutted <- df.crn.all[which(df.crn.all$site_code %in% clim.dataset$site_code),]

df.crn.cutted <- unify.categories(df.crn.cutted)

curve.slopes<-rbind(calculate.curve.slopes(subset(df.crn.cutted, species=="ABAL")),
                    calculate.curve.slopes(subset(df.crn.cutted, species=="PCAB")),
                    calculate.curve.slopes(subset(df.crn.cutted, species=="PISY")),
                    calculate.curve.slopes(subset(df.crn.cutted, species=="FASY")),
                    calculate.curve.slopes(subset(df.crn.cutted, species=="QUSP")))

curve.slopes <- merge(curve.slopes, site.list[, c("site_code", "elevation")], by = "site_code")

curve.slopes <- unify.categories(curve.slopes)

## -----------------------------------------------------------------------------
## Data analysis ####
## Calculating and printing differences in variance ####
print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("--------------------------- AOV a between site -------------------------")
print("------------------------------------------------------------------------")
dunn.results<-dunnTest(cv.std.mid~species, data=df.crn[which(df.crn$species %in% c("ABAL","PCAB","PISY","FASY","QUSP")),])$res
pvals<-rep("non-significant", nrow(dunn.results))
pvals[which(dunn.results$P.adj<0.05)]<-"significant"
dunn.results$P.adj<-pvals

round(mean(df.crn[which(df.crn$species %in% c("ABAL")),]$cv.std.mid),2) ; round(sd(df.crn[which(df.crn$species %in% c("ABAL")),]$cv.std.mid),2)
round(mean(df.crn[which(df.crn$species %in% c("PCAB")),]$cv.std.mid),2) ; round(sd(df.crn[which(df.crn$species %in% c("PCAB")),]$cv.std.mid),2)
round(mean(df.crn[which(df.crn$species %in% c("PISY")),]$cv.std.mid),2) ; round(sd(df.crn[which(df.crn$species %in% c("PISY")),]$cv.std.mid),2)
round(mean(df.crn[which(df.crn$species %in% c("FASY")),]$cv.std.mid),2) ; round(sd(df.crn[which(df.crn$species %in% c("FASY")),]$cv.std.mid),2)
round(mean(df.crn[which(df.crn$species %in% c("QUSP")),]$cv.std.mid),2) ; round(sd(df.crn[which(df.crn$species %in% c("QUSP")),]$cv.std.mid),2)

print(dunn.results)


print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("--------------------------- AOV a within site --------------------------")
print("------------------------------------------------------------------------")
dunn.results<-dunnTest(cv_RWI~species, data=clim.dataset[which(clim.dataset$species %in% c("ABAL","PCAB","PISY","FASY","QUSP")),])$res
pvals<-rep("non-significant", nrow(dunn.results))
pvals[which(dunn.results$P.adj<0.05)]<-"significant"
dunn.results$P.adj<-pvals

round(mean(clim.dataset[which(clim.dataset$species %in% c("ABAL")),]$cv_RWI),2) ; round(sd(clim.dataset[which(clim.dataset$species %in% c("ABAL")),]$cv_RWI),2)
round(mean(clim.dataset[which(clim.dataset$species %in% c("PCAB")),]$cv_RWI),2) ; round(sd(clim.dataset[which(clim.dataset$species %in% c("PCAB")),]$cv_RWI),2)
round(mean(clim.dataset[which(clim.dataset$species %in% c("PISY")),]$cv_RWI),2) ; round(sd(clim.dataset[which(clim.dataset$species %in% c("PISY")),]$cv_RWI),2)
round(mean(clim.dataset[which(clim.dataset$species %in% c("FASY")),]$cv_RWI),2) ; round(sd(clim.dataset[which(clim.dataset$species %in% c("FASY")),]$cv_RWI),2)
round(mean(clim.dataset[which(clim.dataset$species %in% c("QUSP")),]$cv_RWI),2) ; round(sd(clim.dataset[which(clim.dataset$species %in% c("QUSP")),]$cv_RWI),2)

print(dunn.results)


## Calculating and printing trends ####
print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("---------------------- Trendu within a between site --------------------")
print("------------------------------------------------------------------------")
sumarize.trends.results("ABAL")
sumarize.trends.results("PCAB")
sumarize.trends.results("PISY")
sumarize.trends.results("FASY")
sumarize.trends.results("QUSP")


print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("--------------------- Proportion trends within site --------------------")
print("------------------------------------------------------------------------")
sumarize.trends.directionss("ABAL")
sumarize.trends.directionss("PCAB")
sumarize.trends.directionss("PISY")
sumarize.trends.directionss("FASY")
sumarize.trends.directionss("QUSP")


print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("----------------- Correlation within and between site ------------------")
print("------------------------------------------------------------------------")
sumarize.correlationss("ABAL")
sumarize.correlationss("PCAB")
sumarize.correlationss("PISY")
sumarize.correlationss("FASY")
sumarize.correlationss("QUSP")


## -----------------------------------------------------------------------------
## Figure output ####
figure<-ggarrange(ggarrange(plot.boxplot.variance.betweensite(df.crn),
                            plot.boxplot.variance.withinsite(clim.dataset),
                            nrow=1,ncol=2,labels=LETTERS[1:2],align="hv"),
                  ggarrange(plot.panel.2("ABAL", clim.dataset),
                            plot.panel("PCAB", clim.dataset),
                            plot.panel("PISY", clim.dataset),
                            plot.panel("FASY", clim.dataset),
                            plot.panel("QUSP", clim.dataset),
                            nrow=5,ncol=1,labels=LETTERS[3:7],align="hv",
                            common.legend=T, legend="bottom"),
                  nrow=2,ncol=1,labels=LETTERS[1],heights=c(0.2,0.8))




















