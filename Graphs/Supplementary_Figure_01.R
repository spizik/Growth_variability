## Functions ####
## Data preparation and calculations
# Builds a site-level chronology dataset by selecting single-species sites and restricting years to 1961–2017.
# Returns a long-format data frame with site, species, year, standardized and residual indices (std, res), and placeholder coordinate columns.
prepare.crn.data<-function(input){
  
  # input=df.crn.cutted
  
  output.df <- data.frame(site = character(),
                          species = character(),
                          year = numeric(),
                          std = numeric(),
                          res = numeric(),
                          x = numeric(),
                          y = numeric())
  
  for(i in unique(input$site_code)){
    
    # i=names(input)[1]
    # print(i)
    
    if(length(site.list$species[which(site.list$site_code == i)]) == 1){
      sub<-subset(input, site_code == i)
      sub<-sub[which(sub$year %in% c(1961:2017)),]
      
      temp <- data.frame(site = i,
                         species = site.list$species[which(site.list$site_code == i)],
                         year = sub$year,
                         std = sub$std,
                         res = sub$res#,
                         # x = site.list$lon[which(site.list$site_code == i)],
                         # y = site.list$lat[which(site.list$site_code == i)]
      )
      
      output.df <- rbind(output.df, temp)
    }
  }
  
  return(output.df)
}

# Filters the input climate dataset to the target analysis period and keeps only key variability and location variables.
# Returns a data frame with site code, species, year, RWI variability metrics (sd_RWI, cv_RWI), and site coordinates.
prepare.variability.data<-function(input){
  
  # input = clim.dataset
  
  input <- subset(input, year > 1960 & year < 2024)
  output.df <- input[,c("site_code", "species", "year", "sd_RWI", "cv_RWI", "x", "y")]
  
  return(output.df)
}

## Figure functions
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
plot.crn.data<-function(df){
  
  ## Testing arguments
  # df=subset(variogram_crn_data, species=="PISY")
  
  df<-merge(df, site.list[,c("site_code", "elevation")], by.x="site", by.y="site_code")
  df$elev_cat <- cut(df$elevation, breaks=c(0, 300, 500, 700, 900, 1100, 1600), labels=c("<300", "300-500", "500-700", "700-900", "900-1100", ">1100"))
  
  booted.df<-data.frame(calendar_year = c(min(df$year):max(df$year)),
                        min = NA, mid = NA, max = NA)
  
  for(i in booted.df$calendar_year){
    sub <- subset(df,year== i)
    if(nrow(sub)>=5){
      booted.df[which(booted.df$calendar_year == i), c("min", "mid", "max")] <- quantile(apply(replicate(1000, sample(sub$std, nrow(sub), T)), 2, mean), probs=c(0.025, 0.500, 0.975))
    }
  }
  
  gls_model_trend <- gls(mid ~ calendar_year, correlation = corAR1(form = ~ calendar_year), data = booted.df)
  trend_line <- data.frame(calendar_year = booted.df$calendar_year,
                           trend = predict(gls_model_trend))
  
  g<-ggplot(df)
  g<-g+geom_line(aes(x=year, y=std, group=site, colour=elev_cat), 
                 linewidth=0.5, 
                 alpha=0.25)
  g<-g+geom_ribbon(data= booted.df, 
                   mapping= aes(x=calendar_year, ymin=min, ymax=max), 
                   linewidth=0.5, 
                   alpha=0.5, 
                   inherit.aes=F, 
                   fill=cols.species[unique(df$species)])
  g<-g+geom_line(data = trend_line, 
                 mapping = aes(x=calendar_year, y=trend,), 
                 linewidth = 1.2, 
                 linetype = "dashed", 
                 colour=cols.species[unique(df$species)], 
                 inherit.aes = F)
  g<-g+scale_colour_manual(breaks=names(elev_colors), values=elev_colors)
  g<-g+scale_x_continuous("Calendar year", limits=c(1960,2015), breaks=seq(0,3000,10))
  g<-g+scale_y_continuous("Std. chronology", limits=c(0.00,2.25), breaks=seq(0,5,0.25), labels=formatC(seq(0,5,0.25),format="f",digits=2))
  g<-my.theme(g,"none")
  g
}
plot.cv.data<-function(df){
  
  ## Testing arguments
  # df=subset(variogram_site_data, species=="PISY")
  
  df<-merge(df, site.list[,c("site_code", "elevation")], by.x="site_code", by.y="site_code")
  df$elev_cat <- cut(df$elevation, breaks=c(0, 300, 500, 700, 900, 1100, 1600), labels=c("<300", "300-500", "500-700", "700-900", "900-1100", ">1100"))
  
  booted.df<-data.frame(calendar_year = c(min(df$year):max(df$year)),
                        min = NA, mid = NA, max = NA)
  
  for(i in booted.df$calendar_year){
    sub <- subset(df,year== i)
    if(nrow(sub)>=5){
      booted.df[which(booted.df$calendar_year == i), c("min", "mid", "max")] <- quantile(apply(replicate(1000, sample(sub$cv_RWI, nrow(sub), T)), 2, mean), probs=c(0.025, 0.500, 0.975))
    }
  }
  
  gls_model_trend <- gls(mid ~ calendar_year, correlation = corAR1(form = ~ calendar_year), data = booted.df)
  trend_line <- data.frame(calendar_year = booted.df$calendar_year,
                           trend = predict(gls_model_trend))
  
  g<-ggplot(df)
  g<-g+geom_line(aes(x=year, y=cv_RWI, group=site_code, colour=elev_cat), 
                 linewidth=0.5, 
                 alpha=0.25)
  g<-g+geom_ribbon(data= booted.df, 
                   mapping= aes(x=calendar_year, ymin=min, ymax=max), 
                   linewidth=0.5, 
                   alpha=0.5, 
                   inherit.aes=F, 
                   fill=cols.species[unique(df$species)])
  g<-g+geom_line(data = trend_line, 
                 mapping = aes(x=calendar_year, y=trend,), 
                 linewidth = 1.2, 
                 linetype = "dashed", 
                 colour=cols.species[unique(df$species)], 
                 inherit.aes = F)
  g<-g+scale_colour_manual(breaks=names(elev_colors), values=elev_colors)
  g<-g+scale_x_continuous("Calendar year", limits=c(1960,2015), breaks=seq(0,3000,10))
  g<-g+scale_y_continuous("Withinsite variability", limits=c(0.0,1.5), breaks=seq(0,5,0.2), labels=formatC(seq(0,5,0.2),format="f",digits=1))
  g<-my.theme(g,"none")
  g
}

## Calculations ####
df.crn.cutted <- df.crn.all[which(df.crn.all$site_code %in% clim.dataset$site_code),]
variogram_crn_data<-prepare.crn.data(df.crn.cutted)
variogram_site_data<-prepare.variability.data(clim.dataset)

variogram_crn_data<-subset(variogram_crn_data, year>1960 & year<2017)
variogram_site_data<-subset(variogram_site_data, year>1960 & year<2017)
variogram_crn_data <- variogram_crn_data[which(variogram_crn_data$site %in% variogram_site_data$site_code),]

## Figure ####
figure<-ggarrange(plot.crn.data(df=subset(variogram_crn_data, species=="ABAL")), plot.cv.data(df=subset(variogram_site_data, species=="ABAL")),
                  plot.crn.data(df=subset(variogram_crn_data, species=="PCAB")), plot.cv.data(df=subset(variogram_site_data, species=="PCAB")),
                  plot.crn.data(df=subset(variogram_crn_data, species=="PISY")), plot.cv.data(df=subset(variogram_site_data, species=="PISY")),
                  plot.crn.data(df=subset(variogram_crn_data, species=="FASY")), plot.cv.data(df=subset(variogram_site_data, species=="FASY")),
                  plot.crn.data(df=subset(variogram_crn_data, species=="QUSP")), plot.cv.data(df=subset(variogram_site_data, species=="QUSP")),
                  nrow=5, ncol=2, align="hv", labels=LETTERS[1:10])


