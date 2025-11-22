

## Functions ####
## Data preparation
# Creates a balanced event-centered dataset by combining pre-event averaged variability with post-event values and adjusting x-coordinates by species offsets.
# Returns a compact dataset suitable for plotting standardized temporal profiles around disturbance events.
omit.data<-function(evt,sites){
  
  ## Testing arguments
  # evt=2003
  # sites=fasy_mod_dataset
  
  temp<-subset(sites,year>=evt-3 & year<=evt-1)[,c("site_code", "species", "cv_RWI")]
  temp<-aggregate(cv_RWI~site_code+species,data=temp,FUN=mean)
  temp$year<-evt-1
  temp<-temp[,c("site_code", "species", "year", "cv_RWI")]
  temp$x<--1
  
  sites<-subset(sites,year>=evt & year<=evt+3)[,c("site_code", "species", "year", "cv_RWI")]
  sites$x<-sites$year-evt
  
  spp<-unique(sites$species)
  
  sites<-rbind(temp, sites)
  
  if(spp=="ABAL") sites$x<-sites$x-0.30
  if(spp=="PCAB") sites$x<-sites$x-0.15
  if(spp=="PISY") sites$x<-sites$x-0.00
  if(spp=="FASY") sites$x<-sites$x+0.15
  if(spp=="QUSP") sites$x<-sites$x+0.30
  
  return(sites)
}

# Generates bootstrap confidence intervals for yearly variability for a single species and re-centers x-coordinates for multi-species plotting.
# Returns a data frame containing min, mid, and max bootstrap envelopes aligned to event-centered time.
boot.line<-function(input.points, spp){
  
  ## Testing arguments
  # input.points=withinsite_2010
  # spp="ABAL"
  
  dta<-subset(input.points,species==spp)
  
  input.line<-data.frame(Species=rep(spp,5),
                         year=c(min(dta$year):max(dta$year)),
                         min=NA,mid=NA,max=NA)
  
  for(i in 1:nrow(input.line)){
    sub.data<-subset(dta,year==input.line$year[i])
    
    input.line[i,c("min","mid","max")]<-quantile(apply(replicate(1000,sample(sub.data$cv_RWI, nrow(sub.data), T)),2,mean),probs=c(0.025,0.500,0.975))
    
  }
  
  input.line$x<-input.line$year-min(input.line$year)-1
  
  if(spp=="ABAL") input.line$x<-input.line$x-0.30
  if(spp=="PCAB") input.line$x<-input.line$x-0.15
  if(spp=="PISY") input.line$x<-input.line$x-0.00
  if(spp=="FASY") input.line$x<-input.line$x+0.15
  if(spp=="QUSP") input.line$x<-input.line$x+0.30
  
  return(input.line)
}

# Normalizes bootstrap confidence intervals by subtracting the first-year median to express deviations relative to the pre-event state.
# Returns updated min, mid, and max fields centered around zero.
calc.differences<-function(input){
  
  ## Testing arguments
  # input=boot.line(input.points,"ABAL")
  
  output<-input
  output$min<-output$min-output$mid[1]
  output$max<-output$max-output$mid[1]
  output
}

# Caps the lower and upper bounds of the 'min' and 'max' columns to the specified limits.
# Returns the modified data frame with values truncated to the given range.
cut.our.range.vals<-function(input,min_val=-0.2,max_val=0.2){
  # input=input.points
  # min_val=-0.4
  # max_val=0.4
  
  input$min[which(input$min<=min_val)]<-min_val
  input$max[which(input$max>=max_val)]<-max_val
  
  return(input)
}


# Selects CRN data for the given species and ±3 years around the event, computing pre-event means.
# Builds a summary data frame and calculates differences using calc.differences().
prepare.crn.data<-function(evt,sp){
  
  ## Testing arguments
  # evt=2003
  # sp="PCAB"
  
  crn<-subset(df.crn,species==sp)
  crn<-subset(crn,year>=evt-3 & year<=evt+3)
  
  show.crn<-data.frame(year=c((evt-1):(evt+3)),
                       x=c(1:5),
                       min=c(mean(crn$cv.std.min[1:3]),crn$cv.std.min[4:nrow(crn)]),
                       mid=c(mean(crn$cv.std.mid[1:3]),crn$cv.std.mid[4:nrow(crn)]),
                       max=c(mean(crn$cv.std.max[1:3]),crn$cv.std.max[4:nrow(crn)]),
                       species=sp)
  
  show.crn<-calc.differences(show.crn)
  
  return(show.crn)
}

## Tests
# Identifies years showing significant deviations in variability range by testing whether the bootstrapped CI excludes zero.
# Returns a table marking years with an asterisk where either lower or upper CI indicates a meaningful shift.
test.range.differences<-function(input){
  
  ## Testing arguments
  # input=recalculate.range.differences(evt.10.ommited)
  
  output<-input[,c("species","year")]
  output$sig<-NA
  output$sig[which(input$min>0 | input$max<0)]<-"*"
  return(output)
}

# Tests intra-site temporal differences by bootstrapping variability envelopes for all species and checking if CI excludes zero.
# Returns a species–year table with asterisk-marked shifts indicating significant departures from the pre-event baseline.
test.intra.site.differences<-function(input.points){
  
  ## Testing arguments
  # input.points=evt.03.ommited
  
  line.data<-rbind(calc.differences(boot.line(input.points,"ABAL")),
                   calc.differences(boot.line(input.points,"PCAB")),
                   calc.differences(boot.line(input.points,"PISY")),
                   calc.differences(boot.line(input.points,"FASY")),
                   calc.differences(boot.line(input.points,"QUSP")))
  line.data<-cut.our.range.vals(line.data)
  
  output<-line.data[,c("Species","year")]
  output$sig<-NA
  output$sig[which(line.data$min>0 | line.data$max<0)]<-"*"
  return(output)
}

# Detects inter-site deviations by evaluating whether bootstrapped CRN variability envelopes differ from the pre-event baseline.
# Returns a table marking years with significant departures based on CI ranges clipped by user-defined truncation limits.
test.inter.site.differences<-function(input.points){
  
  ## Testing arguments
  # input.points=evt.10
  
  input.points<-cut.our.range.vals(input.points)
  
  output<-input.points[,c("species","year")]
  output$sig<-NA
  output$sig[which(input.points$min>0 | input.points$max<0)]<-"*"
  return(output)
}

## Figure making
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
                     # axis.text.y = element_text(colour="black"),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
plot.site.range.differences<-function(input.points){
  
  ## Testing arguments
  # input.points=evt.03.ommited
  
  line.data<-rbind(calc.differences(boot.line(input.points,"ABAL")),
                   calc.differences(boot.line(input.points,"PCAB")),
                   calc.differences(boot.line(input.points,"PISY")),
                   calc.differences(boot.line(input.points,"FASY")),
                   calc.differences(boot.line(input.points,"QUSP")))
  line.data<-cut.our.range.vals(line.data)
  
  years<-as.numeric(unique(input.points$year))
  years<-years[order(years)]
  years<-c("before",years[2:length(years)])
  
  g<-ggplot(line.data)
  g<-g+geom_hline(yintercept=0,linetype="dotted",colour="#000000")
  g<-g+geom_ribbon(mapping=aes(x=year,ymin=min,ymax=max,group=Species,fill=Species),alpha=0.1)
  g<-g+geom_point(mapping=aes(x=year,y=mid,group=Species,colour=Species),shape=1)
  g<-g+geom_line(mapping=aes(x=year,y=mid,group=Species,colour=Species))
  g<-g+scale_x_continuous("Calendar year", limits=c(min(line.data$year),max(line.data$year)), breaks=c(min(line.data$year):max(line.data$year)), labels=years)
  g<-g+scale_y_continuous("Difference in within-site variability", limits=c(-0.20, 0.20),breaks=seq(-2,2,0.1),labels=formatC(seq(-2,2,0.1),format="f",digits=2))
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g<-my.theme(g)
  g
} # 
plot.chronology.range.differences<-function(input.points){
  
  ## Testing arguments
  # input.points=evt.10
  
  input.points<-cut.our.range.vals(input.points)
  
  years<-as.numeric(unique(input.points$year))
  years<-years[order(years)]
  years<-c("before",years[2:length(years)])
  
  g<-ggplot(input.points)
  g<-g+geom_hline(yintercept=0,linetype="dotted",colour="#000000")
  g<-g+geom_ribbon(mapping=aes(x=year,ymin=min,ymax=max,group=species,fill=species),alpha=0.1)
  g<-g+geom_point(mapping=aes(x=year,y=mid,group=species,colour=species),shape=1)
  g<-g+geom_line(mapping=aes(x=year,y=mid,group=species,colour=species))
  g<-g+scale_x_continuous("Calendar year", limits=c(min(input.points$year),max(input.points$year)), breaks=c(min(input.points$year):max(input.points$year)), labels=years)
  g<-g+scale_y_continuous("Difference in between-site variability", limits=c(-0.20, 0.20),breaks=seq(-2,2,0.1),labels=formatC(seq(-2,2,0.1),format="f",digits=2))
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_colour_manual(values=cols.species, breaks=names(cols.species))
  g<-my.theme(g)
  g
} 

## Data preparation ####
## 1992
evt<-1992
withinsite_1992<-rbind(omit.data(evt,abal_mod_dataset),
                       omit.data(evt,pcab_mod_dataset),
                       omit.data(evt,pisy_mod_dataset),
                       omit.data(evt,fasy_mod_dataset),
                       omit.data(evt,qusp_mod_dataset))

crn_diff_1992<-rbind(prepare.crn.data(evt,"ABAL"),
                     prepare.crn.data(evt,"PCAB"),
                     prepare.crn.data(evt,"PISY"),
                     prepare.crn.data(evt,"FASY"),
                     prepare.crn.data(evt,"QUSP"))

## 1995
evt<-1995
withinsite_1995<-rbind(omit.data(evt,abal_mod_dataset),
                       omit.data(evt,pcab_mod_dataset),
                       omit.data(evt,pisy_mod_dataset),
                       omit.data(evt,fasy_mod_dataset),
                       omit.data(evt,qusp_mod_dataset))

crn_diff_1995<-rbind(prepare.crn.data(evt,"ABAL"),
                     prepare.crn.data(evt,"PCAB"),
                     prepare.crn.data(evt,"PISY"),
                     prepare.crn.data(evt,"FASY"),
                     prepare.crn.data(evt,"QUSP"))

## 2003
evt<-2003
withinsite_2003<-rbind(omit.data(evt,abal_mod_dataset),
                       omit.data(evt,pcab_mod_dataset),
                       omit.data(evt,pisy_mod_dataset),
                       omit.data(evt,fasy_mod_dataset),
                       omit.data(evt,qusp_mod_dataset))

crn_diff_2003<-rbind(prepare.crn.data(evt,"ABAL"),
                     prepare.crn.data(evt,"PCAB"),
                     prepare.crn.data(evt,"PISY"),
                     prepare.crn.data(evt,"FASY"),
                     prepare.crn.data(evt,"QUSP"))

## 2010
evt<-2010
withinsite_2010<-rbind(omit.data(evt,abal_mod_dataset),
                       omit.data(evt,pcab_mod_dataset),
                       omit.data(evt,pisy_mod_dataset),
                       omit.data(evt,fasy_mod_dataset),
                       omit.data(evt,qusp_mod_dataset))

crn_diff_2010<-rbind(prepare.crn.data(evt,"ABAL"),
                     prepare.crn.data(evt,"PCAB"),
                     prepare.crn.data(evt,"PISY"),
                     prepare.crn.data(evt,"FASY"),
                     prepare.crn.data(evt,"QUSP"))



## Tests ####
print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("--------------- between-site variability differences -------------------")
print("------------------------------------------------------------------------")
print("----------------------------- 1992 -------------------------------------")
print(test.inter.site.differences(crn_diff_1992))
print("----------------------------- 2003 -------------------------------------")
print(test.inter.site.differences(crn_diff_2003))
print("----------------------------- 1995 -------------------------------------")
print(test.inter.site.differences(crn_diff_1995))
print("----------------------------- 2010 -------------------------------------")
print(test.inter.site.differences(crn_diff_2010))

print("------------------------------------------------------------------------")
print("------------------------------------------------------------------------")
print("--------------- within-site variability differences --------------------")
print("------------------------------------------------------------------------")
print("----------------------------- 1992 -------------------------------------")
print(test.intra.site.differences(withinsite_1992))
print("----------------------------- 2003 -------------------------------------")
print(test.intra.site.differences(withinsite_2003))
print("----------------------------- 1995 -------------------------------------")
print(test.intra.site.differences(withinsite_1995))
print("----------------------------- 2010 -------------------------------------")
print(test.intra.site.differences(withinsite_2010))

## Figure making ####
figure<-ggarrange(plot.chronology.range.differences(crn_diff_1992),
                  plot.site.range.differences(withinsite_1992),
                  
                  plot.chronology.range.differences(crn_diff_2003),
                  plot.site.range.differences(withinsite_2003),
                  
                  plot.chronology.range.differences(crn_diff_1995),
                  plot.site.range.differences(withinsite_1995),
                  
                  plot.chronology.range.differences(crn_diff_2010),
                  plot.site.range.differences(withinsite_2010),
                  
                  ncol=2,nrow=4,align="hv",labels=LETTERS[1:12], 
                  common.legend=T,legend="bottom")



