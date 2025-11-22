## Functions ####

## Plots area of distribution for specific species
#
# sp_name - species name (abbreviation)
# sp_dist - species distribution data (shapefile)
#
plot.sp.dist <- function(sp_name, sp_dist){
  g<-ggplot()
  # State borders
  g <- g + geom_sf(data = europe_plus, fill = NA, color = "black", size = 0.3)
  g <- g + geom_sf(data = cz, fill = "black", color = "black", linewidth = 0.5)
  # Area of distribution
  g <- g + geom_sf(data = sp_dist, fill = alpha(cols.species[sp_name], 0.5), color = cols.species[sp_name], linewidth = 0.1)
  # Graph adjusting
  g <- g + scale_x_continuous("", limits = c(-10.5, 45.0), breaks = seq(-10, 45, 5))
  g <- g + scale_y_continuous("", limits = c(33.0, 71.5), breaks = seq(35, 70, 5))
  g <- g + coord_sf(xlim = c(-10.5, 45.0), ylim = c(33.0, 71.5), expand = FALSE)
  g <- g + theme_classic(base_size = 12) 
  g <- g + theme(
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "right") 
  g
}

## Map of the region ####

# Load elevation raster
hill <- rast("Input_data/Hillshade/eudem_hlsd_3035_europe.tif")

# Define map range in WGS84
bbox_wgs84 <- st_as_sfc(st_bbox(c(xmin = 10, xmax = 20, ymin = 48, ymax = 52), crs = 4326))
# Adjusting coordinate systems (EPSG:3035)
bbox_3035 <- st_transform(bbox_wgs84, crs = crs(hill))
bbox_vect <- vect(bbox_3035)
ext_3035 <- ext(bbox_vect)
# Cut raster (to decrease calculation time)
hill_crop <- crop(hill, ext_3035)
hill_crop <- aggregate(hill_crop, fact = 8, fun = mean)
# Reprojection of the elevation raster to WGS84 EPSG:4326
hill_wgs84 <- project(hill_crop, "EPSG:4326")
# Transfer raster to data.frame for ggplot2 plotting
hill_df <- as.data.frame(hill_wgs84, xy = TRUE)
names(hill_df)[3] <- "shade"
# Get elevations from raster (keep EPSG:4326)
bbox_wgs84_sp <- as(st_as_sf(bbox_wgs84), "Spatial")
elev <- get_elev_raster(locations = st_as_sf(bbox_wgs84), z = 9, clip = "bbox")
# Decrease resolution terra::rast and tranformation to data.frame
elev_terra <- rast(elev)
elev_crop <- crop(elev_terra, ext(10, 20, 48, 52))
elev_crop <- aggregate(elev_crop, fact = 8, fun = mean)

elev_df <- as.data.frame(elev_crop, xy = TRUE)
names(elev_df)[3] <- "elev"

europe <- ne_countries(continent = "Europe", scale = "large")

g<-ggplot(site.list)
# Elevation as background
g <- g + geom_raster(data = elev_df, aes(x = x, y = y, fill = elev), alpha = 0.35)
g <- g + scale_fill_gradientn(colours = terrain.colors(7), guide = "none")
g <- g + geom_raster(data = hill_df, aes(x = x, y = y), fill = "black", alpha = 0.20)
# New colour scale for points
g <- g + ggnewscale::new_scale_fill()
# State borders
g <- g + geom_sf(data = europe, fill = NA, color = "black", size = 0.3)
# Plotting points (everlapping points are jittered)
g <- g + geom_jitter(aes(x = lon, y = lat, fill = species),
                     shape = 21, colour = "black", stroke = 0.3,
                     size = 2, width = 0.01, height = 0.01)
g <- g + scale_fill_manual(values = cols.species, name = "Species",
                           breaks = names(cols.species))
# Graph scaling, theme adjusting, adding necessary things, etc.
g <- g+scale_x_continuous("",limits=c(11.5,19.5),breaks=seq(0,20,1))
g <- g+scale_y_continuous("",limits=c(48.25,51.25),breaks=seq(0,100,1))
g <- g + coord_sf(xlim = c(11.5, 19.25), ylim = c(48.25, 51.25), expand = FALSE)
g <- g + theme_classic(base_size = 12) 
g <- g + theme(
  panel.grid.major = element_line(color = "grey80", size = 0.2),
  panel.grid.minor = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_text(color = "black"),
  legend.position = "right") 
# North arrow
g <- g + annotation_north_arrow(
  location = "tr",  # "tr" = top right
  which_north = "true",
  style = north_arrow_minimal(),
  height = unit(1.2, "cm"),
  width = unit(1.2, "cm"),
  pad_x = unit(0.2, "cm"),
  pad_y = unit(0.2, "cm"))
g <- g + annotation_scale(
  location = "bl",                     
  width_hint = 0.25,                  
  line_width = 0.8,
  text_cex = 0.8,
  pad_x = unit(0.2, "cm"),             
  pad_y = unit(0.2, "cm"),            
  style = "bar",                       
  box_fill = alpha("white", 0.5),      
  box_colour = "White"              
)
mapcr <- g

# deleting unnecessary files
rm(g, 
   hill, hill_crop, hill_wgs84, hill_df,
   elev, elev_terra, elev_crop, elev_df,
   bbox_wgs84, bbox_wgs84_sf,
   bbox_3035, bbox_vect, ext_3035)

## Species distributions ####
# Data loadings
world <- ne_countries(scale = "large", returnclass = "sf")
europe_plus <- world[world$sovereignt %in% c(
  "Austria", "Belgium", "Bosnia and Herzegovina", "Bulgaria", "Croatia",
  "Czechia", "Denmark", "Estonia", "Finland", "France", "Germany",
  "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Kosovo",
  "Latvia", "Lithuania", "Luxembourg", "Moldova", "Montenegro",
  "Netherlands", "North Macedonia", "Norway", "Poland", "Portugal",
  "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
  "Switzerland", "Ukraine", "United Kingdom", "Albania", "Turkey", 
  "Cyprus", "Russia"
), ]
cz <- world[world$sovereignt %in% c(
  "Czechia"
), ]

# Plots  species distribution maps
dist_maps<-ggarrange(plot.sp.dist("ABAL", st_read("Input_data/Species_distribution/Abies alba/shapefiles/Abies_alba_plg.shp")),
                     plot.sp.dist("PCAB", st_read("Input_data/Species_distribution/Picea abies/shapefiles/Picea_abies_plg_clip.shp")),
                     plot.sp.dist("PISY", st_read("Input_data/Species_distribution/Pinus sylvestris/shapefiles/Pinus_sylvestris_plg_clip.shp")),
                     plot.sp.dist("FASY", st_read("Input_data/Species_distribution/Fagus sylvatica/shapefiles/Fagus_sylvatica_sylvatica_plg_clip.shp")),
                     plot.sp.dist("QUSP", st_union(rbind(st_read("Input_data/Species_distribution/Quercus petraea/shapefiles/Quercus_petraea_plg_clip.shp"),
                                                         st_read("Input_data/Species_distribution/Quercus robur/shapefiles/Quercus_robur_plg_clip.shp")))),
                     nrow=1,ncol=5)



## Graph of Tepmperature and CWB data ####
inds<-names(climadata)[which(substr(names(climadata),8,11) %in% c("QUSP", "FASY", "PISY", "PCAB", "ABAL","QURO","QUPE","QUsp"))]
c.data<-climadata[inds]
c.data<-do.call(rbind,c.data)
c.data$site<-substr(rownames(c.data),1,11)
c.data<-aggregate(Temp~Year+site,data=c.data,FUN=mean)

c.data<-c.data[which(c.data$site %in% considered_sites),]

show.df<-data.frame(Year=c(1961:2017),
                    minT=NA,midT=NA,maxT=NA,
                    minSPEI=NA,midSPEI=NA,maxSPEI=NA)

for(i in show.df$Year){
  # print(i)
  sub.clim<-subset(c.data,Year==i)
  show.df[which(show.df$Year==i),c("minT","midT","maxT")]<-quantile(apply(replicate(1000,sample(sub.clim$Temp,nrow(sub.clim),replace=T)),2,mean),probs=c(0.025,0.50,0.975))
  
  sub.spei<-subset(bals,year==i)
  sub.spei[sapply(sub.spei, is.infinite)] <- NA
  sub.spei<-na.omit(sub.spei)
  show.df[which(show.df$Year==i),c("minSPEI","midSPEI","maxSPEI")]<-quantile(apply(replicate(1000,sample(sub.spei$mean_bal,nrow(sub.spei),replace=T)),2,mean),probs=c(0.025,0.50,0.975))
}

l.temp<-lm(midT~Year, data=subset(show.df,Year>=1960 & Year<2018)) ; pred.temp<-data.frame(year=c(1961:2017),temp=predict(l.temp))
l.spei<-lm(midSPEI~Year, data=subset(show.df,Year>=1960 & Year<2018)) ; pred.spei<-data.frame(year=c(1961:2017),spei=predict(l.spei))

g<-ggplot(show.df)
g<-g+annotate("rect",xmin=1971,xmax=1992,ymin=0,ymax=10,fill="#BBBBBB",alpha=0.25)
g<-g+geom_rect_pattern(
  xmin = 1992, xmax = 1998, ymin = 0, ymax = 10,
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
g<-g+geom_vline(xintercept=c(1992, 2003),linetype="dotted",size=1,colour="#D73027")
g<-g+geom_vline(xintercept=c(1995, 2010),linetype="dotted",size=1,colour="#26466D")
g<-g+geom_ribbon(aes(x=Year,ymin=minT,ymax=maxT),fill="#D73027",alpha=0.25)
g<-g+geom_line(aes(x=Year,y=midT),colour="#D73027",size=1.5)

g<-g+geom_ribbon(aes(x=Year,ymin=(minSPEI)/10+5,ymax=(maxSPEI)/10+5),fill="#26466D",alpha=0.25)
g<-g+geom_line(aes(x=Year,y=(midSPEI)/10+5),colour="#26466D",size=1.5)

g<-g+geom_line(data=pred.temp, aes(x=year,y=temp), colour="#D73027", linetype="solid", inherit.aes=F)
g<-g+geom_line(data=pred.spei, aes(x=year,y=(spei/10+5)), colour="#26466D", linetype="solid", inherit.aes=F)
g<-g+scale_y_continuous(
  # Features of the first axis
  name = "Temperature (°C)",
  limits = c(0,10),
  breaks=seq(0,10,1),
  
  # Add a second axis and specify its features
  sec.axis = sec_axis( trans=~.*10-5, name="CWB", breaks=seq(0,100,10), labels=seq(-50,50,10))
)
g<-g+scale_x_continuous("Year",limits=c(1960,2020),breaks=seq(0,3000,10))
g<- g + theme_classic() 
g<- g + theme(axis.line.x = element_line(colour="black"), 
              axis.text.x = element_text(colour="black"), 
              axis.line.y = element_line(colour="black"), 
              axis.text.y = element_text(colour="black"), 
              axis.ticks = element_line(color = "black"), 
              legend.position = "none") 
g.cwb.temp<-g

## MAT and MAP Graph ####

g<-ggplot(clim)
g<-g+geom_point(aes(x=MAT,y=MAP,colour=species),size=1.1)
g<-g+scale_y_continuous("MAP (mm)",limits=c(400,1400),breaks=seq(0,2000,200))
g<-g+scale_x_continuous("MAT (°C)",limits=c(2,10),breaks=seq(-10,30,1))
g<-g+scale_colour_manual(values=cols.species,breaks=cols.species)
g<-g+scale_fill_manual(values=cols.species,breaks=cols.species)
g<- g + theme_classic() 
g<- g + theme(axis.line.x = element_line(colour="black"), 
              axis.text.x = element_text(colour="black"), 
              axis.line.y = element_line(colour="black"), 
              axis.text.y = element_text(colour="black"), 
              axis.ticks = element_line(color = "black"), 
              legend.position = "right") 
g.dots<-g


## Boxplots of CWB in the "extreme" years ####

cwb.1<-bals
cwb.1<-subset(cwb.1, year<2017)

cwb.1$year<-"all"
cwb.1<-rbind(cwb.1,
               bals[which(bals$year %in% c(1992, 1995, 2003,2010)),])
cwb.1$x<-1
cwb.1$x[which(cwb.1$year==2010)]<-5
cwb.1$x[which(cwb.1$year==2003)]<-4
cwb.1$x[which(cwb.1$year==1995)]<-3
cwb.1$x[which(cwb.1$year==1992)]<-2

cwb.1$event<-"all"
cwb.1$event[which(cwb.1$year==2010)]<-"moist"
cwb.1$event[which(cwb.1$year==2003)]<-"dry"
cwb.1$event[which(cwb.1$year==1995)]<-"moist"
cwb.1$event[which(cwb.1$year==1992)]<-"dry"

cwb.1$species<-substr(cwb.1$site,8,11)
cwb.1$species[which(cwb.1$species %in% c("QUSP","QURO","QUPE","QUsp"))]<-"QUSP"
cwb.1<-cwb.1[which(cwb.1$species %in% c("QUSP", "FASY", "PISY", "PCAB", "ABAL")),]

cwb.1<-cwb.1[which(cwb.1$site %in% considered_sites),]
  
g<-ggplot(cwb.1)
g<-g+geom_boxplot(aes(x=as.factor(x),y=mean_bal ,fill=event),colour="#000000",size=0.66)
g<-g+scale_x_discrete("Year",labels=c("1961-2017","1992","1995","2003","2010"))
g<-g+scale_y_continuous("CWB",limits=c(-100,150),breaks=seq(-100,150,25))
g<-g+scale_fill_manual(values=c("#D73027","#BBBBBB","#26466D"),breaks=c("all","dry","moist"))
g<- g + theme_classic() 
g<- g + theme(axis.line.x = element_line(colour="black"), 
              axis.text.x = element_text(colour="black"), 
              axis.line.y = element_line(colour="black"), 
              axis.text.y = element_text(colour="black"), 
              axis.ticks = element_line(color = "black"), 
              legend.position = "bottom") 
g.box.cwb.sucho<-g

## Boxplots of Temperature in the "extreme" years ####

temps.1<-temps
temps.1<-subset(temps.1, year<2017)

temps.1$year<-"all"
temps.1<-rbind(temps.1,
               temps[which(temps$year %in% c(1992, 1995, 2003,2010)),])
temps.1$x<-1
temps.1$x[which(temps.1$year==2010)]<-5
temps.1$x[which(temps.1$year==2003)]<-4
temps.1$x[which(temps.1$year==1995)]<-3
temps.1$x[which(temps.1$year==1992)]<-2

temps.1$event<-"all"
temps.1$event[which(temps.1$year==2010)]<-"moist"
temps.1$event[which(temps.1$year==2003)]<-"dry"
temps.1$event[which(temps.1$year==1995)]<-"moist"
temps.1$event[which(temps.1$year==1992)]<-"dry"

temps.1$species<-substr(temps.1$site,8,11)
temps.1$species[which(temps.1$species %in% c("QUSP","QURO","QUPE","QUsp"))]<-"QUSP"
temps.1<-temps.1[which(temps.1$species %in% c("QUSP", "FASY", "PISY", "PCAB", "ABAL")),]

temps.1<-temps.1[which(temps.1$site %in% considered_sites),]

g<-ggplot(temps.1)
g<-g+geom_boxplot(aes(x=as.factor(x),y=mean_temp ,fill=event),colour="#000000",size=0.66)
g<-g+scale_x_discrete("Year",labels=c("1961-2017","1992","1995","2003","2010"))
g<-g+scale_y_continuous("Temperature (°C)",limits=c(0,20),breaks=seq(0,100,2))
g<-g+scale_fill_manual(values=c("#D73027","#BBBBBB","#26466D"),breaks=c("all","dry","moist"))
g<- g + theme_classic() 
g<- g + theme(axis.line.x = element_line(colour="black"), 
              axis.text.x = element_text(colour="black"), 
              axis.line.y = element_line(colour="black"), 
              axis.text.y = element_text(colour="black"), 
              axis.ticks = element_line(color = "black"), 
              legend.position = "bottom") 
g.box.temp.sucho<-g


## Final graph assembling ####
figure<-ggarrange(dist_maps,
                  mapcr,
                  g.cwb.temp,
                  ggarrange(g.dots,g.box.cwb.sucho,g.box.temp.sucho,nrow=1,ncol=3,labels=c("D","E","F"),align="hv", widths=c(0.4,0.30,0.30)),
                  labels=c("A","B","C",""),nrow=4,ncol=1, heights=c(0.14,0.46,0.20,0.20))

# bottom_row <- ggarrange(
#   g.dots, g.box.cwb.sucho, g.box.temp.sucho,
#   nrow = 1, ncol = 3,
#   labels = c("D","E","F"),
#   widths = c(0.4, 0.30, 0.30)   # bez align
# )
# 
# figure <- ggarrange(
#   dist_maps,
#   mapcr,
#   g.cwb.temp,
#   labels  = c("A","B","C"),
#   nrow    = 3,
#   ncol    = 1,
#   heights = c(0.14, 0.46, 0.20)
# )
# 
# figure <- ggarrange(
#   dist_maps,
#   mapcr,
#   g.cwb.temp,
#   bottom_row,
#   labels  = c("A","B","C",""),
#   nrow    = 4,
#   ncol    = 1,
#   heights = c(0.14, 0.46, 0.20, 0.20)
# )