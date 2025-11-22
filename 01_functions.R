## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
# libraries ####
## used linraries
library("ape")
library("car")
library("caret") 
library("corrplot")
library("cowplot")
library("dplR")
library("DHARMa")
library("dplyr")
library("effects")
library("elevatr")
library("FSA")
library("iml")
library("ggcorrplot")
library("ggpattern")
library("ggplot2")
library("ggpubr")
library("gratia")
library("ggspatial")
library("gstat")
library("Hmisc")
library("lattice")
library("lme4")
library("lmerTest")
library("lubridate")
library("mgcv")
library("MuMIn")
library("ncdf4")
library("nlme")
library("openxlsx")
library("patchwork")
library("qgam")
library("RColorBrewer")
library("reshape")
library("rnaturalearth")
library("rtrend")
library("sf")
library("sfheaders")
library("sjPlot")
library("spacetime")
library("sp")
library("spdep")
library("SPEI")
library("splines")
library("terra")
library("tidyr")
library("vip")  
library("viridis")
library("zoo")

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ---------------------------- Common Functions ---------------------------- ####
## ---------- load.files ####
## Loads RWL files from a single folder into a list
## Optimized for .rwl and .xlsx
## Each plot is stored as a list item
## Plot name = name of the loaded file
#
# folder - the folder where the files to be loaded are stored
# types  - the file types to be loaded
load.files<-function(folder,
                     types="rwl",
                     separator=";",
                     decimal="."){
  
  ## Testing arguments
  # folder="All_data/TRW_data"
  # types="rwl"
  
  files<-list.files(folder) # creates a vector with the names of all files in the given directory
  files<-files[grep(types,files)] # selects only the files with the given extension (xlsx or rwl)
  output<-list() # creates output data.frame
  
  # loop loading individual files
  for(i in 1:length(files)){  
    if(types=="txt")output[[i]]<-read.table(paste0(folder,"/",files[i]),sep=";",row.names=1)
    if(types=="xlsx")output[[i]]<-read.xlsx(paste0(folder,"/",files[i]))
    if(types=="rwl")output[[i]]<-read.rwl(paste0(folder,"/",files[i]))
    if(types=="csv")output[[i]]<-read.table(paste0(folder,"/",files[i]),sep=separator,dec=decimal,row.names=1,header=T)
  }
  names(output)<-gsub(paste0(".",types),"",files) # naming list elements
  output
}

## ---------- unify.categories ####
## Unify catogiries inherited from metadata file
## Universal function working for all files
#
# input - inpud data.frame
# 
unify.categories<-function(input){
  
  ## Testing argiments
  # input=site.list
  
  composition<-rep("Mixed",nrow(input))
  composition[which(input$species_composition %in% c("Monoculture","MON","mon"))]<-"Monoculture"
  composition[which(input$species_composition == "NULL")]<-"No_data"
  
  management<-rep("Managed",nrow(input))
  management[which(input$management_type %in% c( "manpa","manp"))]<-"ManPa"
  management[which(input$management_type %in% c( "non-man","unm"))]<-"Unmanaged"
  management[which(input$management_type == "NULL")]<-"No_data"
  
  input$management_type<-management
  input$species_composition<-composition
  
  species<-input$species 
  species[which(species %in% c("QUPE","QURO","QUSP","QUsp"))]<-"QUSP"
  input$species<-species
  
  input<-input[which(input$species %in% c("ABAL", "FASY", "PCAB", "PISY", "QUSP")),]
  
  return(input)
}

## ---------- prepare.dataset ####
## Prepares a dataset for individual plots
## Each plot (Species/site according to the TAČR database) is one item of the list (output)
## Requires the reshape function
#
# site.list – file with metadata for individual sites
# min.trees – minimum number of trees required per plot
# min.out.year – dataset length constraint (if the series does not reach 2004 (or another specified year), it is not included)
#
prepare.dataset<-function(site.list,
                          tree.core,
                          min.trees=5){
  
  ## Testing arguments
  # site.list=site.list
  # tree.core=tree.core
  # min.trees=5
  # i="M000101QUSP"
  # i="P810312QURO_trw_data"
  # i="RTS1"
  # i="Certoryje"
  
  # výstupní list
  output.main<-list()
  
  # loop beginning, number of itteration correspond with number of sites in tree.core file
  for(i in unique(tree.core$site)){
    # data subset
    sub.data<-subset(tree.core,site==i)
    
    # check! - minimal number of trees per site
    if(length(unique(sub.data$tree))>=min.trees){
      
      # check! - if TRW is in wrong format, it is corrected
      if(max(sub.data$TRW,na.rm=T)>100) sub.data$TRW<-sub.data$TRW/100
      
      # Calculates mean TRW (if multiple cores is avaiable)
      output.site<-aggregate(TRW~year+tree,sub.data,mean,na.rm=T)
      
      # getting site_code
      site.name<-"some site"
      species.name<-"some species"
      ind<-which(site.list$site_code==i)
      if(length(ind)>0){
        site.name<-site.list$site_name[ind]
        species.name<-site.list$species[ind]
      }
      
      # output data.frame per site
      output.temp<-data.frame(Sitecode=as.character(unique(sub.data$site)), # Site name - inherited from database
                              Sitename=site.name, # original site name
                              Species=species.name, # species name
                              TreeID=output.site$tree, # tree ID
                              calendar.year=output.site$year, # callendar zear
                              cambial.age=NA, # cambial age
                              cum.DBH=NA, # cumulative DBH
                              cum.BA=NA, # cumilative BAI
                              TRW=output.site$TRW, # TRW
                              BAI=NA) # BAI
      
      # Calculating remaining data
      for(j in unique(output.temp$TreeID)){
        inds<-which(output.temp$TreeID==j & !is.na(output.temp$TRW)) # tree index in output data.frame
        
        for.bai.trw<-data.frame(TRW=output.temp$TRW[inds])
        for.bai.p.o<-data.frame(series=c("TRW"),po=0)
        output.temp$BAI[inds]<-bai.in(for.bai.trw,for.bai.p.o)$TRW
        
        output.temp$cambial.age[inds]<-c(1:length(inds)) # add cambial age
        output.temp$cum.DBH[inds]<-(cumsum(output.temp$TRW[inds])/10)*2 # add cumulative DBH (cm)
        output.temp$cum.BA[inds]<-round(cumsum(output.temp$BAI[inds])*0.01,2) # add cumulative BAI (cm^2)
      }
      output.main[[i]]<-output.temp
    }
  }
  return(output.main)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ------------------------------ Data detrending ------------------------------ ####
## ---------- detrend.qgam ####
## Performs qGAM detrending
#
# input - data.frame with TRW data
#
detrend.qgam <- function(input){
  # input=temp_trw
  
  output<-data.frame(matrix(NA,nrow=nrow(input),ncol=ncol(input)),row.names=rownames(input))
  names(output)<-names(input)
  
  for(i in 1:ncol(input)){
    
    # print(i)
    
    sub<-data.frame(year=as.numeric(as.character(rownames(input))),
                    trw=input[,i])
    sub<-na.omit(sub)
    sub$pred <- predict(lm(trw ~ year, data = sub))
    
    tryCatch({
      fit <- NULL
      invisible(
        suppressMessages(
          suppressWarnings(
            capture.output(
              fit <- qgam(trw ~ s(year), data = sub, qu = 0.60)
            )
          )
        )
      )
      
      sub$pred<-predict(fit)
      
    }, error = function(e) {
      
    })
    
    sub$rwi<-sub$trw/sub$pred
    
    output[which(as.numeric(as.character(rownames(output))) %in% sub$year),i]<-sub$rwi
  }
  
  return(output)
}

## ---------- detrend.gam ####
## Performs GAM detrending
#
# input - data.frame with TRW data
#
detrend.gam <- function(input){
  # input=temp_trw
  
  output<-data.frame(matrix(NA,nrow=nrow(input),ncol=ncol(input)),row.names=rownames(input))
  names(output)<-names(input)
  
  for(i in 1:ncol(input)){
    
    # print(i)
    
    sub<-data.frame(year=as.numeric(as.character(rownames(input))),
                    trw=input[,i])
    sub<-na.omit(sub)
    sub$pred <- predict(lm(trw ~ year, data = sub))
    
    tryCatch({
      fit <- NULL
      invisible(
        suppressMessages(
          suppressWarnings(
            capture.output(
              fit <- gam(trw ~ s(year), data = sub)
            )
          )
        )
      )
      
      sub$pred<-predict(fit)
      
    }, error = function(e) {
      
    })
    
    sub$rwi<-sub$trw/sub$pred
    
    output[which(as.numeric(as.character(rownames(output))) %in% sub$year),i]<-sub$rwi
  }
  
  return(output)
}

## ---------- create.rwi ####
## Detrend data in input file
#
# input.data - data.frame with TRW data
# detrend - function used for detrending (GAM - qGAM, spline or mean)
#
create.rwi<-function(input.data, detrend="mean"){
  
  # input.data=prepared.data.tree.means$P000230PCAB
  
  years<-as.numeric(as.character(unique(input.data$calendar.year)[order(unique(input.data$calendar.year))]))
  temp_trw<-data.frame(matrix(NA,
                              nrow=length(years),
                              ncol=length(unique(input.data$TreeID))))
  rownames(temp_trw)<-unique(input.data$calendar.year)[order(unique(input.data$calendar.year))]
  names(temp_trw)<-unique(input.data$TreeID)
  
  for(i in names(temp_trw)){
    sub<-subset(input.data,
                TreeID==i)
    year_index<-which(as.numeric(rownames(temp_trw)) %in% sub$calendar.year)
    temp_trw[year_index,i]<-sub$TRW
  }
  
  temp_trw <- temp_trw[,which(apply(temp_trw,2,function(x) length(na.omit(x)))>5)]
  
  ## Spline detrending
  if(detrend=="Spline"){
    
    binary_trw <- apply(temp_trw, 2, function(x) {
      # Find indices of non-NA values
      non_na_indices <- which(!is.na(x))
      
      if (length(non_na_indices) > 0) {
        # Set the range of non-NA values to 1
        x[min(non_na_indices):max(non_na_indices)] <- 1
      }
      
      return(x)
    })
    
    temp_trw <- temp_trw %>%
      mutate(across(everything(), ~ ifelse(is.na(.x) & cumsum(!is.na(.x)) > 0, 0, .x)))
    
    temp_trw <- temp_trw * binary_trw
    
    trw_rwi <- detrend(rwl = temp_trw, method = "Spline")
    
  }
  
  ## GAM detrending
  if(detrend=="GAM"){
    trw_rwi <- detrend.gam(temp_trw)
  }
  
  ## qGAM detrending
  if(detrend=="qGAM"){
    trw_rwi <- detrend.gam(temp_trw)
  }
  
  ## Mean detrending
  if(detrend=="mean"){
    trw_rwi <- temp_trw / colMeans(temp_trw, na.rm=T)
  }
  
  return(trw_rwi)
}

## ---------- calc.rwi.data ####
## Calculates detrended ring-width indices (RWI) for all sites in a dataset
#
# input - list of data.frames, each containing TRW data for one site
#         (typically the output of prepared.data.tree.means)
#
calc.rwi.data<-function(input){
  
  ## Testing arguments
  # input=prepared.data.tree.means
  
  output <- list()
  for(i in names(input)){
    # print(i)
    output[[i]]<-create.rwi(prepared.data.tree.means[[i]])
  }
  return(output)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ------------------------- Climate data calculations ------------------------- ####
## ---------- calculate.bal ####
## Calculates annual water balance statistics for a single site
#
# input - data.frame with monthly climate data (at least Year, Month, Temp, Prec)
# lat   - site latitude in decimal degrees (used for PET calculation via Thornthwaite)
# 
calculate.bal<-function(input,lat){
  
  ## Testing arguments
  # input=climate.data$V002425PCAB_clim
  # lat=subset(site_list,site_code=="V002425PCAB")$lat
  
  input$Pet <- thornthwaite(input$Temp, lat)
  input$Bal <- input$Prec - input$Pet
  
  input<-subset(input, Month>=first_month & Month<=last_month) # časové rozpětí pro výpoty
  
  out<-data.frame(year=unique(input$Year))
  
  out$min_bal<-aggregate(Bal~Year,data=input,min)$Bal
  out$mean_bal<-aggregate(Bal~Year,data=input,mean)$Bal
  
  return(out)
}

## ---------- calculate.spei ####
## Calculates annual SPEI-based drought indices for a single site
#
# input - data.frame with monthly climate data (at least Year, Month, Temp, Prec)
# lat   - site latitude in decimal degrees (used for PET calculation via Thornthwaite)
#
calculate.spei<-function(input,lat){
  
  ## Testing arguments
  # input=climate.data$V002425PCAB_clim
  # lat=subset(site_list,site_code=="V002425PCAB")$lat
  
  input$Pet <- thornthwaite(input$Temp, lat)
  input$Bal <- input$Prec - input$Pet
  input$spei1 <- spei(input[, "Bal"], 1)$fitted
  input$spei3 <- spei(input[, "Bal"], 3)$fitted
  
  input<-subset(input, Month>=first_month & Month<=last_month) # considerring only given months
  
  
  out<-data.frame(year=unique(input$Year))
  
  out$min_spei1<-aggregate(spei1~Year,data=input,min)$spei1
  out$min_spei3<-aggregate(spei3~Year,data=input,min)$spei3
  
  out$mean_spei1<-aggregate(spei1~Year,data=input,mean)$spei1
  out$mean_spei3<-aggregate(spei3~Year,data=input,mean)$spei3
  
  return(out)
}

## ---------- calculate.temp ####
## Calculates annual temperature statistics (growing season April–September)
#
# input - data.frame with monthly climate data (at least Year, Month, Temp)
#
calculate.temp<-function(input){
  
  ## Testing arguments
  # input=climate.data$V002425PCAB_clim
  
  input<-subset(input, Month>=first_month & Month<=last_month) # considerring only given months
  
  
  out<-data.frame(year=unique(input$Year))
  
  out$mean_temp<-aggregate(Temp~Year,data=input,mean)$Temp
  out$resid_temp<-aggregate(Temp~Year,data=input,mean)$Temp-mean(out$mean_temp)
  
  return(out)
}

## ---------- calculate.bals ####
## Calculates water balance statistics for all sites in a climate data list
#
# input - list of data.frames with monthly climate data for each site
# meta  - data.frame with site metadata (including columns site_code and lat)
#
calculate.bals<-function(input,meta){
  
  ## Testing arguments
  # input=climate.data
  # lat=site_list
  
  output<-list()
  
  for(i in names(input)){
    
    # print(i)
    
    latitude<-subset(site_list,site_code==substr(i,1,11))$lat
    
    if(length(latitude)==1){
      output[[i]]<-calculate.bal(climate.data[[i]], latitude)
      output[[i]]$site<-substr(i,1,11)
    } 
    if(is.null(latitude)){
      latitude<-50
      output[[i]]<-calculate.spei(climate.data[[i]], latitude)
      output[[i]]$site<-substr(i,1,11)
    } 
  }
  
  return(output)
}

## ---------- calculate.speis ####
## Calculates SPEI-based drought indices for all sites in a climate data list
#
# input - list of data.frames with monthly climate data for each site
# meta  - data.frame with site metadata (including columns site_code and lat)
#
calculate.speis<-function(input,meta){
  
  ## Testing arguments
  # input=climate.data
  # lat=site_list
  
  output<-list()
  
  for(i in names(input)){
    
    # print(i)
    
    latitude<-subset(site_list,site_code==substr(i,1,11))$lat
    
    if(length(latitude)==1){
      output[[i]]<-calculate.spei(climate.data[[i]], latitude)
      output[[i]]$site<-substr(i,1,11)
    } 
    if(is.null(latitude)){
      latitude<-50
      output[[i]]<-calculate.spei(climate.data[[i]], latitude)
      output[[i]]$site<-substr(i,1,11)
    } 
  }
  
  return(output)
}

## ---------- calculate.temps ####
## Calculates growing-season temperature statistics for all sites in a climate data list
#
# input - list of data.frames with monthly climate data for each site
#
calculate.temps<-function(input){
  
  ## Testing arguments
  # input=climate.data
  
  output<-list()
  
  for(i in names(input)){
    # print(i)
    
    output[[i]]<-calculate.temp(climate.data[[i]])
    output[[i]]$site<-substr(i,1,11)
    
  }
  
  return(output)
}

## ---------- assign.speis ####
## Assigns SPEI-based drought indices to the main dataset
#
# input - data.frame with core data (must contain columns: site_code, year)
# spei  - data.frame with annual SPEI indices (must contain columns:
#         site, year, mean_spei1; optionally also min_spei1, min_spei3, mean_spei3)
#
assign.speis<-function(input, spei){
  
  ## Testing arguments
  # input=df.all
  # spei=speis
      
  input$mean_spei1<-NA
  
  for(i in unique(input$site_code)){
    
    # print(i)
    
    indexes<-which(input$site_code==i & input$year %in% spei$year)
    
    sub.speis<-subset(speis, site==i)
    
    if(nrow(sub.speis)>0){
      
      sub.speis<-sub.speis[which(sub.speis$year %in% input$year[indexes]),]
          
      input$mean_spei1[indexes]<-sub.speis$mean_spei1  
      
    }
    
  }
  
  # output<-na.omit(input)
  output<-input
  
  return(output)
}
## ---------- assign.temps ####
## Assigns growing-season temperature indices to the main dataset
#
# input - data.frame with core data (must contain columns: site_code, year)
# temps - data.frame with annual temperature statistics (must contain columns:
#         site, year, mean_temp; optionally resid_temp)
#
assign.temps<-function(input, temps){
  
  ## Testing arguments
  # input=df.all
  # temps=temps
  
  input$mean_temp<-NA    
  
  for(i in unique(input$site_code)){
    
    # print(i)
    
    indexes<-which(input$site_code==i & input$year %in% temps$year)
    
    sub.temps<-subset(temps, site==i)
    
    if(nrow(sub.temps)>0){
      
      sub.temps<-sub.temps[which(sub.temps$year %in% input$year[indexes]),]
      
      input$mean_temp[indexes]<-sub.temps$mean_temp         
      
    }
    
  }
  
  # output<-na.omit(input)
  output<-input
  
  return(output)
}
## ---------- assign.cwb ####
## Assigns climatic water balance (CWB) indices to the main dataset
#
# input - data.frame with core data (must contain columns: site_code, year)
# temps - (not used directly) – climatic water balance is taken from the object 'bals'
#         which must be available in the workspace and contain columns:
#         site, year, mean_bal
#
assign.cwb<-function(input, temps){
  
  ## Testing arguments
  # input=df.all
  # bals=bals
  
  input$cwb<-NA      
  
  for(i in unique(input$site_code)){
    
    # print(i)
    
    indexes<-which(input$site_code==i & input$year %in% temps$year)
    
    sub.bals<-subset(bals, site==i)
    
    if(nrow(sub.bals)>0){
      
      sub.bals<-sub.bals[which(sub.bals$year %in% input$year[indexes]),]
      
      input$cwb[indexes]<-sub.bals$mean_bal       
      
    }
    
  }
  
  # output<-na.omit(input)
  output<-input
  
  return(output)
}
## ---------- assign.resid ####
## Assigns residual temperature (anomalies) to the main dataset
#
# input - data.frame with core data (must contain columns: site_code, year)
# temps - data.frame with annual temperature statistics (must contain columns:
#         site, year, resid_temp)
#
assign.resid<-function(input, temps){
  
  ## Testing arguments
  # input=df.all
  # temps=temps
  
  input$resid_temp<-NA      
  
  for(i in unique(input$site_code)){
    
    # print(i)
    
    indexes<-which(input$site_code==i & input$year %in% temps$year)
    
    sub.temps<-subset(temps, site==i)
    
    if(nrow(sub.temps)>0){
      
      sub.temps<-sub.temps[which(sub.temps$year %in% input$year[indexes]),]
      
      input$resid_temp[indexes]<-sub.temps$resid_temp         
      
    }
    
  }
  
  # output<-na.omit(input)
  output<-input
  
  return(output)
}

## ---------- add.sajrajty  ####
## Adds site- and year-specific external data (e.g. SOx, Ndep) to the modelling dataset
#
# dta_mod  - data.frame with core data used for modelling
#            (must contain columns: site_code, year; plus the target column given by col_name)
# sajrajty - data.frame in wide format with one row per site and multiple year columns
#            (must contain column site_code and yearly values in columns 6:123)
# col_name - name of the column in dta_mod where the extracted values will be stored
#            (e.g. "sox" or another variable name)
#
add.sajrajty<-function(dta_mod, sajrajty, col_name){
  
  ## Testing argumets
  # dta_mod=mod_dataset
  # sajrajty=extracted.sox
  # col_name="sox"
  
  start_index<-4
  if(col_name=="sox") start_index<-3
  
  for(i in unique(sajrajty$site_code)){
    # i=sajrajty$site_code[1]
    sub_sajrajty<-data.frame(site_code=i,
                             year=as.numeric(as.character(substr(names(sajrajty)[6:123],start_index,10))),
                             sajrajt=as.numeric(sajrajty[which(sajrajty$site_code==i),6:123]))
    
    sub_sajrajty<-sub_sajrajty[which(paste0(sub_sajrajty$site_code,"_",sub_sajrajty$year) %in% paste0(dta_mod$site_code,"_",dta_mod$year)),]
    
    dta_mod[which(paste0(dta_mod$site_code,"_",dta_mod$year) %in% paste0(sub_sajrajty$site_code,"_",sub_sajrajty$year)),col_name]<-sub_sajrajty$sajrajt
  }
  
  return(dta_mod)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ------------------------- Chronologies calculations ------------------------- ####
## ---------- create.chronology ####
## Creates a site-level ring-width chronology from RWI data
#
# input.data - data.frame with detrended RWI values for a single site
#              (columns = tree series, rows = calendar years)
#
create.chronology<-function(input.data){
  
  ## Testing arguments
  # input.data=rwi.data[[unique(output$site_code)]]
  
  trw_rwi <- input.data
  
  trw_crn <- chron(trw_rwi, biweight = TRUE, prewhiten = TRUE, nyrs=30)
  return(trw_crn)
}

## ---------- create.chronologies ####
## Creates ring-width chronologies for all sites in a list of RWI datasets
#
# input.data - list of data.frames, each containing detrended RWI values
#              (typically output from calc.rwi.data)
#
create.chronologies<-function(input.data){
  
  ## Testing arguments
  # input.data=rwi.data[[unique(output$site_code)]]
  
  output <- list()
  
  for(i in names(input.data)){
    output[[i]] <- create.chronology(input.data[[i]])
  }
  
  input.data
  
  return(output)
}

## ---------- bootstrap.chronologies ####
## Estimates inter-annual variability of standard and residual chronologies
## using non-parametric bootstrap across trees
#
# input  - data.frame with annual chronology values for individual trees
#          (must contain columns: species, year, std, res)
# sp     - species code to be analysed (character; matches input$species)
# n_boot - number of bootstrap replicates (default = 1000)
#
bootstrap.chronologies<-function(input, sp, n_boot=1000){
  
  ## Testing arguments
  # input<-dataset
  # sp<-"ABAL"
  # n_boot=1000
  
  input<-subset(input, species==sp)
  
  out<-data.frame(species=sp,
                  year=c(1961:2020),
                  std.mean=NA, std.min=NA, std.mid=NA, std.max=NA,
                  res.mean=NA, res.min=NA, res.mid=NA, res.max=NA,
                  
                  sd.std.mean=NA, sd.std.min=NA, sd.std.mid=NA, sd.std.max=NA, 
                  sd.res.mean=NA, sd.res.min=NA, sd.res.mid=NA, sd.res.max=NA,
                  
                  cv.std.mean=NA, cv.std.min=NA, cv.std.mid=NA, cv.std.max=NA, 
                  cv.res.mean=NA, cv.res.min=NA, cv.res.mid=NA, cv.res.max=NA)
  
  out$std.mean<-merge(out, aggregate(std~year,input,mean), by.x="year", by.y="year")$std
  out$res.mean<-merge(out, aggregate(res~year,input,mean), by.x="year", by.y="year")$res
  
  out$sd.std.mean<-merge(out, aggregate(std~year,subset(input, species==sp & year>=1961 & year<=2020),FUN=sd), by.x="year", by.y="year")$std
  out$sd.res.mean<-merge(out, aggregate(res~year,subset(input, species==sp & year>=1961 & year<=2020),FUN=sd), by.x="year", by.y="year")$res
  
  out$cv.std.mean<-merge(out, aggregate(std~year,subset(input, species==sp & year>=1961 & year<=2020),FUN=function(x){sd(x)/mean(x)}), by.x="year", by.y="year")$std
  out$cv.res.mean<-merge(out, aggregate(res~year,subset(input, species==sp & year>=1961 & year<=2020),FUN=function(x){sd(x)/mean(x)}), by.x="year", by.y="year")$res
  
  for(i in 1:nrow(out)){
    
    sub.data<-subset(input,year==out$year[i])
    sub.data<-na.omit(sub.data)
    
    if(nrow(sub.data)>=5){
      cor_boot_std <- numeric(n_boot) 
      cor_boot_res <- numeric(n_boot) 
      cor_boot_sd_std <- numeric(n_boot) 
      cor_boot_sd_res <- numeric(n_boot) 
      cor_boot_cv_std <- numeric(n_boot) 
      cor_boot_cv_res <- numeric(n_boot) 
      for (j in 1:n_boot) {
        # Sample rows with replacement
        sample_indices <- sample(1:nrow(sub.data), size = nrow(sub.data), replace = TRUE)
        sample_data <- sub.data[sample_indices, ]
        
        # Calculate correlation for the sample
        cor_boot_std[j] <- mean(sample_data$std)
        cor_boot_res[j] <- mean(sample_data$res)
        cor_boot_sd_std[j] <- sd(sample_data$std)
        cor_boot_sd_res[j] <- sd(sample_data$res)
        cor_boot_cv_std[j] <- sd(sample_data$std)/mean(sample_data$std)
        cor_boot_cv_res[j] <- sd(sample_data$res)/mean(sample_data$res)
      }
      
      out[i,c("std.min","std.mid","std.max")]<-quantile(cor_boot_std, c(0.025, 0.500, 0.975))
      out[i,c("res.min","res.mid","res.max")]<-quantile(cor_boot_res, c(0.025, 0.500, 0.975))
      out[i,c("sd.std.min","sd.std.mid","sd.std.max")]<-quantile(cor_boot_sd_std, c(0.025, 0.500, 0.975))
      out[i,c("sd.res.min","sd.res.mid","sd.res.max")]<-quantile(cor_boot_sd_res, c(0.025, 0.500, 0.975))
      out[i,c("cv.std.min","cv.std.mid","cv.std.max")]<-quantile(cor_boot_cv_std, c(0.025, 0.500, 0.975))
      out[i,c("cv.res.min","cv.res.mid","cv.res.max")]<-quantile(cor_boot_cv_res, c(0.025, 0.500, 0.975))
    }
  }
  
  return(out)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ----------------------- Setting up final data.frames ------------------------ ####
## ---------- init.data.frame ####
## Initializes the output data.frame for a single site and adds metadata
## together with basic stand-level growth characteristics
#
# input.data - data.frame with tree-level data for one site
#              (must contain at least: calendar.year, TreeID, cambial.age,
#               cum.DBH, cum.BA)
# metadata   - data.frame (single row) with site metadata, including:
#              site_code, site_name, management_type, species_composition,
#              site_category, elevation, lon, lat, species
#
init.data.frame<-function(input.data,
                          metadata,
                          years=c(1960:2020)){
  
  ## Testing data.frame
  # input.data=prepared.data.tree.means$P810114PCAB
  # metadata=subset(site.list,
  #                 site_code=="P810114PCAB")
  
  new.df<-data.frame(site_code=metadata$site_code,
                     site_name=metadata$site_name,
                     management_type=metadata$management_type,
                     species_composition=metadata$species_composition,
                     site_category=metadata$site_category,
                     elevation=metadata$elevation,
                     x=metadata$lon,
                     y=metadata$lat,
                     species=metadata$species,
                     
                     year=years,
                     sample_depth=NA,
                     cambial_age=NA,
                     max_cambial_age=NA,
                     cambial_age_range=NA,
                     cum_DBH=NA,
                     cum_BA=NA,
                     
                     cv_TRW=NA,
                     mid_TRW=NA,
                     std_chron=NA,
                     res_chron=NA,
                     last1=NA,
                     last3=NA,
                     last5=NA,
                     
                     sd_RWI=NA,
                     cv_RWI=NA,
                     last1_RWI=NA,
                     last3_RWI=NA,
                     last5_RWI=NA
  )
  
  temp<-aggregate(TreeID~calendar.year,input.data,length)
  new.df$sample_depth[which(new.df$year %in% temp$calendar.year)]<-temp$TreeID[which(temp$calendar.year %in% years)]
  
  temp<-aggregate(cambial.age~calendar.year,input.data,mean)
  new.df$cambial_age[which(new.df$year %in% temp$calendar.year)]<-round(temp$cambial.age[which(temp$calendar.year %in% years)],0)
  
  temp<-aggregate(cambial.age~calendar.year,input.data,function(x){quantile(x,probs=c(0.95))})
  new.df$max_cambial_age[which(new.df$year %in% temp$calendar.year)]<-round(temp$cambial.age[which(temp$calendar.year %in% years)],0)
  
  
  ## Setting calculations for age range
  temp<-aggregate(cambial.age~calendar.year,input.data,function(x){quantile(x,probs=c(0.95))})
  temp$cambial.age<-temp$cambial.age-aggregate(cambial.age~calendar.year,input.data,function(x){quantile(x,probs=c(0.05))})$cambial.age
  new.df$cambial_age_range[which(new.df$year %in% temp$calendar.year)]<-round(temp$cambial.age[which(temp$calendar.year %in% years)],0)
  
  temp<-aggregate(cum.DBH~calendar.year,input.data,mean)
  new.df$cum_DBH[which(new.df$year %in% temp$calendar.year)]<-round(temp$cum.DBH[which(temp$calendar.year %in% years)],1)
  
  temp<-aggregate(cum.BA~calendar.year,input.data,mean)
  new.df$cum_BA[which(new.df$year %in% temp$calendar.year)]<-round(temp$cum.BA[which(temp$calendar.year %in% years)],1)

  return(new.df)
}

## ---------- calculate.site.data ####
## Prepare annual stand-level data for a single site
#
# input    - data.frame with tree-level data for one site
#            (must contain at least: calendar.year, TreeID, TRW)
# metadata - data.frame (single row) with site metadata
#            (passed to init.data.frame; see init.data.frame for details)
#
calculate.site.data<-function(input,
                              metadata){
  
  ## Testing arguments
  # input=prepared.data.tree.means$C004104FASY
  # metadata=subset(site.list,
  #                 site_code=="C004104FASY")
  
  output<-init.data.frame(input,metadata)
  rwi_data<-rwi.data[[unique(output$site_code)]]
  chronology<-crn.data[[unique(output$site_code)]]
  
  for(j in c(output$year)){
    # print(j)
    sub.crn<-na.omit(chronology[which(as.numeric(as.character(rownames(chronology)))==j),])
    sub.trw<-na.omit(subset(input,calendar.year==j))
    sub.rwi<-na.omit(as.numeric(rwi_data[which(as.numeric(rownames(rwi_data))==j),]))
    
    if(nrow(sub.trw)>=5 & nrow(sub.crn)==1 & length(sub.rwi)>=5){
      
      sub.trw.1<-na.omit(subset(input,calendar.year==(j-1)))
      sub.trw.3<-na.omit(subset(input,calendar.year>=(j-3) & calendar.year<(j)))
      sub.trw.5<-na.omit(subset(input,calendar.year>=(j-5) & calendar.year<(j)))
      
      sub.rwi.1<-na.omit(as.numeric(rwi_data[which(as.numeric(rownames(rwi_data))==(j-1)),]))
      sub.rwi.3<-na.omit(as.numeric(colMeans(rwi_data[which(as.numeric(rownames(rwi_data))>=(j-3) & as.numeric(rownames(rwi_data))<(j)),])))
      sub.rwi.5<-na.omit(as.numeric(colMeans(rwi_data[which(as.numeric(rownames(rwi_data))>=(j-5) & as.numeric(rownames(rwi_data))<(j)),])))
      
      output[which(output$year==j),c("std_chron")]<-sub.crn$std
      output[which(output$year==j),c("res_chron")]<-sub.crn$res
      
      
      output[which(output$year==j),c("cv_TRW")]<-sd(sub.trw$TRW)/mean(sub.trw$TRW)
      output[which(output$year==j),c("mid_TRW")]<-mean(sub.trw$TRW)
      
      output[which(output$year==j),c("sd_RWI")]<-sd(sub.rwi)
      output[which(output$year==j),c("cv_RWI")]<-sd(sub.rwi)/mean(sub.rwi)
      
      if(length(unique(sub.trw.1$calendar.year))==1) output[which(output$year==j),c("last1")]<-mean(sub.trw.1$TRW)
      if(length(unique(sub.trw.3$calendar.year))==3) output[which(output$year==j),c("last3")]<-mean(aggregate(TRW~TreeID,data=sub.trw.3,FUN=mean)$TRW)
      if(length(unique(sub.trw.5$calendar.year))==5) output[which(output$year==j),c("last5")]<-mean(aggregate(TRW~TreeID,data=sub.trw.5,FUN=mean)$TRW)

      if(length(unique(sub.trw.1$calendar.year))==1) output[which(output$year==j),c("last1_RWI")]<-mean(sub.rwi.1)
      if(length(unique(sub.trw.3$calendar.year))==3) output[which(output$year==j),c("last3_RWI")]<-mean(sub.rwi.3)
      if(length(unique(sub.trw.5$calendar.year))==5) output[which(output$year==j),c("last5_RWI")]<-mean(sub.rwi.5)
    }
  }
  return(output)
}

## ---------- calculate.all.sites.data ####
## Calculates stand-level growth and chronology metrics for all sites
#
# input    - list of data.frames with tree-level data for all sites
#            (each element corresponds to one site; e.g. prepared.data.tree.means)
# metadata - data.frame with site metadata for all sites
#            (must contain at least: site_code and other columns required by
#             calculate.site.data.20ths / init.data.frame)
#
calculate.all.sites.data<-function(input,
                                   metadata){
  
  ## Testing arguments
  # input=prepared.data.tree.means
  # metadata=site.list
  # i="P810114PCAB"
  
  complete.output<-list()
  for(i in names(input)){
    
    # print(i)
    
    complete.output[[i]]<-calculate.site.data(prepared.data.tree.means[[i]],
                                              subset(site.list,
                                                     site_code==i))
  }
  return(complete.output)
}

## ---------- prepare.data ####
## Filters and prepares tree-ring datasets based on minimum age, sample depth,  
## and required temporal coverage
#
# input  - a *list* of data.frames, one per site (e.g. raw_all)
#          Each data.frame must contain at least the columns:
#          - year           ... calendar year
#          - cambial_age    ... age of cambium for each observation
#          - sample_depth   ... number of trees available in the given year
#
# min_age           - minimum required cambial age for an observation to be retained  
#                     (default = 20 years)
# min_trees         - minimum number of trees in a given year required for that year/site  
#                     to be included (default = 5)
# min_final_year    - site must have data reaching at least this year (default = 2015)
# min_starting_year - site must have data beginning at or before this year (default = 2015)
#
prepare.data<-function(input,
                       min_age=20,
                       min_trees=5,
                       min_final_year=2015,
                       min_starting_year=2015){
  
  ## Testing arguments
  # input=raw_all
  # min_age=20
  # min_trees=5
  
  output<-list()
  for(i in names(input)){
    
    sub<-na.omit(input[[i]])
    
    sub<-sub[which(sub$cambial_age>=min_age & sub$sample_depth>=min_trees),]
    
    if(max(sub$year)>=min_final_year & min(sub$year)<=min_starting_year) output[[i]]<-sub
    
  }
  output<-do.call(rbind,output)
  return(output)
}
## ---------------------------------------------------------------- scale.variables ####
## Scales numerical predictor variables (centering and standardizing) for modelling
#
# input - data.frame containing model predictors; must include the variables listed below.
#
scale.variables<-function(input){
  
  ## Testing arguments
  # input=mod_dataset
  
  input$mid_TRW<-scale(input$mid_TRW)
  input$elevation<-scale(input$elevation)
  input$max_cambial_age<-scale(input$max_cambial_age)
  input$max_age<-scale(input$max_age)
  input$median_age<-scale(input$median_age)
  input$cambial_age_range<-scale(input$cambial_age_range)
  input$max_range<-scale(input$max_range)
  input$median_range<-scale(input$median_range)
  input$cambial_age<-scale(input$cambial_age)
  input$mean_temp<-scale(input$mean_temp)
  input$normal_temp<-scale(input$normal_temp)
  input$resid_temp<-scale(input$resid_temp)
  input$mean_cwb<-scale(input$mean_cwb)
  input$normal_cwb<-scale(input$normal_cwb)
  input$resid_cwb<-scale(input$resid_cwb)
  input$last5<-scale(input$last5)
  
  input$sox<-scale(input$sox)
  input$sothf<-scale(input$sothf)
  input$nox<-scale(input$nox)
  input$nh<-scale(input$nh)
  input$nsum<-scale(input$nsum)
  
  input$pH_FH<-scale(input$pH_FH)
  input$pH_L1<-scale(input$pH_L1)
  input$C.N_FH<-scale(input$C.N_FH)
  input$C.N_L1<-scale(input$C.N_L1)
  
  return(input)
}


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ---------------------- Calculates site characteristics ---------------------- ####
## ---------- calculate.site.statistics ####
## Calculates basic time-series statistics for a single site
#
# input - data.frame with tree-level ring-width data for one site
#         (must contain columns: Sitecode, Species, TreeID, calendar.year)
#
calculate.site.statistics <- function(input){
  
  ## Testing arguments
  # input=input$P810114PCAB
  
  site <- unique(input$Sitecode)
  species <- unique(input$Species)
  number_of_trees <- length(unique(input$TreeID))
  number_of_trees_since61 <- length(unique(subset(input, calendar.year>1960)$TreeID))
  start_years <- aggregate(calendar.year ~ TreeID, input, min)
  end_years <- aggregate(calendar.year ~ TreeID, input, max)
  series_lengths <- aggregate(calendar.year ~ TreeID, input, length)
  
  start_years<-start_years[order(start_years$calendar.year),]
  end_years<-end_years[order(end_years$calendar.year, decreasing = T),]
  
  output <- data.frame(site_code = site,
                       species = species,
                       n_total = number_of_trees,
                       n_since61 = number_of_trees_since61,
                       start_oldest = min(start_years$calendar.year),
                       start_5 = start_years$calendar.year[5],
                       start_25perc = round(quantile(start_years$calendar.year, probs=0.25), 0), #zajisti min 1 strom
                       start_median = round(median(start_years$calendar.year), 0),
                       end_max = max(end_years$calendar.year),
                       end_5 =  end_years$calendar.year[5],
                       end_25perc = round(quantile(end_years$calendar.year, probs=0.75), 0),
                       end_median = round(median(end_years$calendar.year), 0))
  
  return(output)
}

## ---------- calculate.site.statistics.allsites ####
## Calculates basic time-series statistics for all sites in a list
#
# input - list of data.frames with tree-level ring-width data for multiple sites
#         (e.g. prepared.data.tree.means; each element = one site)
# min_n - minimum number of trees required for a site to be included (default = 5)
#
calculate.site.statistics.allsites <- function(input, min_n = 5){
  
  ## Testing arguments
  # input = prepared.data.tree.means
  # min_n = 5
  
  output <- data.frame(site_code = character(),
                       species = character(),
                       n_total = numeric(),
                       n_since61 = numeric(),
                       start_oldest = numeric(),
                       start_5 = numeric(),
                       start_25perc = numeric(), 
                       start_median = numeric(),
                       end_max = numeric(),
                       end_5 =  numeric(),
                       end_25perc = numeric(),
                       end_median = numeric())
  
  for(i in names(prepared.data.tree.means)){
    sub_data <- prepared.data.tree.means[[i]]
    
    if(length(unique(sub_data$TreeID)) >= min_n){
      output <- rbind(output, calculate.site.statistics(sub_data))
    }
    
  }
  
  return(output)
}






