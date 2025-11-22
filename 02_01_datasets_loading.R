## Funkce
source("01_functions.R")

## Loading common datasets ####
site.list<-read.table("Input_data/Site_meta/sites_all.csv",sep=",",header=T)
site.list<-unify.categories(site.list)

## Tree.core list ####
tree.core<-read.table("Input_data/Database_file/sampleTable.csv",sep=",", header=T)

tree.core<-data.frame(site=tree.core[,1],
                      species=substr(tree.core[,1],8,11),
                      tree=formatC(tree.core[,2],width=4,flag="0"),
                      core=tree.core[,3],
                      year=tree.core[,4],
                      TRW=as.numeric(as.character(tree.core[,5]))/100)

## Mean TRW data ####
prepared.data.tree.means<-prepare.dataset(site.list, tree.core)

## Site summary ####
source("Data_preparation/10_site_summary.R")
# site_stats <- subset(site_stats, start_25perc<=1950 & end_25perc>=2010)
# site_stats <- site_stats[which(site_stats$species %in% c("ABAL", "FASY", "PCAB", "PISY", "QUSP")),]
# 
# site.list <- site.list[which(site.list$site_code %in% site_stats$site_code),]
# prepared.data.tree.means <- prepared.data.tree.means[site_stats$site_code]

# Loading chronology data
df.crn.all<-read.table("Calculated_datasets/Bootstrapped_crn_variability/chronologies_individual.txt",sep=";",dec=".",header=T)
df.crn<-read.table("Calculated_datasets/Bootstrapped_crn_variability/chronologies_data.txt",sep=";",dec=".",header=T)

# Loading site data
df.all<-read.table("Calculated_datasets/Finalized_datasets/dataset.txt",sep=";",dec=".",header=T)
df.all<-unify.categories(df.all)

## Saving sites suitable for the study
considered_sites<-unique(df.all$site_code)

## Site climate data ####
files<-list.files("Calculated_datasets/Site_climate")
climadata<-list()
for(i in files){
  climadata[[gsub("_clim.csv","",i)]]<-read.table(paste0("Calculated_datasets/Site_climate/",i),sep=",",header=T)
}

clim<-data.frame(site=names(climadata), MAT=NA, GsT=NA, MAP=NA, GsP=NA, elevation=NA, species=NA)

for(i in clim$site){
  
  elev<-site.list$elevation[which(site.list$site_code==i)]
  species<-site.list$species [which(site.list$site_code==i)]
  if(length(elev)==1) clim$elevation[which(clim$site==i)]<-site.list$elevation[which(site.list$site_code==i)]
  if(length(elev)==1) clim$species[which(clim$site==i)]<-site.list$species[which(site.list$site_code==i)]
  
  sub.dta<-climadata[[i]]
  clim$MAT[which(clim$site==i)]<-mean(aggregate(Temp~Year,mean,data=sub.dta)$Temp)
  clim$MAP[which(clim$site==i)]<-mean(aggregate(Prec~Year,sum,data=sub.dta)$Prec)
  
  
  sub.dta<-sub.dta[which(sub.dta$Month %in% c(4:9)),]
  clim$GsT[which(clim$site==i)]<-mean(aggregate(Temp~Year,mean,data=sub.dta)$Temp)
  clim$GsP[which(clim$site==i)]<-mean(aggregate(Prec~Year,sum,data=sub.dta)$Prec)
}

clim$species[which(clim$species %in% c("QUSP","QURO","QUPE","QUsp"))]<-"QUSP"
clim<-clim[which(clim$species %in% c("QUSP", "FASY", "PISY", "PCAB", "ABAL")),]
clim<-na.omit(clim)
considered_sites<-unique(df.all$site_code)
clim<-clim[which(clim$site %in% considered_sites),]

## Loading and adjusting climate.dataset
bals<-read.table("Calculated_datasets/Recalculated_climate/bals.txt", dec=".", sep=";")
speis<-read.table("Calculated_datasets/Recalculated_climate/speis.txt", dec=".", sep=";")
temps<-read.table("Calculated_datasets/Recalculated_climate/temps.txt", dec=".", sep=";")

clim.dataset<-subset(df.all,year>1960 & year<2018)
clim.dataset<-assign.cwb(clim.dataset, bals)
clim.dataset<-assign.speis(clim.dataset, speis)
clim.dataset<-assign.temps(clim.dataset, temps)
clim.dataset<-assign.resid(clim.dataset, temps)

## Loading model datasets ####
main_mod_dataset<-read.table("Calculated_datasets/Finalized_datasets/mod_dataset_main.txt", sep=";", dec=".", header=T)
abal_mod_dataset<-read.table("Calculated_datasets/Finalized_datasets/mod_dataset_ABAL.txt", sep=";", dec=".", header=T)
fasy_mod_dataset<-read.table("Calculated_datasets/Finalized_datasets/mod_dataset_FASY.txt", sep=";", dec=".", header=T)
pcab_mod_dataset<-read.table("Calculated_datasets/Finalized_datasets/mod_dataset_PCAB.txt", sep=";", dec=".", header=T)
pisy_mod_dataset<-read.table("Calculated_datasets/Finalized_datasets/mod_dataset_PISY.txt", sep=";", dec=".", header=T)
qusp_mod_dataset<-read.table("Calculated_datasets/Finalized_datasets/mod_dataset_QUSP.txt", sep=";", dec=".", header=T)

main_mod_dataset_scaled<-read.table("Calculated_datasets/Finalized_datasets/scaled_mod_dataset_main.txt", sep=";", dec=".", header=T)
abal_mod_dataset_scaled<-read.table("Calculated_datasets/Finalized_datasets/scaled_mod_dataset_ABAL.txt", sep=";", dec=".", header=T)
fasy_mod_dataset_scaled<-read.table("Calculated_datasets/Finalized_datasets/scaled_mod_dataset_FASY.txt", sep=";", dec=".", header=T)
pcab_mod_dataset_scaled<-read.table("Calculated_datasets/Finalized_datasets/scaled_mod_dataset_PCAB.txt", sep=";", dec=".", header=T)
pisy_mod_dataset_scaled<-read.table("Calculated_datasets/Finalized_datasets/scaled_mod_dataset_PISY.txt", sep=";", dec=".", header=T)
qusp_mod_dataset_scaled<-read.table("Calculated_datasets/Finalized_datasets/scaled_mod_dataset_QUSP.txt", sep=";", dec=".", header=T)

## Loading calculated models ####
# NOTE:
# The fitted models used in the manuscript are not included in the repository.
# They cannot be estimated on the reduced demo dataset and therefore are not loaded here.

# mod_main <- readRDS("Calculated_datasets/Calculated_models/model_MAIN.rds")
# mod_abal <- readRDS("Calculated_datasets/Calculated_models/model_ABAL.rds")
# mod_fasy <- readRDS("Calculated_datasets/Calculated_models/model_FASY.rds")
# mod_pcab <- readRDS("Calculated_datasets/Calculated_models/model_PCAB.rds")
# mod_pisy <- readRDS("Calculated_datasets/Calculated_models/model_PISY.rds")
# mod_qusp <- readRDS("Calculated_datasets/Calculated_models/model_QUSP.rds")

## Model formula ####
form_main <- log_cv_RWI ~ 
  
  species * (median_age + 
               median_range + 
               mid_TRW + 
               mean_temp + 
               mean_cwb + 
               sox + 
               pH_L1 +
               C.N_FH
  ) +
  species * (median_age:mean_temp +
               median_age:mean_cwb +
               median_range:mean_temp +
               median_range:mean_cwb +
               pH_L1 : mid_TRW +
               pH_L1 : mean_temp +
               sox : mid_TRW +
               sox : mean_temp)

form_species <- log_cv_RWI ~ 
  
  median_age + 
  median_range + 
  mid_TRW + 
  mean_temp + 
  mean_cwb + 
  sox + 
  pH_L1 +
  C.N_FH +
  median_age:mean_temp +
  median_age:mean_cwb +
  median_range:mean_temp +
  median_range:mean_cwb +
  pH_L1 : mid_TRW +
  pH_L1 : mean_temp +
  sox : mid_TRW +
  sox : mean_temp



