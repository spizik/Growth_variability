## load new data ####
df.all<-read.table("Calculated_datasets/Finalized_datasets/dataset.txt",sep=";",dec=".",header=T)
df.all<-unify.categories(df.all)

extracted.sox<-read.xlsx("Calculated_datasets/Depozitions/extracted_SO.xlsx")
extracted.sothf<-read.xlsx("Calculated_datasets/Depozitions/extracted_SOthf.xlsx")
extracted.nox<-read.xlsx("Calculated_datasets/Depozitions/extracted_NO.xlsx")
extracted.nh<-read.xlsx("Calculated_datasets/Depozitions/extracted_NH.xlsx")
extracted.nsum<-read.xlsx("Calculated_datasets/Depozitions/extracted_Nsum.xlsx")

soil_dataset_FH<-read.table("Calculated_datasets/Soil_data/soil_FH.txt", sep=";", header=T)
soil_dataset_L1<-read.table("Calculated_datasets/Soil_data/soil_lay1.txt", sep=";", header=T)
# soil_dataset<-read.table("Calculated_datasets/Soil_data/soil_FH_and_lay1.txt", sep=";", header=T)

bals<-read.table("Calculated_datasets/Recalculated_climate/bals.txt", dec=".", sep=";")
speis<-read.table("Calculated_datasets/Recalculated_climate/speis.txt", dec=".", sep=";")
temps<-read.table("Calculated_datasets/Recalculated_climate/temps.txt", dec=".", sep=";")

# climate.dataset ####
clim.dataset<-subset(df.all,year>1960 & year<2018)
clim.dataset<-assign.cwb(clim.dataset, bals)
clim.dataset<-assign.speis(clim.dataset, speis)
clim.dataset<-assign.temps(clim.dataset, temps)
clim.dataset<-assign.resid(clim.dataset, temps)

## Basic dataset ####
mod_dataset<-clim.dataset[,c("site_code", 
                             "species", 
                             "elevation",
                             "x",
                             "y",
                             "management_type", 
                             "site_category",       
                             "sample_depth", 
                             "year", 
                             "max_cambial_age", 
                             "cambial_age_range", 
                             "cambial_age",
                             "cv_TRW",
                             "sd_RWI",
                             "cv_RWI",
                             "mid_TRW", 
                             "mean_temp",
                             "resid_temp",
                             "cwb",
                             "last5")]

## spočítá dlouhodobé průměry Teplot a CW ####
mod_dataset<-merge(mod_dataset,aggregate(mean_temp~site_code,mod_dataset,mean),by.x="site_code",by.y="site_code")
mod_dataset<-merge(mod_dataset,aggregate(cwb~site_code,mod_dataset,mean),by.x="site_code",by.y="site_code")

mod_dataset<-merge(mod_dataset,aggregate(max_cambial_age~site_code,mod_dataset,max),by.x="site_code",by.y="site_code")
mod_dataset<-merge(mod_dataset,aggregate(cambial_age_range~site_code,mod_dataset,max),by.x="site_code",by.y="site_code")

mod_dataset<-merge(mod_dataset,aggregate(cambial_age~site_code,mod_dataset,median),by.x="site_code",by.y="site_code")
mod_dataset<-merge(mod_dataset,aggregate(cambial_age_range.x~site_code,mod_dataset,median),by.x="site_code",by.y="site_code")

mod_dataset<-add.sajrajty(mod_dataset, extracted.sox, "sox")
mod_dataset<-add.sajrajty(mod_dataset, extracted.sothf, "sothf")
mod_dataset<-add.sajrajty(mod_dataset, extracted.nox, "nox")
mod_dataset<-add.sajrajty(mod_dataset, extracted.nh, "nh")
mod_dataset<-add.sajrajty(mod_dataset, extracted.nsum, "nsum")

mod_dataset<-merge(mod_dataset,soil_dataset_FH[,c("site_code", "pH", "C.N")],by.x="site_code",by.y="site_code")
mod_dataset<-merge(mod_dataset,soil_dataset_L1[,c("site_code", "pH", "C.N")],by.x="site_code",by.y="site_code")

## Renaming columns ####
names(mod_dataset)<-c("site_code", 
                      "species", 
                      "elevation",
                      "coord_x",
                      "coord_y",
                      "management_type", 
                      "site_category",       
                      "sample_depth", 
                      "year", 
                      "max_cambial_age", 
                      "cambial_age_range", 
                      "cambial_age",
                      "cv_TRW", 
                      "sd_RWI",
                      "cv_RWI",
                      "mid_TRW", 
                      "mean_temp",
                      "resid_temp",
                      "mean_cwb",
                      "last5",
                      "normal_temp",
                      "normal_cwb",
                      "max_age",
                      "max_range",
                      "median_age",
                      "median_range",
                      "sox",
                      "sothf",
                      "nox",
                      "nh",
                      "nsum",
                      "pH_FH",  
                      "C.N_FH",
                      "pH_L1",  
                      "C.N_L1")

mod_dataset$resid_cwb<-mod_dataset$normal_cwb-mod_dataset$mean_cwb

mod_dataset<-na.omit(mod_dataset)

## Data transformation - to get normal distribution ####
mod_dataset$log_cv_TRW<-base::log(mod_dataset$cv_TRW)
mod_dataset$log_sd_RWI<-base::log(mod_dataset$sd_RWI)
mod_dataset$log_cv_RWI<-base::log(mod_dataset$cv_RWI)

## Splitting datasets ####
main_mod_dataset<-mod_dataset
abal_mod_dataset<-subset(mod_dataset,species=="ABAL")
fasy_mod_dataset<-subset(mod_dataset,species=="FASY")
pcab_mod_dataset<-subset(mod_dataset,species=="PCAB")
pisy_mod_dataset<-subset(mod_dataset,species=="PISY")
qusp_mod_dataset<-subset(mod_dataset,species=="QUSP")

main_mod_dataset_scaled<-scale.variables(main_mod_dataset)
abal_mod_dataset_scaled<-scale.variables(abal_mod_dataset)
fasy_mod_dataset_scaled<-scale.variables(fasy_mod_dataset)
pcab_mod_dataset_scaled<-scale.variables(pcab_mod_dataset)
pisy_mod_dataset_scaled<-scale.variables(pisy_mod_dataset)
qusp_mod_dataset_scaled<-scale.variables(qusp_mod_dataset)

## ---------------------------------------------------------------------------- ####
## Data saving ####
write.table(main_mod_dataset, "Calculated_datasets/Finalized_datasets/mod_dataset_main.txt", sep=";", dec=".", row.names=F)
write.table(abal_mod_dataset, "Calculated_datasets/Finalized_datasets/mod_dataset_ABAL.txt", sep=";", dec=".", row.names=F)
write.table(fasy_mod_dataset, "Calculated_datasets/Finalized_datasets/mod_dataset_FASY.txt", sep=";", dec=".", row.names=F)
write.table(pcab_mod_dataset, "Calculated_datasets/Finalized_datasets/mod_dataset_PCAB.txt", sep=";", dec=".", row.names=F)
write.table(pisy_mod_dataset, "Calculated_datasets/Finalized_datasets/mod_dataset_PISY.txt", sep=";", dec=".", row.names=F)
write.table(qusp_mod_dataset, "Calculated_datasets/Finalized_datasets/mod_dataset_QUSP.txt", sep=";", dec=".", row.names=F)

write.table(main_mod_dataset_scaled, "Calculated_datasets/Finalized_datasets/scaled_mod_dataset_main.txt", sep=";", dec=".", row.names=F)
write.table(abal_mod_dataset_scaled, "Calculated_datasets/Finalized_datasets/scaled_mod_dataset_ABAL.txt", sep=";", dec=".", row.names=F)
write.table(fasy_mod_dataset_scaled, "Calculated_datasets/Finalized_datasets/scaled_mod_dataset_FASY.txt", sep=";", dec=".", row.names=F)
write.table(pcab_mod_dataset_scaled, "Calculated_datasets/Finalized_datasets/scaled_mod_dataset_PCAB.txt", sep=";", dec=".", row.names=F)
write.table(pisy_mod_dataset_scaled, "Calculated_datasets/Finalized_datasets/scaled_mod_dataset_PISY.txt", sep=";", dec=".", row.names=F)
write.table(qusp_mod_dataset_scaled, "Calculated_datasets/Finalized_datasets/scaled_mod_dataset_QUSP.txt", sep=";", dec=".", row.names=F)

