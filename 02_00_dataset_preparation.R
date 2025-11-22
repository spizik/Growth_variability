## Funkce
source("01_functions.R")

## Loading site metadata ####
site.list<-read.table("Input_data/Site_meta/sites_all.csv",sep=",",header=T)
site.list<-unify.categories(site.list)

## Creates tree.core.list
tree.core<-read.table("Input_data/Database_file/sampleTable.csv",sep=",", header=T)

tree.core<-data.frame(site=tree.core[,1],
                      species=substr(tree.core[,1],8,11),
                      tree=formatC(tree.core[,2],width=4,flag="0"),
                      core=tree.core[,3],
                      year=tree.core[,4],
                      TRW=as.numeric(as.character(tree.core[,5]))/100)

## Meant TRW data
prepared.data.tree.means<-prepare.dataset(site.list, tree.core)

## Preparing data ####
source("Data_preparation/01_rwi_calc.R")

## Creating individual environmental datasets on site level
source("Data_preparation/02_climate_extractions.R")
source("Data_preparation/03_climate_calculations.R")
source("Data_preparation/04_atmospheric_depozitions.R")
source("Data_preparation/05_soil_data_preparation.R")

## Creating dataset of tree growth and growth variability
rwi.data<-load.files("Calculated_datasets/rwi_data", "txt")
crn.data<-create.chronologies(rwi.data)

source("Data_preparation/06_cv_chronologies.R")
source("Data_preparation/07_site_data_calculations.R")
source("Data_preparation/08_model_dataset_preparation.R")


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

## Model calculation ####
# NOTE:
# The fitted models used in the manuscript are not included in the repository.
# They cannot be estimated on the reduced demo dataset and therefore are not loaded here.

# source("Data_preparation/09_models_calculation.R")
