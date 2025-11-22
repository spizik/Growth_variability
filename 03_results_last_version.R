## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## Functions ####
source("01_functions.R")

## Data preparation ####
# Dataset recalculation...
# detrendig_method <- "GAM"
# detrendig_method <- "qGAM"
detrendig_method <- "mean"
# detrendig_method <- "Spline"
# source("02_00_dataset_preparation.R")

## Data loading ####
source("02_01_datasets_loading.R")

## Setting up colours ####
cols.management <- c(
  "Managed"   = "#0072B2", 
  "ManPa"     = "#E69F00",  
  "Unmanaged" = "#CC79A7"   
)

cols.category <- c(
  "moist"    = "#4477AA",  
  "rich"     = "#228833", 
  "moderate" = "#EE6677",
  "acid"     = "#AA3377", 
  "extreme"  = "#66CCEE"   
)


cols.species<-brewer.pal(n = 8, name = 'Dark2')
cols.species<-c(cols.species,"#777777")
names(cols.species)<-c("ABAL","PCAB","FASY","QUSP","FREX","PISY","LADE","PSME","none")

cols.species <- c(
  "ABAL" = "#117733",  # tmavě zelená – jedle
  "PCAB" = "#CC6677",  # lososová – smrk
  "PISY" = "#DDCC77",  # hořčicově žlutá – borovice
  "FASY" = "#332288",  # temná modrofialová – buk
  "QUSP" = "#88CCEE",  # světle modrá – dub
  "none" = "#999999"   # neutrální šedá
)

cols.elevation<-brewer.pal(n = 7, name = 'YlOrRd')
names(cols.elevation)<-c("300","500","700","900","1100","1300","1600")

elev_colors <- c(
  "<400"     = "#3B4CC0",  # tmavě modrá
  "400-600"  = "#2C7BB6",
  "600-800"  = "#00A6CA",
  "800-1000" = "#00CCBC",
  "1000-1200"= "#90D743",
  "2000"     = "#FDE725"
)

cols.class<-c("#26466D","#BBBBBB","#D73027")
names(cols.class)<-c("young","none","old")

val_cols<-c("#26466D","#BBBBBB","#D73027")
names(val_cols)<-c("negative", "none","positive")

cols.methods <- c(
  "mean"   = "#009E73",  # strong green
  "spline" = "#D55E00",  # orange-red
  "GAM"    = "#882255",  # dark purple
  "qGAM"   = "#6699CC"   # steel blue
)

## Elevation zones ####
min_sites <- 5
# elevation_breaks <- c(0, 300, 600, 900, 1200, 2000)
# elevation_labels <- c("<300", "300-600", "600-900", "900-1200", ">1200")
elevation_breaks <- c(0, 400, 600, 800, 1000, 1200, 2000)
elevation_labels <- c("<400", "400-600", "600-800", "800-1000", "1000-1200", "2000")

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## Main results
## -----------------------------------------------------------------------------
## --------------
## Methods ####
print("------------------------------------------------------------------------")
print(paste0("Elevation range: ", min(df.all$elevation), "- ", max(df.all$elevation)))

print("------------------------------------------------------------------------")
print(paste0("Temperature range: ", round(min(clim$MAT),1), "- ", round(max(clim$MAT),1)))

print("------------------------------------------------------------------------")
print(paste0("Precipitation range: ", round(min(clim$MAP),0), "- ", round(max(clim$MAP),0)))


print("------------------------------------------------------------------------")
print(paste0("Pocet situ: ", length(unique(df.all$site_code))))
print(paste0("Pocet stromu: ", sum(aggregate(sample_depth~site_code, df.all, max)$sample_depth)))

print("------------------------------------------------------------------------")
print(paste0("Pocet situ ABAL: ", length(unique(subset(df.all, species=="ABAL")$site_code))))
print(paste0("Pocet stromu ABAL: ", sum(aggregate(sample_depth~site_code, subset(df.all, species=="ABAL"), max)$sample_depth)))

print("------------------------------------------------------------------------")
print(paste0("Pocet situ PCAB: ", length(unique(subset(df.all, species=="PCAB")$site_code))))
print(paste0("Pocet stromu PCAB: ", sum(aggregate(sample_depth~site_code, subset(df.all, species=="PCAB"), max)$sample_depth)))

print("------------------------------------------------------------------------")
print(paste0("Pocet situ PISY: ", length(unique(subset(df.all, species=="PISY")$site_code))))
print(paste0("Pocet stromu PISY: ", sum(aggregate(sample_depth~site_code, subset(df.all, species=="PISY"), max)$sample_depth)))

print("------------------------------------------------------------------------")
print(paste0("Pocet situ FASY: ", length(unique(subset(df.all, species=="FASY")$site_code))))
print(paste0("Pocet stromu FASY: ", sum(aggregate(sample_depth~site_code, subset(df.all, species=="FASY"), max)$sample_depth)))

print("------------------------------------------------------------------------")
print(paste0("Pocet situ QUSP: ", length(unique(subset(df.all, species=="QUSP")$site_code))))
print(paste0("Pocet stromu QUSP: ", sum(aggregate(sample_depth~site_code, subset(df.all, species=="QUSP"), max)$sample_depth)))

## --------------
## --------------
## Figure 1 ####

# teploty a srazky jed situ
# srovnani prumerneho SPEI a SPEI vybranzch suchych roku
# ddto pro teploty

# figure making
source("Graphs/Figure_1.R")

# figure
w=22 ; h=32
ggsave("outputs/main_figure_1.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/main_figure_1.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Figure 2 ####
## figure making
source("Graphs/Figure_2.R")

# figure output
w=20 ; h=36
ggsave("outputs/main_figure_2.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/main_figure_2.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Figure 3 ####
## figure making
source("Graphs/Figure_3.R")

# figure output
w=10 ; h=8
ggsave("outputs/main_figure_3.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/main_figure_3.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Figure 4 ####
source("Graphs/Figure_4.R")

# figure output
w=20 ; h=36
ggsave("outputs/main_figure_4.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/main_figure_4.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Figure 5 ####
# NOTE:
# The fitted models used in the manuscript are not included in the repository.
# They cannot be estimated on the reduced demo dataset and therefore are not loaded here.

# source("Graphs/Figure_5.R")

# source("Data_preparation/11_testing_calculated_models.R")

# Model rsqr
# r.squaredGLMM(mod_main)
# r.squaredGLMM(mod_abal)
# r.squaredGLMM(mod_pcab)
# r.squaredGLMM(mod_pisy)
# r.squaredGLMM(mod_fasy)
# r.squaredGLMM(mod_qusp)
# 
# r.squaredGLMM(mod_main)[2] - r.squaredGLMM(mod_main)[1]
# r.squaredGLMM(mod_abal)[2] - r.squaredGLMM(mod_abal)[1]
# r.squaredGLMM(mod_pcab)[2] - r.squaredGLMM(mod_pcab)[1]
# r.squaredGLMM(mod_pisy)[2] - r.squaredGLMM(mod_pisy)[1]
# r.squaredGLMM(mod_fasy)[2] - r.squaredGLMM(mod_fasy)[1]
# r.squaredGLMM(mod_qusp)[2] - r.squaredGLMM(mod_qusp)[1]

# spatial autocorelation in model
# calculate_sp_autocorrel_in_model(mod_main, main_mod_dataset_scaled)
# calculate_sp_autocorrel_in_model(mod_abal, abal_mod_dataset_scaled)
# calculate_sp_autocorrel_in_model(mod_pcab, pcab_mod_dataset_scaled)
# calculate_sp_autocorrel_in_model(mod_pisy, pisy_mod_dataset_scaled)
# calculate_sp_autocorrel_in_model(mod_fasy, fasy_mod_dataset_scaled)
# calculate_sp_autocorrel_in_model(mod_qusp, qusp_mod_dataset_scaled)

# figure output
# w=52 ; h=28
# ggsave("outputs/main_figure_5.png",figure,dpi=600,width=w,height=h,units="cm")
# ggsave("outputs/main_figure_5.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## --------------
## Supplementary figure 1 ####
source("Graphs/Supplementary_Figure_01.R")

w=24 ; h=40
ggsave("outputs/supplementary_figure_01.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_01.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 2 ####
source("Graphs/Supplementary_Figure_02.R")

w=24 ; h=40
ggsave("outputs/supplementary_figure_02.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_02.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 3 ####
source("Graphs/Supplementary_Figure_03.R")

w=24 ; h=40
# w=36 ; h=40
ggsave("outputs/supplementary_figure_03.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_03.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 4 ####
source("Graphs/Supplementary_Figure_04.R") 

# figure output
w=24 ; h=32
ggsave("outputs/supplementary_figure_04.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_04.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 5 ####
source("Graphs/Supplementary_Figure_05.R")

# figure output
# w=24 ; h=8
w=30 ; h=40
ggsave("outputs/supplementary_figure_05.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_05.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 6 ####
source("Graphs/Supplementary_Figure_06.R")

# figure output
w=24 ; h=8
# w=14 ; h=10
ggsave("outputs/supplementary_figure_06.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_06.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 7 ####
source("Graphs/Supplementary_Figure_07.R")

# figure output
# w=24 ; h=8
w=24 ; h=35
ggsave("outputs/supplementary_figure_07.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_07.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 8 ####
source("Graphs/Supplementary_Figure_08.R")

# figure output
# w=24 ; h=8
w=24 ; h=24
ggsave("outputs/supplementary_figure_08.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_08.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 9 ####
source("Graphs/Supplementary_Figure_09.R")

# figure output
# w=24 ; h=8
w=24 ; h=30
ggsave("outputs/supplementary_figure_09.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_09.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 10 ####
source("Graphs/Supplementary_Figure_10.R")

# figure output
w=12 ; h=24
# w=24 ; h=30
ggsave("outputs/supplementary_figure_10.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_10.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 11 ####
source("Graphs/Supplementary_Figure_11.R")

# figure output
# w=24 ; h=8
w=24 ; h=30
ggsave("outputs/supplementary_figure_11.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_11.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 12 ####
source("Graphs/Supplementary_Figure_12.R")

# figure output
# w=24 ; h=8
w=24 ; h=36
ggsave("outputs/supplementary_figure_12.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_12.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)

## --------------
## Supplementary figure 13 ####
source("Graphs/Supplementary_Figure_13.R")

# figure output
# w=24 ; h=8
w=24 ; h=16
ggsave("outputs/supplementary_figure_13.png",figure,dpi=600,width=w,height=h,units="cm")
ggsave("outputs/supplementary_figure_13.pdf",figure,dpi=600,width=w/2.54,height=h/2.54)



