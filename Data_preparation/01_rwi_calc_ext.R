


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
                      # TRW=tree.core[,5]/100,
                      TRW=as.numeric(as.character(tree.core[,5]))/100)

## Meant TRW data
prepared.data.tree.means<-prepare.dataset(site.list, tree.core)


## calculates RWI
for(i in names(prepared.data.tree.means)){
  
  # print(i)
  write.table(create.rwi(prepared.data.tree.means[[i]], "Spline"), paste0("Calculated_datasets/rwi_data_Spline/",i,".txt"), dec = ".", sep = ";")
  write.table(create.rwi(prepared.data.tree.means[[i]], "GAM"), paste0("Calculated_datasets/rwi_data_GAM/",i,".txt"), dec = ".", sep = ";")
}