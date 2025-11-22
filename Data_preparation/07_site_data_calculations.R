## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ------------------------------------------------------- event calculations ####
## raw datasets ####
raw_all<-calculate.all.sites.data(prepared.data.tree.means,site.list)

res_all<-prepare.data(raw_all,
                      min_age=20,
                      min_trees=5,
                      min_final_year=2015,
                      min_starting_year=2015)

## -------------------------------------------------------------- data saving ####
write.table(res_all,"Calculated_datasets/Finalized_datasets/dataset.txt",dec=".",sep=";")





























