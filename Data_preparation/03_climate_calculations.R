## Data loading ####
climate.data<-load.files(folder="Calculated_datasets/Site_climate",types="csv",separator=",",decimal=".")
site_list<-site.list

first_month<-4
last_month<-9

## Calculations ####
bals<-do.call(rbind, calculate.bals(climate.data))
speis<-do.call(rbind, calculate.speis(climate.data, site_list))
temps<-do.call(rbind, calculate.temps(climate.data))

write.table(bals, "Calculated_datasets/Recalculated_climate/bals.txt", dec=".", sep=";")
write.table(speis, "Calculated_datasets/Recalculated_climate/speis.txt", dec=".", sep=";")
write.table(temps, "Calculated_datasets/Recalculated_climate/temps.txt", dec=".", sep=";")


