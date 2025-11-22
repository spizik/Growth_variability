## -----------------------------------------------------------------------------
## ----------------------------------------------------------data preparation ####
dataset<-data.frame(std=numeric(), res=numeric(), samp.depth=numeric(), site_code=character(), species=character(), year=numeric())
for(i in names(crn.data)){
  temp<-crn.data[[i]]
  temp$site_code<-i
  temp$species<-substr(i,8,11)
  temp$year<-as.numeric(as.character(rownames(temp)))
  temp$site_age<-c(1:nrow(temp))
  dataset<-rbind(dataset,temp)
}

sub <- subset(dataset, year == 1961)
sub <- subset(sub, samp.depth >=5)$site_code
dataset <- dataset[which(dataset$site_code %in% sub),]
dataset <- subset(dataset, site_age >= 20)

dataset <- subset(dataset, year>=1961 & year<=2020)
dataset <- unify.categories(dataset)
## ------------------------------------------------bootstrapping chronologies ####

dta_abal<-bootstrap.chronologies(dataset,"ABAL")
dta_pcab<-bootstrap.chronologies(dataset,"PCAB")
dta_pisy<-bootstrap.chronologies(dataset,"PISY")
dta_qusp<-bootstrap.chronologies(dataset,"QUSP")
dta_fasy<-bootstrap.chronologies(dataset,"FASY")

## ---------------------------------------------------------------data saving ####
write.table(dataset,
            "Calculated_datasets/Bootstrapped_crn_variability/chronologies_individual.txt",
            dec=".",sep=";")

write.table(rbind(dta_abal,
                  dta_pcab,
                  dta_pisy,
                  dta_fasy,
                  dta_qusp),
            "Calculated_datasets/Bootstrapped_crn_variability/chronologies_data.txt",
            dec=".",sep=";")
