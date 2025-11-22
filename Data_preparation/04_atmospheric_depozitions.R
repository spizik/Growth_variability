library("terra")
library("openxlsx")

## Data loading ----------------------------------------------------------------
site_list<-site.list
site_list<-site_list[,c("site_code", "species", "elevation", "lon", "lat")]
names(site_list)<-c("site_code", "species", "elevation", "x", "y")

## Extracting SOX --------------------------------------------------------------
sox<-rast("Input_data/Depozitions/deposits_SO.tif")
extracted.sox<-site_list

for(i in names(sox)){
  # i=names(sox)[1]
  extracted.sox[,i] <- terra::extract(sox[[i]], site_list[,c("x","y")])[,2]
}

## Extracting SOthf ------------------------------------------------------------
sothf<-rast("Input_data/Depozitions/deposits_SOTHF.tif")
extracted.sothf<-site_list

for(i in names(sothf)){
  # i=names(sox)[1]
  extracted.sothf[,i] <- terra::extract(sothf[[i]], site_list[,c("x","y")])[,2]
}

## Extracting NOx --------------------------------------------------------------
nox<-rast("Input_data/Depozitions/deposits_NO.tif")
extracted.nox<-site_list

for(i in names(nox)){
  # i=names(sox)[1]
  extracted.nox[,i] <- terra::extract(nox[[i]], site_list[,c("x","y")])[,2]
}

## Extracting Nh ---------------------------------------------------------------
nh<-rast("Input_data/Depozitions/deposits_NH.tif")
extracted.nh<-site_list

for(i in names(nh)){
  # i=names(sox)[1]
  extracted.nh[,i] <- terra::extract(nh[[i]], site_list[,c("x","y")])[,2]
}

## Extracting Nsum -------------------------------------------------------------
nsum<-rast("Input_data/Depozitions/Deposits_NSUM.tif")
extracted.nsum<-site_list

for(i in names(nsum)){
  # i=names(sox)[1]
  extracted.nsum[,i] <- terra::extract(nsum[[i]], site_list[,c("x","y")])[,2]
}

## Exporting data --------------------------------------------------------------
write.xlsx(extracted.sox, "Calculated_datasets/Depozitions/extracted_SO.xlsx")
write.xlsx(extracted.sothf, "Calculated_datasets/Depozitions/extracted_SOthf.xlsx")
write.xlsx(extracted.nox, "Calculated_datasets/Depozitions/extracted_NO.xlsx")
write.xlsx(extracted.nh, "Calculated_datasets/Depozitions/extracted_NH.xlsx")
write.xlsx(extracted.nsum, "Calculated_datasets/Depozitions/extracted_Nsum.xlsx")



