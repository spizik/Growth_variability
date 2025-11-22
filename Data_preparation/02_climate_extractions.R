## allows uploading a text file with coordinates of any number of points and extracting time series for them
## currently set to use grids covering the period 1961–2020; a minor update will be needed when new grids become available

## uploading coordinates: the file should contain the columns CODE, Latitude, and Longitude
site_list <- site.list
site_list <- site_list[, c("site_code", "lat", "lon")]
names(site_list) <- c("CODE", "Latitude", "Longitude")

temp <- nc_open("Input_data/Clima_grids/t_wgs_month.nc")
prec <- nc_open("Input_data/Clima_grids/sra_wgs_month.nc")

STORAGE_T <- data.frame(TIMESTAMP = seq(from = as.Date("19610101", tryFormats = c("%Y%m%d")), to = as.Date("20201201", tryFormats = c("%Y%m%d")), by = "month"))
STORAGE_P <- data.frame(TIMESTAMP = seq(from = as.Date("19610101", tryFormats = c("%Y%m%d")), to = as.Date("20201201", tryFormats = c("%Y%m%d")), by = "month"))

for (i in c(1:nrow(site_list))) {
  
  # Loading coordinates
  site1x <- site_list[i, "Longitude"]; site1y <- site_list[i, "Latitude"]
  
  # selecting specific pixel
  Order.X <- data.frame(ORDER = c(1:temp$dim$x$len), GRID = ncvar_get(nc = temp, varid = "x"))
  Order.Y <- data.frame(ORDER = c(1:temp$dim$y$len), GRID = ncvar_get(nc = temp, varid = "y"))
  Order.X$DIFFERENCE <- abs(Order.X$GRID - site1x)
  Order.Y$DIFFERENCE <- abs(Order.Y$GRID - site1y)
  site1x.order <- Order.X[Order.X$DIFFERENCE == min(Order.X$DIFFERENCE),"ORDER"] # selecting closest pixel
  site1y.order <- Order.Y[Order.Y$DIFFERENCE == min(Order.Y$DIFFERENCE),"ORDER"]
  
  # Data extraction
  STORAGE_T[,site_list[i,"CODE"]] <- ncvar_get(nc = temp, varid = "T", start = c(site1x.order, site1y.order, 1), count = c(1,1,-1)) 
  STORAGE_P[,site_list[i,"CODE"]] <- ncvar_get(nc = prec, varid = "SRA", start = c(site1x.order, site1y.order, 1), count = c(1,1,-1)) 
  
}


for(i in site_list$CODE){
  
  # print(i)
  
  sub_T <- STORAGE_T[,c("TIMESTAMP", i)]
  sub_P <- STORAGE_P[,c("TIMESTAMP", i)]
  
  output <- data.frame(Date = sub_T$TIMESTAMP,
                       Year = as.numeric(as.character(substr(sub_T$TIMESTAMP, 1, 4))),
                       Month = as.numeric(as.character(substr(sub_T$TIMESTAMP, 6, 7))),
                       Temp = sub_T[,i],
                       Prec = sub_P[,i])
  
  write.table(output, paste0("Calculated_datasets/Site_climate/",i,"_clim.csv"), sep=",", dec=".")
  
}