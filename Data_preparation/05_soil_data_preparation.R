library("openxlsx")
library("sf")

## Data loading and ommiting ####
site_list <-site.list
soil_data<-read.xlsx("Input_data/Soil_data/Pudy_rajonizace.xlsx")

site_list <- site_list[, c("site_code", "species", "site_category", "lat", "lon")]


## Transform WGS to s-jtsk ####
# Convert to sf object
sf_points <- st_as_sf(site_list, coords = c("lon", "lat"), crs = 4326) # WGS84

# Transform to S-JTSK (EPSG:5514)
sf_points_sjtsk <- st_transform(sf_points, 5514)

# Convert back to data.frame
df_sjtsk <- as.data.frame(st_coordinates(sf_points_sjtsk))

# Combine with original data
site_list <- cbind(site_list, df_sjtsk)


## assining soil data to points ####
output_FH<-NULL
output_FH_layer1_mean<-NULL
output_layer1<-NULL
output_layer2<-NULL

cols_to_pick <- c("pHH2O", "C", "N", "Ptot", "Ca", "K", "Mg", "C.N", "RC1", "RC2", "RC3")
cols_to_rename <- c("pH", "C", "N", "P", "Ca", "K", "Mg", "C.N", "RC1", "RC2", "RC3")

for(i in 1:nrow(site_list)){
  
  # print(i)
  
  sub.site <- site_list[i,]
  
  soil_data_temp <- soil_data
  soil_data_temp$X_diff <- soil_data_temp$JTSK.x - sub.site$X
  soil_data_temp$Y_diff <- soil_data_temp$JTSK.y - sub.site$Y
  soil_data_temp$distance <- sqrt(soil_data_temp$X_diff^2 + soil_data_temp$Y_diff^2)
  soil_data_temp <- soil_data_temp[order(soil_data_temp$distance),]
  
  soil_data_temp <- soil_data_temp[which(soil_data_temp$distance == min(soil_data_temp$distance)),]
  
  temp_output_FH <- soil_data_temp[which(soil_data_temp$Horizon == "FH"), c(cols_to_pick)]
  temp_output_layer1 <- soil_data_temp[which(soil_data_temp$Horizon == "0-30 cm"), c(cols_to_pick)]
  temp_output_layer2 <- soil_data_temp[which(soil_data_temp$Horizon == "30-80 cm"), c(cols_to_pick)]
  
  temp_output_FH_layer1_mean <- t(as.data.frame(colMeans(soil_data_temp[which(soil_data_temp$Horizon %in% c("FH", "0-30 cm")),c(cols_to_pick)])))
  
  names(temp_output_FH) <- cols_to_rename
  names(temp_output_layer1) <- cols_to_rename
  names(temp_output_layer2) <- cols_to_rename
  names(temp_output_FH_layer1_mean) <- cols_to_rename
  
  output_FH <- rbind(output_FH, cbind(sub.site, temp_output_FH))
  output_layer1 <- rbind(output_layer1, cbind(sub.site, temp_output_layer1))
  output_layer2 <- rbind(output_layer2, cbind(sub.site, temp_output_layer2))
  output_FH_layer1_mean <- rbind(output_FH_layer1_mean, cbind(sub.site, temp_output_FH_layer1_mean))
  
}

## Writing the data ####
write.table(output_FH,"Calculated_datasets/Soil_data/soil_FH.txt",dec=".",sep=";")
write.table(output_layer1,"Calculated_datasets/Soil_data/soil_lay1.txt",dec=".",sep=";")
write.table(output_layer2,"Calculated_datasets/Soil_data/soil_lay2.txt",dec=".",sep=";")
write.table(output_FH_layer1_mean,"Calculated_datasets/Soil_data/soil_FH_and_lay1.txt",dec=".",sep=";")



