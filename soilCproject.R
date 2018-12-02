#TO-DO correct for bag weights [DONE]
#double-check entered data [DONE]
#working up of bulk density data [DONE]
#re-vist how low volume of coarse fragments are handled
mainDir <- 'C:/Users/smdevine/Desktop/rangeland project'
results <- 'C:/Users/smdevine/Desktop/rangeland project/results'
vol_0_10cm <- pi * (4.7 / 2)^2 * 10 #id of soil core liner was 47 mm
vol_10_30cm <- 2 * vol_0_10cm
vol_0_5cm <- 0.5 * vol_0_10cm
BD_data <- read.csv(file.path('C:/Users/smdevine/Desktop/rangeland project/soils_data', 'soils_datasheet.7.23.18.csv'), stringsAsFactors = FALSE)
BD_data$sample_code
summary(BD_data$frags_g[BD_data$frags_g > 0 & !is.na(BD_data$frags_mL)] / BD_data$frags_mL[BD_data$frags_g > 0 & !is.na(BD_data$frags_mL)])
avgBD_frags <- mean(BD_data$frags_g[BD_data$frags_g > 0] / BD_data$frags_mL[BD_data$frags_g > 0], na.rm = TRUE)
BD_data$frags_mL[is.na(BD_data$frags_mL) & BD_data$frags_g > 0] <- BD_data$frags_g[is.na(BD_data$frags_mL) & BD_data$frags_g > 0] / avgBD_frags
BD_data$bag_weight <- ifelse(BD_data$bag_type==1, 4.8, ifelse(BD_data$bag_type==2, 3.9, ifelse(BD_data$bag_type==3, 6.2, 6.4))) #bag type 4 was 6.4 g
summary(BD_data$bag_weight)
BD_data$fines_g_AD <- BD_data$dry_wt_g - BD_data$frags_g - BD_data$bag_weight #subtract out frag and bag mass
##need to correct samples 2-1 and 3-1 because they were only 0-5 cm samples
BD_data$depth_code[BD_data$sample_code=='2_1' | BD_data$sample_code=='3_1'] <- 3
BD_data$core_volume_mL <- ifelse(BD_data$depth_code==1, vol_0_10cm, ifelse(BD_data$depth_code==2, vol_10_30cm, vol_0_5cm))
BD_data$fines_mL <- BD_data$core_volume_mL - BD_data$frags_mL
BD_data$frags_vol_perc <- 100 * (BD_data$frags_mL / BD_data$core_volume_mL)
BD_data$g_H2O_AD <- BD_data$moist_wt_g - BD_data$dry_wt_g #moist_wt_g and dry_wt_g included bag weights
#read in oven-dry moisture measurements before calculating volumetric moisture
OD_data <- read.csv(file.path('C:/Users/smdevine/Desktop/rangeland project/soils_data', 'airdry_soilmoisture.csv'), stringsAsFactors = FALSE)
colnames(OD_data)[1] <- 'sample_code'
OD_data$sample_code <- gsub('-1', '_1', OD_data$sample_code)
OD_data$sample_code <- gsub('-2', '_2', OD_data$sample_code)
BD_data <- merge(BD_data, OD_data, by = 'sample_code')
BD_data$fines_g_OD <- BD_data$fines_g_AD - (BD_data$percent.H2O.AD.basis / 100)  * BD_data$fines_g_AD
BD_data$g_H2O_total <- BD_data$g_H2O_AD + (BD_data$percent.H2O.AD.basis / 100)  * BD_data$fines_g_AD
BD_data$H2O_gravimetric <- BD_data$g_H2O_total / BD_data$fines_g_OD #percent.H2O.AD.basis the percentage of air-dry moisture on an air-dry basis after drying at 100 C for 24 hrs, assuming coarse fragments did not contribute moisture weight as part of pre-sieved mass
BD_data$H2O_volumetric <- BD_data$g_H2O_total / BD_data$core_volume_mL
hist(BD_data$H2O_gravimetric)
hist(BD_data$H2O_volumetric)
BD_data$bulk_density_g_cm3 <- BD_data$fines_g_OD / BD_data$fines_mL
# BD_data$whole_soil_g_AD <- BD_data$dry_wt_g - BD_data$bag_weight
# BD_data$whole_soil_g_OD <- BD_data$whole_soil_g_AD - (BD_data$percent.H2O.AD.basis / 100)  * BD_data$fines_g_AD #subtract out water associated with fine fraction and assume coarse fraction had no water after air-drying inside bags before sieving
# BD_data$bulk_density_g_cm3_whole_soil <- BD_data$whole_soil_g_OD / BD_data$core_volume_mL
dim(BD_data)
summary(BD_data$bulk_density_g_cm3)
hist(BD_data$bulk_density_g_cm3)
hist(BD_data$bulk_density_g_cm3_whole_soil)
summary(BD_data$bulk_density_g_cm3_whole_soil - BD_data$bulk_density_g_cm3)
hist(BD_data$frags_vol_perc)
summary(BD_data$frags_vol_perc)
write.csv(BD_data, file.path(results, 'soil_data', paste0('soilBD_H2O_frags_', Sys.Date(), '.csv')), row.names=FALSE)
#read-in BD_data
BD_data <- read.csv(file.path(results, 'soil_data', 'soilBD_H2O_frags_2018-09-17.csv'), stringsAsFactors = FALSE)
#split dataset by depth
BD_data_0_10cm <- BD_data[BD_data$depth_code==1 | BD_data$depth_code==3, ]
BD_data_10_30cm <- BD_data[BD_data$depth_code==2,]

#read in shapefile of points and merge with two dataframes above
library(raster)
sampling_pts <- shapefile(file.path(mainDir, 'sampling points 2018', 'soil_sampling_points.shp'))
biomass_Apr2017 <- raster(file.path(mainDir, 'sampling strategy April 2018', 'Biomass_2017-04-10_APAR.tif'))
plot(biomass_Apr2017)
sampling_pts$point_no <- as.integer(gsub('point', '', sampling_pts$Comment))
sampling_pts$point_no
plot(sampling_pts, pch=16, col = 'brown')
text(sampling_pts, sampling_pts$point_no)
sampling_pts$Apr2017biomass <- extract(biomass_Apr2017, coordinates(sampling_pts)[,1:2], buffer=1.5, fun=mean) #this calculates a mean of all 30 x 30 cm cells within 1.5 m of soil sampling point
point_data_0_10cm <- merge(sampling_pts, BD_data_0_10cm, 'point_no')
point_data_0_10cm <- point_data_0_10cm[order(point_data_0_10cm$point_no), ]
plot(point_data_0_10cm, cex=point_data_0_10cm$Apr2017biomass/1000, pch=1)
plot(point_data_0_10cm, cex=1.5 * point_data_0_10cm$bulk_density_g_cm3, pch=1, main='0-10')
plot(point_data_0_10cm, cex=point_data_0_10cm$frags_vol_perc / 3, pch=1, main='0-10 coarse fragment % by vol')
plot(point_data_0_10cm, cex=point_data_0_10cm$H2O_volumetric * 20, pch=1)
plot(point_data_0_10cm$H2O_gravimetric, point_data_0_10cm$Apr2017biomass)
plot(point_data_0_10cm$bulk_density_g_cm3, point_data_0_10cm$Apr2017biomass)
plot(point_data_0_10cm$frags_vol_perc, point_data_0_10cm$Apr2017biomass)
plot(point_data_0_10cm$H2O_volumetric, point_data_0_10cm$Apr2017biomass) #all this should be done against 2018 biomass
summary(lm(Apr2017biomass ~ bulk_density_g_cm3, data = point_data_0_10cm))
summary(lm(Apr2017biomass ~ bulk_density_g_cm3_whole_soil, data = point_data_0_10cm))
summary(lm(Apr2017biomass ~ H2O_volumetric, data = point_data_0_10cm)) #should be done against 2018 biomass
summary(lm(Apr2017biomass ~ frags_vol_perc, data = point_data_0_10cm))
summary(lm(Apr2017biomass ~ H2O_gravimetric, data = point_data_0_10cm))
plot(point_data_0_10cm$frags_vol_perc, point_data_0_10cm$H2O_volumetric)
plot(point_data_0_10cm$frags_vol_perc, point_data_0_10cm$OD.H2O.perc)

#now the 10-30 cm data
point_data_10_30cm <- merge(sampling_pts, BD_data_10_30cm, 'point_no')
#plot(point_data_10_30cm, cex=point_data_10_30cm$Apr2017biomass/1000, pch=1)
point_data_10_30cm <- point_data_10_30cm[order(point_data_10_30cm$point_no), ]
plot(point_data_10_30cm, cex=1.5 * point_data_10_30cm$bulk_density_g_cm3, pch=1)
plot(point_data_10_30cm, cex=point_data_10_30cm$frags_vol_perc / 3, pch=1, main='10-30 cm coarse fragment % by vol')
plot(point_data_10_30cm, cex=point_data_10_30cm$H2O_volumetric * 20, pch=1)
plot(point_data_10_30cm$bulk_density_g_cm3, point_data_10_30cm$Apr2017biomass)
plot(point_data_10_30cm$frags_vol_perc, point_data_10_30cm$Apr2017biomass)
plot(point_data_10_30cm$H2O_volumetric, point_data_10_30cm$Apr2017biomass) #should be done against 2018 biomass
summary(lm(Apr2017biomass ~ bulk_density_g_cm3, data = point_data_10_30cm))
summary(lm(Apr2017biomass ~ bulk_density_g_cm3_whole_soil, data = point_data_10_30cm))
summary(lm(Apr2017biomass ~ frags_vol_perc, data = point_data_10_30cm))
summary(lm(Apr2017biomass ~ H2O_volumetric, data = point_data_10_30cm)) #should be done against 2018 biomass

#multiple linear regression
summary(lm(Apr2017biomass ~ bulk_density_g_cm3 + frags_vol_perc, data = point_data_10_30cm))
summary(lm(Apr2017biomass ~ bulk_density_g_cm3 + frags_vol_perc + H2O_volumetric, data = point_data_0_10cm))
summary(lm(Apr2017biomass ~ OD.H2O.perc + frags_vol_perc, data = point_data_0_10cm))

#write to shapefiles
fname <- file.path('C:/Users/smdevine/Desktop/RangelandProject/sampling points 2018', 'soil_data_0_10cm.shp')
fname <- gsub("\\\\", "/", fname)
shapefile(point_data_0_10cm, file.path(results, 'shapefiles', 'point_data_0_10cm.shp'), overwrite=TRUE)
shapefile(point_data_10_30cm, file.path(results, 'shapefiles', 'point_data_10_30cm.shp'), overwrite=TRUE)

#read in sensor locations
sensor.locations <- shapefile(file.path(mainDir, 'soil.moisture.sensors', 'sensor_forage2017.shp'))
plot(sensor.locations, add=TRUE, col='red')
text(sensor.locations, labels=sensor.locations$location, offset=0.2, pos=1)
plot(sensor.locations, cex=sensor.locations$clp041017/1000, pch=1)

#see data at point 88 and 56
as.data.frame(point_data_0_10cm[point_data_0_10cm$sample_code=='88_1',])
100*(6.3+247.1*0.042)/vol_0_10cm #9.29% vol H2O
as.data.frame(point_data_0_10cm[point_data_0_10cm$sample_code=='56_1',])
100*(5.5+256.8*0.032)/179.45 #7.64% vol H2O
