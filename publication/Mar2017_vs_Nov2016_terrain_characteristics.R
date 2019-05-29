terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
library(raster)
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

list.files(file.path(terrainDir, 'filtered_Hogan'))
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan', '3m_filtered'), full.names = TRUE))
Mar2017_terrain_3m
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
Mar2017_terrain_3m$annual_kwh.m2 <- solrad_raster
pts_terrain_summary_3m_Mar2017 <- extract(Mar2017_terrain_3m, soil_0_30cm_shp, df=TRUE)
colnames(pts_terrain_summary_3m_Mar2017)[1] <- 'location'

#Nov 2016 DSM
list.files(file.path(terrainDir, 'filtered_Grace'))
Nov2016_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Grace'), full.names = TRUE))
Nov2016_terrain_3m
names(Nov2016_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster_Nov2016 <- raster(file.path(solradDir, 'solrad_3m_filtered_Nov2016.tif'))
solrad_raster_Nov2016 <- solrad_raster_Nov2016 / 1000
Nov2016_terrain_3m$annual_kwh.m2 <- solrad_raster_Nov2016
pts_terrain_summary_3m_Nov2016 <- extract(Nov2016_terrain_3m, soil_0_30cm_shp, df=TRUE)
pts_terrain_summary_3m_Nov2016
colnames(pts_terrain_summary_3m_Nov2016)[1] <- 'location'
#write.csv(pts_terrain_summary_3m_Nov2016, file.path(results, 'terrain_characteristics', 'terrain_3m_filtered_Nov2016.csv'), row.names = FALSE)

#make comparisions
summary(lm(pts_terrain_summary_3m_Mar2017$elevation ~ pts_terrain_summary_3m_Nov2016$elevation)) #0.9942
summary(lm(pts_terrain_summary_3m_Mar2017$aspect ~ pts_terrain_summary_3m_Nov2016$aspect)) #0.6657
summary(lm(pts_terrain_summary_3m_Mar2017$curvature_mean ~ pts_terrain_summary_3m_Nov2016$curvature_mean)) #0.9348
summary(lm(pts_terrain_summary_3m_Mar2017$curvature_plan ~ pts_terrain_summary_3m_Nov2016$curvature_plan)) #0.9108
summary(lm(pts_terrain_summary_3m_Mar2017$curvature_profile ~ pts_terrain_summary_3m_Nov2016$curvature_profile)) #0.9432
summary(lm(pts_terrain_summary_3m_Mar2017$slope ~ pts_terrain_summary_3m_Nov2016$slope)) #0.9864
summary(lm(pts_terrain_summary_3m_Mar2017$annual_kwh.m2 ~ pts_terrain_summary_3m_Nov2016$annual_kwh.m2)) #0.9974
