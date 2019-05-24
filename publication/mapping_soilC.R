terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
FiguresDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/Figures'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
library(raster)

#read-in intermediate results for 0-30 cm, 0-10 cm, and 10-30 cm
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_0_10cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_10_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

#read-in terrain properties
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan', '3m_filtered'), full.names = TRUE))
names(Mar2017_terrain_3m)
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
Mar2017_terrain_3m$annual_kwh.m2 <- solrad_raster

NDVI_stack <- stack(list.files(NDVIDir, full.names = TRUE))
names(NDVI_stack)
# cellStats(NDVI_stack, 'mean')
NDVI_2017 <- NDVI_stack[[c(1,3,5,7,9)]]
NDVI_2018 <- NDVI_stack[[c(2,4,6,8)]]
meanNDVI_2017 <- calc(NDVI_2017, fun=mean)#, filename=file.path(FiguresDir, 'meanNDVI2017.tif'))
meanNDVI_2018 <- calc(NDVI_2018, fun=mean)#, filename=file.path(FiguresDir, 'meanNDVI2018.tif'))

#create grid for spatial predictions and upscaling NDVI & reflectance data to same grid
r <- raster(Mar2017_terrain_3m)
e <- extent(meanNDVI_2017)
r <- crop(r, e)
Mar2017_terrain_3m_cropped <- crop(Mar2017_terrain_3m, r)
Mar2017_terrain_3m_cropped$NDVI_2017mean_1m <- resample(meanNDVI_2017, Mar2017_terrain_3m_cropped)

#map organic carbon 0-30 cm
#5 var model
#Fig 4a
lm_terrain5_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
summary(lm_terrain5_0_30cm)
kgOrgC.m2_terrain5_0_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_0_30cm)
#writeRaster(kgOrgC.m2_terrain5_0_30cm, filename=file.path(FiguresDir, 'Fig4a.tif'))
round(quantile(kgOrgC.m2_terrain5_0_30cm), 2)
#  0%  25%  50%  75% 100% 
#0.55 3.16 3.51 3.89 5.57  
# plot(kgOrgC.m2_terrain5_0_30cm)

#2 var model
#Fig 4b
lm_terrain2_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
kgOrgC.m2_terrain2_0_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain2_0_30cm) 
writeRaster(kgOrgC.m2_terrain2_0_30cm, filename=file.path(FiguresDir, 'Fig4b.tif'))
round(quantile(kgOrgC.m2_terrain2_0_30cm), 2)
#   0%   25%   50%   75%  100% 
#-0.44  3.24  3.60  3.92  5.01
# plot(kgOrgC.m2_terrain2_0_30cm)
sum(kgOrgC.m2_terrain2_0_30cm < 0) #not sure why this doesn't work anymore
predictions_2var <- getValues(kgOrgC.m2_terrain2_0_30cm)
sum(predictions_2var < 0, na.rm = TRUE) #only one negative value

#do the same mapping for 0-10 cm layer
lm_terrain5_0_10cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_10cm_shp)
summary(lm_terrain5_0_10cm)
kgOrgC.m2_terrain5_0_10cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_0_10cm, filename=file.path(FiguresDir, 'Appendices', 'kgOrgC_m2_MLRbest5var_0_10cm_FINAL.tif'))
# plot(kgOrgC.m2_terrain5_0_10cm)


#do the same mapping for 10-30 cm layer
#need to get NIR into covariate stack
lm_terrain5_10_30cm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + annual_kwh.m2 + NIR_meanGS2017 + NDVI_2017mean_1m, data =  soil_10_30cm_shp)
kgOrgC.m2_terrain5_10_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_10_30cm, filename=file.path(FiguresDir, 'Appendices', 'kgOrgC_m2_MLRbest5var_10_30cm_FINAL.tif'))
plot(kgOrgC.m2_terrain5_10_30cm)