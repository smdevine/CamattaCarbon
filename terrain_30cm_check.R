mainDir <- 'C:/Users/smdevine/Desktop/rangeland project'
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data/soil C'
soilDataDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
spatialforageDir <- 'C:/Users/smdevine/Desktop/rangeland project/clip_plots/results'
forageDir <- 'C:/Users/smdevine/Desktop/rangeland project/clip_plots'
#0_30 dataset
modelResults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/model_results'
FiguresDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Figures'
CarbonDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject'
ResultsDir <- 'C:/Users/smdevine/Desktop/rangeland project/results'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
ReflectanceDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Reflectance'
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
normalize_var <- function(x) {
  (x - mean(x)) / sd(x)
}
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
library(raster)
library(car)
library(gstat)
library(spdep)
library(automap)
set.seed(20161203)
library(dismo)
kf <- kfold(1:105, k=20)
#read-in intermediate results for 0-30 cm, 0-10 cm, and 10-30 cm
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_0_10cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_10_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

#fix 0-30 cm weighting scheme for texture, BD
sum(soil_0_10cm_shp$point_no - soil_10_30cm_shp$point_no) #check again that location order is the same in each dataset
sum(soil_0_10cm_shp$point_no - soil_0_30cm_shp$point_no)
soil_0_30cm_shp$orgC.percent <- soil_0_10cm_shp$orgC.percent * (soil_0_10cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD)) + soil_10_30cm_shp$orgC.percent * (soil_10_30cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD))
soil_0_30cm_shp$bulk_density_g_cm3 <- soil_0_10cm_shp$bulk_density_g_cm3 * (soil_0_10cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD)) + soil_10_30cm_shp$bulk_density_g_cm3 * (soil_10_30cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD))
soil_0_30cm_shp$PH <- soil_0_10cm_shp$PH * (soil_0_10cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD)) + soil_10_30cm_shp$PH * (soil_10_30cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD))
#fix textural weighting scheme
soil_0_30cm_shp$sand_wtd <- soil_0_10cm_shp$SAND * (soil_0_10cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD)) + soil_10_30cm_shp$SAND * (soil_10_30cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD))
soil_0_30cm_shp$silt_wtd <- soil_0_10cm_shp$SILT * (soil_0_10cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD)) + soil_10_30cm_shp$SILT * (soil_10_30cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD))
soil_0_30cm_shp$clay_wtd <- soil_0_10cm_shp$CLAY * (soil_0_10cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD)) + soil_10_30cm_shp$CLAY * (soil_10_30cm_shp$fines_g_OD / (soil_0_10cm_shp$fines_g_OD + soil_10_30cm_shp$fines_g_OD))

#read-in 30 cm terrain and solrad data
terrain_stack30cm <- stack(list.files(file.path(terrainDir, 'filtered_Hogan', '30cm_filtered'), full.names = TRUE))
names(terrain_stack30cm) <- c('elevation', 'curvature_mean', 'slope', 'annual_kwh.m2')
soil_0_30cm_df$elevation_1m <- extract(terrain_stack30cm$elevation, coordinates(soil_0_30cm_shp)[,1:2], buffer=1, fun=mean, na.rm=TRUE)
plot(soil_0_30cm_df$elevation, soil_0_30cm_df$elevation_1m)
summary(lm(soil_0_30cm_df$elevation_1m ~ soil_0_30cm_df$elevation)) #r2=0.9996
soil_0_30cm_df$curvature_mean_1m <- extract(terrain_stack30cm$curvature_mean, coordinates(soil_0_30cm_shp)[,1:2], buffer=1, fun=mean, na.rm=TRUE)
hist(soil_0_30cm_df$curvature_mean_1m)
plot(soil_0_30cm_df$curvature_mean, soil_0_30cm_df$curvature_mean_1m)
summary(lm(soil_0_30cm_df$curvature_mean_1m ~ soil_0_30cm_df$curvature_mean)) #r2=0.02347; just as before, very little correspondence
soil_0_30cm_df$slope_1m <- extract(terrain_stack30cm$slope, coordinates(soil_0_30cm_shp)[,1:2], buffer=1, fun=mean, na.rm=TRUE)
hist(soil_0_30cm_df$slope_1m)
plot(soil_0_30cm_df$slope, soil_0_30cm_df$slope_1m)
text(soil_0_30cm_df$slope, soil_0_30cm_df$slope_1m, labels=soil_0_30cm_df$point_no)
summary(lm(soil_0_30cm_df$slope_1m ~ soil_0_30cm_df$slope)) #r2=0.9258
soil_0_30cm_df$annual_khw.m2_1m <- extract(terrain_stack30cm$annual_kwh.m2, coordinates(soil_0_30cm_shp)[,1:2], buffer=1, fun=mean, na.rm=TRUE) / 1000 #exact same
hist(soil_0_30cm_df$annual_khw.m2_1m)
plot(soil_0_30cm_df$annual_khw.m2, soil_0_30cm_df$annual_khw.m2_1m)
text(soil_0_30cm_df$annual_khw.m2, soil_0_30cm_df$annual_khw.m2_1m, labels=soil_0_30cm_df$point_no)
summary(lm(soil_0_30cm_df$annual_khw.m2_1m ~ soil_0_30cm_df$annual_khw.m2))

#check relationship with curvature
summary(lm(soil_0_30cm_df$kgOrgC.m2 ~ soil_0_30cm_df$curvature_mean_1m)) #NS

mean_curv_9m <- raster(list.files(file.path(terrainDir, 'filtered_Hogan', '9m_filtered'), full.names = TRUE))
soil_0_30cm_df$curvature_mean_9m <- extract(mean_curv_9m, coordinates(soil_0_30cm_shp)[,1:2])
hist(soil_0_30cm_df$curvature_mean_9m)
plot(soil_0_30cm_df$curvature_mean, soil_0_30cm_df$curvature_mean_9m)
summary(lm(soil_0_30cm_df$curvature_mean_9m ~ soil_0_30cm_df$curvature_mean)) #r2=0.891 ; good correspondence between 3m and 9m
summary(lm(soil_0_30cm_df$kgOrgC.m2 ~ soil_0_30cm_df$curvature_mean_9m)) #same r2 as for 3m
