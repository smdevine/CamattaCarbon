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
calculate_thresholds <- function(x, iterations) {
  less_than_20 <- sum(x > mean(soil_0_30cm_shp$kgOrgC.m2)*0.8 & x < mean(soil_0_30cm_shp$kgOrgC.m2)*1.2) / iterations
  less_than_10 <- sum(x > mean(soil_0_30cm_shp$kgOrgC.m2)*0.9 & x < mean(soil_0_30cm_shp$kgOrgC.m2)*1.1) / iterations
  less_than_5 <- sum(x > mean(soil_0_30cm_shp$kgOrgC.m2)*0.95 & x < mean(soil_0_30cm_shp$kgOrgC.m2)*1.05) / iterations
  c(less_than_20, less_than_10, less_than_5)
}
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
#forage_terrain_energy <- read.csv(file.path(ResultsDir, 'tables', 'forage_terrain_energy_3m_final.csv'), stringsAsFactors = FALSE)
#list.files(file.path(soilCresults, 'shapefiles'))

library(raster)
library(car)
library(gstat)
library(spdep)
library(automap)
library(fpc)
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
soil_0_30cm_shp$kgOrgC.m2_0_10cm <- soil_0_10cm_shp$kgOrgC.m2
soil_0_30cm_shp$kgOrgC.m2_10_30cm <- soil_10_30cm_shp$kgOrgC.m2

corr_test <- function(x, df, y, mtd) {
  test <- cor.test(x, df[[y]], method = mtd)
  result <- data.frame(col.1=test$p.value, col.2=test$estimate)
  colnames(result) <- c(paste0(y, '.p.val.', mtd), paste0(y, if(mtd=='pearson') {'.tau.'} else {'.rho.'}, mtd))
  result
}
hist(soil_0_30cm_shp$kgOrgC.m2)
hist(soil_0_30cm_shp$kgIC.m2)
summary(soil_0_30cm_shp$kgIC.m2)
IC_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), corr_test, df=soil_0_30cm_shp, y='kgIC.m2', mtd='pearson'))
IC_corrs

IC_rank_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), corr_test, df=soil_0_30cm_shp, y='kgIC.m2', mtd='spearman'))
IC_rank_corrs

summary(lm(kgIC.m2 ~ annual_kwh.m2, data = soil_0_30cm_shp)) #r2=0.1
summary(lm(kgIC.m2 ~ curvature_mean + slope + elevation + annual_kwh.m2, data = soil_0_30cm_shp)) #r2=0.15
summary(lm(kgIC.m2 ~ curvature_mean + annual_kwh.m2, data = soil_0_30cm_shp)) #r2=0.14
summary(lm(kgIC.m2 ~ NDVI_2017mean_1m + annual_kwh.m2, data = soil_0_30cm_shp)) #r2=0.13
summary(lm(kgIC.m2 ~ NDVI_2018mean_1m + NDVI_2017mean_1m + elevation + annual_kwh.m2, data = soil_0_30cm_shp))
kgIC.m2_2var.lm <- lm(kgIC.m2 ~ NDVI_2017mean_1m + annual_kwh.m2, data = soil_0_30cm_shp)
plot(kgIC.m2_2var.lm$fitted.values, soil_0_30cm_shp$kgIC.m2)

kgIC.m2_4var.lm <- lm(kgIC.m2 ~ NDVI_2018mean_1m + NDVI_2017mean_1m + elevation + annual_kwh.m2, data = soil_0_30cm_shp)
plot(kgIC.m2_4var.lm$fitted.values, soil_0_30cm_shp$kgIC.m2)
