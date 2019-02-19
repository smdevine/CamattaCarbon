#intermed. results on line 370
#corrected SIC, TN, and P contents for 10-30 cm soils by multiplying by 2 (1/8/19)
#changed P contents to reflect Olsen (HCO3-P) extracts (1/8/19)
##TO-DO
#soil P vs. biomass
#soil N maps even though highly correlated with soil C
#soil N vs. biomass
#soil N, P, and C vs. biomass
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
kf <- kfold(1:105, k=20) #was 20 for all results in v3 of Chp 3


#read-in 0-30 cm and set-up k-fold tags [go to line 370 to read-in 0-30cm data now]
#soil_0_30cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
soil_0_30cm_shp$energy_colors <- ifelse(soil_0_30cm_shp$annual_kwh.m2 <= 1200, 'blue', ifelse(soil_0_30cm_shp$annual_kwh.m2 > 1200 & soil_0_30cm_shp$annual_kwh.m2 < 1410, 'orange2', 'red3'))
soil_0_30cm_shp$kgTC.m2 <- soil_0_30cm_shp$kgOrgC.m2 + soil_0_30cm_shp$kgIC.m2
#plot(soil_0_30cm_shp, cex=soil_0_30cm_shp$kgOrgC.m2/2, pch=20)
sd(soil_0_30cm_shp$kgOrgC.m2)/mean(soil_0_30cm_shp$kgOrgC.m2) #CV is 19.5%
soil_0_30cm_shp$aspect_class <- ifelse(soil_0_30cm_shp$aspect >= 45 & soil_0_30cm_shp$aspect < 135, 'east', ifelse(soil_0_30cm_shp$aspect >= 135 & soil_0_30cm_shp$aspect < 225, 'south', ifelse(soil_0_30cm_shp$aspect >= 225 & soil_0_30cm_shp$aspect < 315, 'west', 'north')))
soil_0_30cm_shp$aspect_class_sp <- ifelse(soil_0_30cm_shp$aspect >= 67.5 & soil_0_30cm_shp$aspect < 112.5, 'east', ifelse(soil_0_30cm_shp$aspect >= 112.5 & soil_0_30cm_shp$aspect < 157.5, 'southeast', ifelse(soil_0_30cm_shp$aspect >= 157.5 & soil_0_30cm_shp$aspect < 202.5, 'south', ifelse(soil_0_30cm_shp$aspect >= 202.5 & soil_0_30cm_shp$aspect < 247.5, 'southwest', ifelse(soil_0_30cm_shp$aspect >= 247.5 & soil_0_30cm_shp$aspect < 292.5, 'west', ifelse(soil_0_30cm_shp$aspect >= 292.5 & soil_0_30cm_shp$aspect < 337.5, 'northwest', ifelse(soil_0_30cm_shp$aspect >= 337.5 & soil_0_30cm_shp$aspect < 22.5, 'north', 'northeast')))))))

#0-10 dataset (modified orgC, TN, clay, IC, and P colnames to match naming conventions for 0-30)
soil_0_10cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
soil_0_10cm_shp$kgTC.m2 <- soil_0_10cm_shp$kgOrgC.m2 + soil_0_10cm_shp$kgIC.m2
#plot(soil_0_10cm_shp, cex=soil_0_10cm_shp$kgOrgC.m2/1.5, pch=20)

#10-30 dataset (modified orgC, TN, clay, IC, and P colnames to match naming conventions for 0-30)
soil_10_30cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
soil_10_30cm_shp$kgTC.m2 <- soil_10_30cm_shp$kgOrgC.m2 + soil_10_30cm_shp$kgIC.m2
#plot(soil_10_30cm_shp, cex=soil_10_30cm_shp$kgOrgC.m2/1.5, pch=20)

#get some summary stats
#check weighted average 0-30 cm %org C
sum(soil_0_10cm_shp$point_no - soil_10_30cm_shp$point_no) #check again that location order is the same in each dataset
sum(soil_0_10cm_shp$point_no - soil_0_30cm_shp$point_no)
orgC_0_30_percent <- soil_0_10cm_shp$orgC.percent * 1/3 + soil_10_30cm_shp$orgC.percent * 2/3
summary(orgC_0_30_percent)
hist(orgC_0_30_percent) #mean is 0.94%
summary(soil_0_10cm_shp$orgC.percent)
sd(soil_0_10cm_shp$orgC.percent)
summary(soil_10_30cm_shp$orgC.percent)
sd(soil_10_30cm_shp$orgC.percent)
hist(soil_0_10cm_shp$orgC.percent)
hist(soil_10_30cm_shp$orgC.percent)

DB_0_30 <- soil_0_10cm_shp$bulk_density_g_cm3 * 1/3 + soil_10_30cm_shp$bulk_density_g_cm3 * 2/3
summary(DB_0_30)
summary(soil_0_10cm_shp$bulk_density_g_cm3)
summary(soil_10_30cm_shp$bulk_density_g_cm3)

summary(soil_0_30cm_shp$kgOrgC.m2)
sd(soil_0_30cm_shp$kgOrgC.m2)
summary(soil_0_10cm_shp$kgOrgC.m2)
summary(soil_10_30cm_shp$kgOrgC.m2)
summary(soil_0_10cm_shp$kgOrgC.m2 / soil_0_30cm_shp$kgOrgC.m2)
sum(soil_0_30cm_shp$kgOrgC.m2 > mean(soil_0_30cm_shp$kgOrgC.m2)*0.8 & soil_0_30cm_shp$kgOrgC.m2 < mean(soil_0_30cm_shp$kgOrgC.m2)*1.2) #75 of 105 samples are within 20% of the mean
test <- sample(1:105, 104)
test[order(test)]
calc_rand_mean <- function(n) {
  mean(soil_0_30cm_shp$kgOrgC.m2[sample(1:105, n)])
}
replicate(3, calc_rand_mean(3))
lapply(1:3, FUN=calc_rand_mean, 3)
iterate_rand_mean <- function(iterations, n) {
  replicate(iterations, calc_rand_mean(n))
}
sample_1 <- iterate_rand_mean(10000, 1)
sample_2 <- iterate_rand_mean(10000, 2)
sample_3 <- iterate_rand_mean(10000, 3)
sample_4 <- iterate_rand_mean(10000, 4)
sample_5 <- iterate_rand_mean(10000, 5)
sample_6 <- iterate_rand_mean(10000, 6)
sample_7 <- iterate_rand_mean(10000, 7)
sample_8 <- iterate_rand_mean(10000, 8)
sample_9 <- iterate_rand_mean(10000, 9)
sample_10 <- iterate_rand_mean(10000, 10)
sample_12 <- iterate_rand_mean(10000, 12)
sample_15 <- iterate_rand_mean(10000, 15)
sample_18 <- iterate_rand_mean(10000, 18)
sample_20 <- iterate_rand_mean(10000, 20)
sample_21 <- iterate_rand_mean(10000, 21)
sample_24 <- iterate_rand_mean(10000, 24)
sample_25 <- iterate_rand_mean(10000, 25)
sample_27 <- iterate_rand_mean(10000, 27)
sample_29 <- iterate_rand_mean(10000, 29)
sample_30 <- iterate_rand_mean(10000, 30)
#sample_35 <- iterate_rand_mean(10000, 35)

calculate_thresholds(sample_2, 10000) #0.8594 0.5556 0.2964
calculate_thresholds(sample_3, 10000) #0.9324 0.6350 0.3484
calculate_thresholds(sample_4, 10000) #0.9639 0.7072 0.4031
calculate_thresholds(sample_5, 10000) #0.9826 0.7614 0.4494
calculate_thresholds(sample_6, 10000) #0.9925 0.8102 0.4864
calculate_thresholds(sample_7, 10000) #0.9953 0.8434 0.5235
calculate_thresholds(sample_8, 10000) #0.9982 0.8755 0.5589
calculate_thresholds(sample_9, 10000) #0.9991 0.8982 0.5822
calculate_thresholds(sample_10, 10000) #random approach: 0.9992 0.9119 0.6017 
calculate_thresholds(sample_15, 10000) #1.0000 0.9737 0.7284
calculate_thresholds(sample_20, 10000) #1.0000 0.9894 0.7996
calculate_thresholds(sample_25, 10000) #1.0000 0.9980 0.8657
calculate_thresholds(sample_29, 10000) #1.0000 0.9993 0.9008
calculate_thresholds(sample_30, 10000) #1.0000 0.9996 0.9091
calculate_thresholds(sample_35, 10000) #1.000 1.000 0.938


#get summary stats for soil properties
summary(soil_0_30cm_shp$kgIC.m2)
summary(soil_0_10cm_shp$kgIC.m2)
summary(soil_10_30cm_shp$kgIC.m2)
summary(soil_0_10cm_shp$inorgC.percent)
sum(soil_0_10cm_shp$inorgC.percent==0)
sd(soil_0_10cm_shp$inorgC.percent)
summary(soil_10_30cm_shp$inorgC.percent)
sd(soil_10_30cm_shp$inorgC.percent)
sum(soil_10_30cm_shp$inorgC.percent==0)
summary(soil_0_30cm_shp$kgIC.m2 / soil_0_30cm_shp$kgTC.m2)
summary(soil_0_10cm_shp$kgIC.m2 / soil_0_10cm_shp$kgTC.m2)
summary(soil_10_30cm_shp$kgIC.m2 / soil_10_30cm_shp$kgTC.m2)
summary(soil_0_30cm_shp$kgTC.m2)
summary(soil_0_10cm_shp$kgTC.m2)
summary(soil_10_30cm_shp$kgTC.m2)

summary(soil_0_30cm_shp$clay_wtd)
summary(soil_0_10cm_shp$CLAY)
sd(soil_0_10cm_shp$CLAY)
summary(soil_10_30cm_shp$CLAY)
sd(soil_10_30cm_shp$CLAY)
summary(soil_0_10cm_shp$WMPD_mm)
summary(soil_10_30cm_shp$WMPD_mm)
summary(as.factor(soil_0_10cm_shp$TEXTURE)) / 105
summary(as.factor(soil_10_30cm_shp$TEXTURE)) / 105

summary(soil_0_30cm_shp$gP.m2)
summary(soil_0_10cm_shp$gP.m2)
summary(soil_10_30cm_shp$gP.m2)
summary(soil_0_10cm_shp$HCO3_P)
sd(soil_0_10cm_shp$HCO3_P)
sum(soil_0_10cm_shp$HCO3_P < 3) #0 very low
sum(soil_0_10cm_shp$HCO3_P >= 3 & soil_0_10cm_shp$HCO3_P < 7) #37 (35%)low
sum(soil_0_10cm_shp$HCO3_P >= 7 & soil_0_10cm_shp$HCO3_P < 13)#57 (54%) medium
sum(soil_0_10cm_shp$HCO3_P >= 13 & soil_0_10cm_shp$HCO3_P < 22) #8 (8%) high
sum(soil_0_10cm_shp$HCO3_P >= 22) #3 (3%)

summary(soil_10_30cm_shp$HCO3_P)
sd(soil_10_30cm_shp$HCO3_P)
sum(soil_10_30cm_shp$HCO3_P < 3) #31 (30%) very low
sum(soil_10_30cm_shp$HCO3_P >= 3 & soil_10_30cm_shp$HCO3_P < 7) #71 (68%)low
sum(soil_10_30cm_shp$HCO3_P >= 7 & soil_10_30cm_shp$HCO3_P < 13)#2 (2%) medium
sum(soil_10_30cm_shp$HCO3_P >= 13 & soil_10_30cm_shp$HCO3_P < 22) #0 (0%) high
sum(soil_10_30cm_shp$HCO3_P >= 22) #1 (1%)

#read-in terrain properties
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan', '3m_filtered'), full.names = TRUE))
names(Mar2017_terrain_3m)
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
#plot(solrad_raster)
#plot(soil_0_30cm_shp, pch=1, add=TRUE)
Mar2017_terrain_3m$annual_kwh.m2 <- solrad_raster
#plot(Mar2017_terrain_3m$curvature_mean)
#plot(soil_0_30cm_shp, pch=1, cex=soil_0_30cm_shp$kgOrgC.m2/2, add=TRUE)
#plot(Mar2017_terrain_3m$elevation)

#read-in NDVI rasters
#list.files(NDVIDir)
NDVI_stack <- stack(list.files(NDVIDir, full.names = TRUE))
#names(NDVI_stack)
#class(NDVI_stack)
cellStats(NDVI_stack, 'mean')
NDVI_2017 <- NDVI_stack[[c(1,3,5,7,9)]]
NDVI_2018 <- NDVI_stack[[c(2,4,6,8)]]
maxNDVI_2017 <- calc(NDVI_2017, fun=max)
meanNDVI_2017 <- calc(NDVI_2017, fun=mean)#, filename=file.path(FiguresDir, 'meanNDVI2017.tif'))
#maxNDVI_2017_1m <- aggregate(maxNDVI_2017, fact=3, fun='mean')
#meanNDVI_2017_1m <- aggregate(meanNDVI_2017, fact=3, fun='mean')
#maxNDVI_2017_2m <- aggregate(maxNDVI_2017, fact=6, fun='mean')
maxNDVI_2018 <- calc(NDVI_2018, fun=max)
meanNDVI_2018 <- calc(NDVI_2018, fun=mean)#, filename=file.path(FiguresDir, 'meanNDVI2018.tif'))
#maxNDVI_2018_1m <- aggregate(maxNDVI_2018, fact=3, fun='mean')
#meanNDVI_2018_1m <- aggregate(meanNDVI_2018, fact=3, fun='mean')
#maxNDVI_2018_2m <- aggregate(maxNDVI_2018, fact=6, fun='mean')
#plot(maxNDVI_2017)
#plot(maxNDVI_2018)
#plot(maxNDVI_2017_1m)
#plot(maxNDVI_2018_1m)

#compare NDVI to SOC
soil_0_30cm_shp$NDVI_2017max <- extract(maxNDVI_2017, soil_0_30cm_shp)
soil_0_30cm_shp$NDVI_2017mean <- extract(meanNDVI_2017, soil_0_30cm_shp)
soil_0_30cm_shp$NDVI_2017max_1m <- extract(maxNDVI_2017, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_2017mean_1m <- extract(meanNDVI_2017, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_2017max_2m <- extract(maxNDVI_2017, soil_0_30cm_shp, buffer=2, fun=mean)
soil_0_30cm_shp$NDVI_2017mean_2m <- extract(meanNDVI_2017, soil_0_30cm_shp, buffer=2, fun=mean)
soil_0_30cm_shp$NDVI_2018max <- extract(maxNDVI_2018, soil_0_30cm_shp)
soil_0_30cm_shp$NDVI_2018mean <- extract(meanNDVI_2018, soil_0_30cm_shp)
soil_0_30cm_shp$NDVI_2018max_1m <- extract(maxNDVI_2018, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_2018mean_1m <- extract(meanNDVI_2018, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_2018max_2m <- extract(maxNDVI_2018, soil_0_30cm_shp, buffer=2, fun=mean)
soil_0_30cm_shp$NDVI_2018mean_2m <- extract(meanNDVI_2018, soil_0_30cm_shp, buffer=2, fun=mean)
#soil_0_30cm_shp$NDVI_avg_1m <- rowMeans(as.data.frame(soil_0_30cm_shp)[,c('NDVI_2017max_1m', 'NDVI_2018max_1m')])
hist(soil_0_30cm_shp$NDVI_2017max)
hist(soil_0_30cm_shp$NDVI_2017max_1m)
hist(soil_0_30cm_shp$NDVI_2017max_2m)
hist(soil_0_30cm_shp$NDVI_2018max)
hist(soil_0_30cm_shp$NDVI_2018max_1m)
hist(soil_0_30cm_shp$NDVI_2018max_2m)
plot(soil_0_30cm_shp$NDVI_2017max_1m, soil_0_30cm_shp$kgOrgC.m2)
abline(lm(kgOrgC.m2 ~ NDVI_2017max_1m, data = soil_0_30cm_shp), lty=2)
plot(lm(kgOrgC.m2 ~ NDVI_2017max_1m, data = soil_0_30cm_shp))
lm_NDVI2017_OrgC <- lm(kgOrgC.m2 ~ NDVI_2017max_1m, data = soil_0_30cm_shp)
summary(lm_NDVI2017_OrgC) #r2=0.25

plot(soil_0_30cm_shp$NDVI_2018max_1m, soil_0_30cm_shp$kgOrgC.m2) 
abline(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_0_30cm_shp), lty=2)

summary(lm(kgOrgC.m2 ~ NDVI_2017max, data = soil_0_30cm_shp)) #r^2=0.17
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m, data = soil_0_30cm_shp)) #r^2=0.25
summary(lm(kgOrgC.m2 ~ NDVI_2017max_2m, data = soil_0_30cm_shp)) #r^2=0.26
summary(lm(kgOrgC.m2 ~ NDVI_2018max, data=soil_0_30cm_shp)) #no relationship: r^2 < 0.01; p-val=0.55
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data=soil_0_30cm_shp)) #no relationship: r^2 < 0.01; p-val=0.65
summary(lm(kgOrgC.m2 ~ NDVI_2017mean, data=soil_0_30cm_shp)) #r2=0.27
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m, data = soil_0_30cm_shp)) #r2=0.31
summary(lm(kgOrgC.m2 ~ NDVI_2018mean_1m, data = soil_0_30cm_shp)) #NS
plot(soil_0_30cm_shp$elevation, soil_0_30cm_shp$NDVI_2017mean_1m)
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation, data = soil_0_30cm_shp)) #r^2=0.47; 5-var model, all params significant; adj r2=0.44
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + NDVI_2018max_1m + curvature_mean + annual_kwh.m2 + slope + elevation, data = soil_0_30cm_shp)) #also r2=0.47 with dry year NDVI NS
vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation, data = soil_0_30cm_shp)) #r2=0.5; NS: NDVI_2018mean_1m
vif(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2, data = soil_0_30cm_shp))


#add NDVI to 0-10 and 10-30 separate datasets
soil_0_10cm_shp$NDVI_2017max_1m <- extract(maxNDVI_2017, soil_0_10cm_shp, buffer=1, fun=mean)
soil_0_10cm_shp$NDVI_2018max_1m <- extract(maxNDVI_2018, soil_0_10cm_shp, buffer=1, fun=mean)
soil_0_10cm_shp$NDVI_2017mean_1m <- extract(meanNDVI_2017, soil_0_10cm_shp, buffer=1, fun=mean)
soil_0_10cm_shp$NDVI_2018mean_1m <- extract(meanNDVI_2018, soil_0_10cm_shp, buffer=1, fun=mean)
soil_10_30cm_shp$NDVI_2017max_1m <- extract(maxNDVI_2017, soil_10_30cm_shp, buffer=1, fun=mean)
soil_10_30cm_shp$NDVI_2018max_1m <- extract(maxNDVI_2018, soil_10_30cm_shp, buffer=1, fun=mean)
soil_10_30cm_shp$NDVI_2017mean_1m <- extract(meanNDVI_2017, soil_10_30cm_shp, buffer=1, fun=mean)
soil_10_30cm_shp$NDVI_2018mean_1m <- extract(meanNDVI_2018, soil_10_30cm_shp, buffer=1, fun=mean)

#read-in reflectance data from 2016
# list.files(ReflectanceDir)
# NIR_Nov2016 <- raster(file.path(ReflectanceDir, 'Nov2016', 'Camatta_11112016_Topoed_nir.tif'))
# NIR_Nov2016_1m <- aggregate(NIR_Nov2016, fact=3, fun='mean')
# soil_0_30cm_shp$NIR_Nov2016_1m <- extract(NIR_Nov2016_1m, soil_0_30cm_shp)
# plot(soil_0_30cm_shp$NIR_Nov2016_1m, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(kgOrgC.m2 ~ NIR_Nov2016_1m, data = soil_0_30cm_shp))
# Red_Nov2016 <- raster(file.path(ReflectanceDir, 'Camatta_11112016_Topoed_red.tif'))
# Red_Nov2016_1m <- aggregate(Red_Nov2016, fact=3, fun='mean')
# soil_0_30cm_shp$Red_Nov2016_1m <- extract(Red_Nov2016_1m, soil_0_30cm_shp)
# plot(soil_0_30cm_shp$Red_Nov2016_1m, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(kgOrgC.m2 ~ Red_Nov2016_1m, data = soil_0_30cm_shp))
# Rededge_Nov2016 <- raster(file.path(ReflectanceDir, 'Camatta_11112016_Topoed_rededge.tif'))
# Rededge_Nov2016_1m <- aggregate(Rededge_Nov2016, fact=3, fun='mean')
# soil_0_30cm_shp$Rededge_Nov2016_1m <- extract(Rededge_Nov2016_1m, soil_0_30cm_shp)
# plot(soil_0_30cm_shp$Rededge_Nov2016_1m, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(kgOrgC.m2 ~ Rededge_Nov2016_1m, data = soil_0_30cm_shp))
# summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_Nov2016_1m + Red_Nov2016_1m + Rededge_Nov2016_1m, data = soil_0_30cm_shp))
# vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_Nov2016_1m + Red_Nov2016_1m + Rededge_Nov2016_1m, data = soil_0_30cm_shp))

#read-in reflectance from Jan 2017
# list.files(file.path(ReflectanceDir, 'Jan2017'))
# NIR_Jan2017 <- raster(file.path(ReflectanceDir, 'Jan2017', 'Camatta_01162017_Topoed_nir.tif'))
# NIR_Jan2017_1m <- aggregate(NIR_Jan2017, fact=3, fun='mean')
# soil_0_30cm_shp$NIR_Jan2017_1m <- extract(NIR_Jan2017_1m, soil_0_30cm_shp)
# plot(soil_0_30cm_shp$NIR_Jan2017_1m, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(kgOrgC.m2 ~ NIR_Jan2017_1m, data = soil_0_30cm_shp))
# Red_Jan2017 <- raster(file.path(ReflectanceDir, 'Jan2017', 'Camatta_01162017_Topoed_red.tif'))
# Red_Jan2017_1m <- aggregate(Red_Jan2017, fact=3, fun='mean')
# soil_0_30cm_shp$Red_Jan2017_1m <- extract(Red_Jan2017_1m, soil_0_30cm_shp)
# plot(soil_0_30cm_shp$Red_Jan2017_1m, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(kgOrgC.m2 ~ Red_Jan2017_1m, data = soil_0_30cm_shp))
# Rededge_Jan2017 <- raster(file.path(ReflectanceDir, 'Camatta_11112016_Topoed_rededge.tif'))
# Rededge_Jan2017_1m <- aggregate(Rededge_Jan2017, fact=3, fun='mean')
# soil_0_30cm_shp$Rededge_Jan2017_1m <- extract(Rededge_Jan2017_1m, soil_0_30cm_shp)
# plot(soil_0_30cm_shp$Rededge_Jan2017_1m, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(kgOrgC.m2 ~ Rededge_Jan2017_1m, data = soil_0_30cm_shp))
# summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_Jan2017_1m + Red_Jan2017_1m + Rededge_Jan2017_1m, data = soil_0_30cm_shp))
# vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_Jan2017_1m + Red_Jan2017_1m + Rededge_Jan2017_1m, data = soil_0_30cm_shp))
# summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Rededge_Jan2017_1m, data = soil_0_30cm_shp))

#function to explore and/or add reflectance data to soils data
add_reflectance_to_df <- function(df, Dirname, fname_date, buffer) {
  NIR <- raster(file.path(ReflectanceDir, Dirname, paste0('Camatta_', fname_date ,'_Topoed_nir.tif')))
  #NIR <- aggregate(NIR, fact=agglevel, fun='mean')
  df[[paste0('NIR_', Dirname)]] <- extract(NIR, df, buffer=buffer, fun=mean)
  plot(df[[paste0('NIR_', Dirname)]], df$kgOrgC.m2, main=paste('NIR band in', Dirname))
  #as.formula(paste(varname, '~ curvature_mean + slope + annual_kwh.m2 + WMPD_mm'))
  print(summary(lm(as.formula(paste0('kgOrgC.m2 ~ NIR_', Dirname)), data = df)))
  Red <- raster(file.path(ReflectanceDir, Dirname, paste0('Camatta_', fname_date , '_Topoed_red.tif')))
  #Red <- aggregate(Red, fact=agglevel, fun='mean')
  df[[paste0('Red_', Dirname)]] <- extract(Red, df, buffer=buffer, fun=mean)
  plot(df[[paste0('Red_', Dirname)]], df$kgOrgC.m2, main=paste('Red band in', Dirname))
  print(summary(lm(as.formula(paste0('kgOrgC.m2 ~ Red_', Dirname)), data = df)))
  Rededge <- raster(file.path(ReflectanceDir, Dirname, paste0('Camatta_', fname_date, '_Topoed_rededge.tif')))
  #Rededge <- aggregate(Rededge, fact=agglevel, fun='mean')
  df[[paste0('Rededge_', Dirname)]] <- extract(Rededge, df, buffer=buffer, fun=mean)
  plot(df[[paste0('Rededge_', Dirname)]], df$kgOrgC.m2, main=paste('Rededge band in', Dirname))
  print(summary(lm(as.formula(paste0('kgOrgC.m2 ~ Rededge_', Dirname)), data = df)))
  print(summary(lm(as.formula(paste0('kgOrgC.m2 ~ Rededge_', Dirname, ' + Red_', Dirname, ' + NIR_', Dirname, ' + NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation')), data = df)))
  df
}

#add 2017 growing season reflectance data
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Nov2016', fname_date = '11112016', buffer = 1) #r2=0.50 all reflectance bands NS
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Jan2017', fname_date = '01162017', buffer = 1) #Multiple R-squared:  0.50; direct sig. association with Red Band (p=0.03)
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Feb2017', fname_date = '02152017', buffer = 1) #Multiple R-squared: 0.52 ; all Red bands NS in MLR but direct sig. association with Red Band (p<0.001; r2=0.1)
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Mar2017', fname_date = '03172017', buffer = 1) #Multiple R-squared:  0.51; all Red bands NS in MLR and  direct sig. association with Red Band (p<0.001; r2=0.19) and NIR (p=0.008); NDVI also now non-sig (p=0.1)
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Apr2017', fname_date = '04062017', buffer = 1) #Multiple R-squared:  0.50; all Red bands NS in MLR and  direct sig. association with Red Band (p<0.001; r2=0.23) and NIR (p<0.001; r2=0.11); NDVI also now non-sig (p=0.07)
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'May2017', fname_date = '04302017', buffer = 1) #Multiple R-squared:  0.53; all Red bands NS in MLR and  direct sig. association with Red Band (p<0.001; r2=0.17); NDVI also now non-sig (p=0.15)
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Jan2018', fname_date = '01162018', buffer = 1) #Multiple R-squared:  0.51 (all Red band p-val NS in MLR); NIR_May_2017  direct association with NIR and Rededge Band (p<0.001); NDVI back to significance
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Feb2018', fname_date = '02152018', buffer = 1) #Multiple R-squared:  0.50 (all Red band p-val NS); direct association with Red bands all sig; NDVI back to significance; reflectance numbers very low in Feb2018
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Mar2018', fname_date = '03192018', buffer = 1) #Multiple R-squared: 0.51 (all Red band p-val NS); direct association with Red bands all sig; NDVI back to significance;
soil_0_30cm_shp <- add_reflectance_to_df(soil_0_30cm_shp, Dirname = 'Apr2018', fname_date = '04142018', buffer = 1) #Multiple R-squared: 0.50   (all Red band p-val NS); direct sig. association with Red and Rededge bands are sig; NDVI back to significance; numbers very low in Feb2018


hist(soil_0_30cm_shp$NIR_Jan2017)
hist(soil_0_30cm_shp$NIR_Feb2017) #wacky
hist(soil_0_30cm_shp$NIR_Mar2017)
hist(soil_0_30cm_shp$NIR_Apr2017)
hist(soil_0_30cm_shp$NIR_May2017)
hist(soil_0_30cm_shp$Red_Jan2017)
hist(soil_0_30cm_shp$Red_Feb2017) #wacky
hist(soil_0_30cm_shp$Red_Mar2017)
hist(soil_0_30cm_shp$Red_Apr2017)
hist(soil_0_30cm_shp$Red_May2017)
hist(soil_0_30cm_shp$NIR_Jan2018)
hist(soil_0_30cm_shp$NIR_Feb2018)#wacky
hist(soil_0_30cm_shp$NIR_Mar2018)#wacky
hist(soil_0_30cm_shp$NIR_Apr2018)
hist(soil_0_30cm_shp$Red_Jan2018)
hist(soil_0_30cm_shp$Red_Feb2018)#wacky
hist(soil_0_30cm_shp$Red_Mar2018)#wacky
hist(soil_0_30cm_shp$Red_Apr2018)

#get mean Red 2017 growing season, excluding Feb
soil_0_30cm_shp$Red_meanGS2017 <- apply(as.data.frame(soil_0_30cm_shp)[,c('Red_Jan2017', 'Red_Mar2017', 'Red_Apr2017', 'Red_May2017')], 1, FUN = mean) #excludes Feb2017 which were wonky
soil_0_30cm_shp$NIR_meanGS2017 <- apply(as.data.frame(soil_0_30cm_shp)[,c('NIR_Jan2017', 'NIR_Mar2017', 'NIR_Apr2017', 'NIR_May2017')], 1, FUN = mean)
soil_0_30cm_shp$Red_meanGS2018 <- apply(as.data.frame(soil_0_30cm_shp)[,c('Red_Jan2018', 'Red_Apr2018')], 1, FUN = mean) #excludes Feb2018 and Mar2018 which were wonky
soil_0_30cm_shp$NIR_meanGS2018 <- apply(as.data.frame(soil_0_30cm_shp)[,c('NIR_Jan2018', 'NIR_Apr2018')], 1, FUN = mean) #excludes Feb2018 and Mar2018 which were wonky

#add these columns to 0-10 and 10-30 subsets
addColumn <- function(df, df_input, column_names) {
  for (i in seq_along(column_names)) {
    df[[column_names[i]]] <- df_input[[column_names[i]]]
  }
  df
}
soil_0_10cm_shp <- addColumn(soil_0_10cm_shp, soil_0_30cm_shp, c('Red_meanGS2017', 'Red_meanGS2018', 'Red_Nov2016', 'Red_May2017', 'Red_Jan2018', 'Red_Apr2018', 'NIR_meanGS2017', 'NIR_meanGS2018', 'NIR_Nov2016', 'NIR_May2017', 'NIR_Jan2018', 'NIR_Apr2018'))
soil_10_30cm_shp <- addColumn(soil_10_30cm_shp, soil_0_30cm_shp, c('Red_meanGS2017', 'Red_meanGS2018', 'Red_Nov2016', 'Red_May2017', 'Red_Jan2018', 'Red_Apr2018', 'NIR_meanGS2017', 'NIR_meanGS2018', 'NIR_Nov2016', 'NIR_May2017', 'NIR_Jan2018', 'NIR_Apr2018'))

#save intermediate results
#0-30 cm
shapefile(soil_0_30cm_shp, file.path(soilCresults, 'intermediate_results', 'soil_0_30cm.shp'), overwrite=TRUE)
write.csv(as.data.frame(soil_0_30cm_shp), file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), row.names=FALSE)

#0-10 cm
shapefile(soil_0_10cm_shp, file.path(soilCresults, 'intermediate_results', 'soil_0_10cm.shp'), overwrite=TRUE)
write.csv(as.data.frame(soil_0_10cm_shp), file.path(soilCresults, 'intermediate_results', 'soil_0_10cm_df.csv'), row.names=FALSE)

#10-30 cm
shapefile(soil_10_30cm_shp, file.path(soilCresults, 'intermediate_results', 'soil_10_30cm.shp'), overwrite=TRUE)
write.csv(as.data.frame(soil_10_30cm_shp), file.path(soilCresults, 'intermediate_results', 'soil_10_30cm_df.csv'), row.names=FALSE)

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

#create summary stats for 105 pts
summary_stats_0_30cm <- as.data.frame(lapply(as.data.frame(soil_0_30cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'clay_wtd', 'sand_wtd', 'silt_wtd', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'elevation', 'slope', 'annual_kwh.m2', 'curvature_mean', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_0_30cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_0_30cm <- summary_stats_0_30cm[,c(ncol(summary_stats_0_30cm), 1:(ncol(summary_stats_0_30cm)-1))]
write.csv(summary_stats_0_30cm, file.path(CarbonDir, 'summaries', 'tables', 'summary_stats_0_30cm.csv'), row.names = FALSE)

summary_stats_0_10cm <- as.data.frame(lapply(as.data.frame(soil_0_10cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_0_10cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_0_10cm <- summary_stats_0_10cm[,c(ncol(summary_stats_0_10cm), 1:(ncol(summary_stats_0_10cm)-1))]
write.csv(summary_stats_0_10cm, file.path(CarbonDir, 'summaries', 'tables', 'summary_stats_0_10cm.csv'), row.names = FALSE)

summary_stats_10_30cm <- as.data.frame(lapply(as.data.frame(soil_10_30cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_10_30cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_10_30cm <- summary_stats_10_30cm[,c(ncol(summary_stats_10_30cm), 1:(ncol(summary_stats_10_30cm)-1))]
write.csv(summary_stats_10_30cm, file.path(CarbonDir, 'summaries', 'tables', 'summary_stats_10_30cm.csv'), row.names = FALSE)

#create grid for spatial predictions and upscaling NDVI & reflectance data to same grid
r <- raster(Mar2017_terrain_3m)
#adjuster <- 30
e <- extent(meanNDVI_2017)
#e <- extent(c(xmin(soil_0_30cm_shp)- adjuster, xmax(soil_0_30cm_shp)+adjuster, ymin(soil_0_30cm_shp) - adjuster, ymax(soil_0_30cm_shp) + adjuster))
r <- crop(r, e)
#r <- raster(xmn=(xmin(soil_0_30cm_shp)-10), xmx=(xmax(soil_0_30cm_shp)+10), ymn=(ymin(soil_0_30cm_shp) - 10), ymx=(ymax(soil_0_30cm_shp) + 10), resolution=3, crs=crs(soil_0_30cm_shp))

summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_May2017 + Red_May2017 + Rededge_May2017, data = soil_0_30cm_shp)) #r2=0.53; adj r2=0.49 but VIF problems with NIR and Rededge
vif(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_May2017 + Red_May2017 + Rededge_May2017, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp)) #r2=0.51
vif(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp))
plot(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_May2017 + Rededge_May2017, data = soil_0_30cm_shp)) #r2=0.51
vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_May2017 + Rededge_May2017, data = soil_0_30cm_shp))

#test with means
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_meanGS2017 + NIR_meanGS2017, data = soil_0_30cm_shp)) #r2=0.50
vif(lm(kgOrgC.m2 ~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_meanGS2017 + NIR_meanGS2017, data = soil_0_30cm_shp)) #not good on this front

#test Nov2016
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_Nov2016 + Red_Nov2016 + Rededge_Nov2016, data = soil_0_30cm_shp)) #r2=0.48; VIF problems with NIR and Rededge again
vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_Nov2016 + Red_Nov2016 + Rededge_Nov2016, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + NIR_Nov2016, data = soil_0_30cm_shp)) #r2=0.48
vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + NIR_Nov2016, data = soil_0_30cm_shp))
plot(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + NIR_Nov2016, data = soil_0_30cm_shp)) #slight VIF with Red and NIR
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Rededge_Nov2016, data = soil_0_30cm_shp)) #r2=0.48
vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Rededge_Nov2016, data = soil_0_30cm_shp))

#combine both dates for reflectance
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp)) #NS (p>0.5): NIR_Nov2016 #r2=0.52
vif(lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp))
lm_orgC_8var <- lm(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp)
plot(lm_orgC_8var$fitted.values, soil_0_30cm_shp$kgOrgC.m2)
abline(0,1,lty=2)

#now do random stratification approach
#standardize predictors first 
soil_0_30cm_shp$curvature_mean_norm <- normalize_var(soil_0_30cm_shp$curvature_mean)
soil_0_30cm_shp$slope_norm <- normalize_var(soil_0_30cm_shp$slope)
soil_0_30cm_shp$annual_kwh.m2_norm <- normalize_var(soil_0_30cm_shp$annual_kwh.m2)
soil_0_30cm_shp$elevation_norm <- normalize_var(soil_0_30cm_shp$elevation)
soil_0_30cm_shp$NDVI_2017mean_1m_norm <- normalize_var(soil_0_30cm_shp$NDVI_2017mean_1m)
# soil_0_30cm_shp$NIR_meanGS2017_norm <- normalize_var(soil_0_30cm_shp$NIR_meanGS2017)
# soil_0_30cm_shp$Red_meanGS2017_norm <- normalize_var(soil_0_30cm_shp$Red_meanGS2017)
# soil_0_30cm_shp$Red_meanGS2017_norm

summary(lm(kgOrgC.m2 ~ curvature_mean_norm + slope_norm + annual_kwh.m2_norm + elevation_norm + NDVI_2017mean_1m_norm, data =  soil_0_30cm_shp[-c(2,5,82),]))
plot(lm(kgOrgC.m2 ~ curvature_mean_norm + slope_norm + annual_kwh.m2_norm + elevation_norm + NDVI_2017mean_1m_norm, data =  soil_0_30cm_shp[-c(2,5,82),]))

Mar2017_terrain_3m_cropped <- crop(Mar2017_terrain_3m, r)
#Mar2017_terrain_3m_cropped$maxNDVI_2017 <- resample(maxNDVI_2017, Mar2017_terrain_3m_cropped)
Mar2017_terrain_3m_cropped$NDVI_2017mean_1m <- resample(meanNDVI_2017, Mar2017_terrain_3m_cropped) #naming of NDVI is only to match soil_0_30cm
#soil_0_30cm_shp$meanNDVI2017test <- extract(Mar2017_terrain_3m_cropped$NDVI_2017mean_1m, soil_0_30cm_shp)
#summary(lm(NDVI_2017mean_1m ~ meanNDVI2017test, data = soil_0_30cm_shp)) #r2=0.90
#plot(soil_0_30cm_shp$meanNDVI2017test, soil_0_30cm_shp$NDVI_2017mean_1m)
plot(Mar2017_terrain_3m_cropped$NDVI_2017mean_1m)
plot(soil_0_30cm_shp, add=TRUE)
terrain_features_3m_df <- as.data.frame(Mar2017_terrain_3m_cropped)
dim(terrain_features_3m_df)
#terrain_features_3m_df <- terrain_features_3m_df[!is.na(terrain_features_3m_df$maxNDVI_2017),] #get rid of NAs
#dim(terrain_features_3m_df)
terrain_features_3m_df <- terrain_features_3m_df[,colnames(terrain_features_3m_df) %in% c('curvature_mean', 'elevation', 'slope', 'annual_kwh.m2', 'NDVI_2017mean_1m')] #'slope', 'annual_kwh.m2'; 3 classes with these works pretty well:'curvature_mean', 'elevation', 'meanNDVI_2017'
dim(terrain_features_3m_df)
terrain_features_3m_df_norm <- as.data.frame(lapply(terrain_features_3m_df, function(x) {(x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE)}))
colnames(terrain_features_3m_df_norm) <- c('curvature_mean_norm', 'elevation_norm', 'slope_norm', 'annual_kwh.m2_norm', 'NDVI_2017mean_1m_norm')

#get weights from lm
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + elevation_norm + annual_kwh.m2_norm + slope_norm, data=soil_0_30cm_shp)) #0.50
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + elevation_norm + annual_kwh.m2_norm, data=soil_0_30cm_shp)) #0.45
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + annual_kwh.m2_norm, data=soil_0_30cm_shp)) #r2=0.44
vif(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + annual_kwh.m2_norm, data=soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm, data=soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + annual_kwh.m2_norm + slope_norm, data=soil_0_30cm_shp))

vif(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + elevation_norm, data=soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm, data=soil_0_30cm_shp)) #r2=0.41
vif(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm, data=soil_0_30cm_shp))
lapply(terrain_features_3m_df_norm, summary)
kmeans_cluster <- function(classes, vars) {
  km.out.norm <- kmeans(na.omit(terrain_features_3m_df_norm[,vars]), classes) #verified this omits all rows where any var is NA
  catch_clusters <- rep(NA, nrow(terrain_features_3m_df_norm))
  catch_clusters[!is.na(terrain_features_3m_df_norm$NDVI_2017mean_1m_norm)] <- km.out.norm$cluster
  #Mar2017_terrain_3m_cropped$climate_cluster <- catch_clusters
  raster_object <- raster(extent(Mar2017_terrain_3m_cropped), resolution=res(Mar2017_terrain_3m_cropped), crs=crs(Mar2017_terrain_3m_cropped))
  catch_clusters <- setValues(raster_object, catch_clusters)
  #writeRaster(catch_clusters, filename = file.path(FiguresDir, 'cluster_rasters', paste0('cluster', classes, '_', length(vars), 'vars.tif')))
  plot(catch_clusters)
  cluster_ID <- extract(catch_clusters, soil_0_30cm_shp)
  cluster_ID
}
soil_0_30cm_shp$SOC_0_30cm_5var.prediction <- lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + elevation_norm + annual_kwh.m2_norm + slope_norm, data=soil_0_30cm_shp)$fitted.values
summary(soil_0_30cm_shp$SOC_0_30cm_5var.prediction)
#quintiles
quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, c(0.2, 0.4, 0.6, 0.8))
soil_0_30cm_shp$class_SOC_0_30cm_5var.quintile <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.215, 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.580, 2, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.812, 3, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 4.086, 4, 5))))
#sectiles
quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, c(1/6, 2/6, 3/6, 4/6, 5/6))
soil_0_30cm_shp$class_SOC_0_30cm_5var.sectile <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.143, 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.487, 2, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.675, 3, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.880, 4, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 4.147, 5, 6)))))

#quartiles
soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.358, 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.675, 2, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 4.000, 3, 4))) #breaks are from actual data
table(soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile, summary)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile, sd)
class_SOC_0_30cm_5var.quartile_sd <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile, sd))
class_SOC_0_30cm_5var.quartile_n <- as.numeric(table(soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile))
class_SOC_0_30cm_5var.quartile_CI_error <- qnorm(0.975) * class_SOC_0_30cm_5var.quartile_sd / sqrt(class_SOC_0_30cm_5var.quartile_n)
class_SOC_0_30cm_5var.quartile_CI_error

#now do terciles
quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, 0.333) #3.485
quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, 0.666) #3.879
soil_0_30cm_shp$class_SOC_0_30cm_5var.terciles <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.485, 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.879, 2, 3))

#now do halves
soil_0_30cm_shp$class_SOC_0_30cm_5var.lower.upper <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < 3.675, 1, 2)

soil_0_30cm_shp$class_2 <- kmeans_cluster(2, c('NDVI_2017mean_1m_norm', 'curvature_mean_norm'))
table(soil_0_30cm_shp$class_2)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2, summary)
class_2_sd <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2, sd))
class_2_n <- as.numeric(table(soil_0_30cm_shp$class_2))
class_2_CI_error <- qnorm(0.975) * class_2_sd / sqrt(class_2_n)
class_2_CI_error

soil_0_30cm_shp$class_3 <- kmeans_cluster(3, c('NDVI_2017mean_1m_norm', 'curvature_mean_norm'))
round(30*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
table(soil_0_30cm_shp$class_3)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_3, summary)
class_3_means <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_3, mean))
class_3_means
class_3_sd <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_3, sd))
class_3_n <- as.numeric(table(soil_0_30cm_shp$class_3))
class_3_CI_error <- qnorm(0.975) * class_3_sd / sqrt(class_3_n)
class_3_CI_error
class_3_upperCI <- class_3_means + qnorm(0.975) * class_3_sd / sqrt(class_3_n)
class_3_lowerCI <- class_3_means - qnorm(0.975) * class_3_sd / sqrt(class_3_n)
class_3_upperCI
class_3_lowerCI
  
soil_0_30cm_shp$class_4 <- kmeans_cluster(4, c('NDVI_2017mean_1m_norm', 'curvature_mean_norm'))
table(soil_0_30cm_shp$class_4)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4, summary) #2 classes don't have different means
tapply(soil_0_30cm_shp$NDVI_2017mean_1m, soil_0_30cm_shp$class_4, summary)
class_4_sd <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4, sd))
class_4_n <- as.numeric(table(soil_0_30cm_shp$class_4))
class_4_CI_error <- qnorm(0.975) * class_4_sd / sqrt(class_4_n)
class_4_CI_error

#most optimum number of classes
test2 <- pamk(na.omit(terrain_features_3m_df_norm[,c('curvature_mean_norm', 'NDVI_2017mean_1m_norm')]), krange=2:5, critout = TRUE) #3 clusters most efficient
test3 <- pamk(na.omit(terrain_features_3m_df_norm[,c('curvature_mean_norm', 'elevation_norm', 'slope_norm', 'annual_kwh.m2_norm', 'NDVI_2017mean_1m_norm')]), krange=2:5, critout = TRUE)

#NDVI classification
soil_0_30cm_shp$class_2_NDVI <- kmeans_cluster(2, 'NDVI_2017mean_1m_norm')
soil_0_30cm_shp$class_3_NDVI <- kmeans_cluster(3, 'NDVI_2017mean_1m_norm')

#5-var classification
soil_0_30cm_shp$class_2_5var <- kmeans_cluster(2, c('curvature_mean_norm', 'elevation_norm', 'slope_norm', 'annual_kwh.m2_norm', 'NDVI_2017mean_1m_norm'))
table(soil_0_30cm_shp$class_2_5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2_5var, summary)
class_2_5var_sd <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2_5var, sd))
class_2_5var_n <- as.numeric(table(soil_0_30cm_shp$class_2_5var))
class_2_5var_CI_error <- qnorm(0.975) * class_2_5var_sd / sqrt(class_2_5var_n)
class_2_5var_CI_error

soil_0_30cm_shp$class_3_5var <- kmeans_cluster(3, c('curvature_mean_norm', 'elevation_norm', 'slope_norm', 'annual_kwh.m2_norm', 'NDVI_2017mean_1m_norm'))
table(soil_0_30cm_shp$class_3_5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_3_5var, summary)

soil_0_30cm_shp$class_4_5var <- kmeans_cluster(4, c('curvature_mean_norm', 'elevation_norm', 'slope_norm', 'annual_kwh.m2_norm', 'NDVI_2017mean_1m_norm'))
table(soil_0_30cm_shp$class_4_5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4_5var, summary)
round(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4_5var, mean), 2)
round(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4_5var, sd), 2)
class_4_5var_sd <- as.numeric(tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4_5var, sd))
class_4_5var_n <- as.numeric(table(soil_0_30cm_shp$class_4_5var))
class_4_5var_CI_error <- qnorm(0.975) * class_4_5var_sd / sqrt(class_4_5var_n)
class_4_5var_CI_error


#try using principal components
test <- prcomp(~ curvature_mean_norm + elevation_norm + slope_norm + annual_kwh.m2_norm + NDVI_2017mean_1m_norm, data=terrain_features_3m_df_norm)
dim(test$x)
summary(test)
plot(test)
biplot(test)

PC_results_5var <- princomp(~ curvature_mean_norm + elevation_norm + slope_norm + annual_kwh.m2_norm + NDVI_2017mean_1m_norm, data=terrain_features_3m_df_norm)
dim(PC_results_5var$scores)
dim(terrain_features_3m_df_norm[!is.na(terrain_features_3m_df_norm$NDVI_2017mean_1m_norm),])
biplot(PC_results_5var)
test <- pamk(PC_results_5var$scores, krange = 2:5, critout = TRUE) #optimum cluster was 3


kmeans_PCcluster <- function(PC_input, classes) {
  km.out.norm <- kmeans(PC_input$scores, classes) #PC input has no NAs
  catch_clusters <- rep(NA, nrow(terrain_features_3m_df_norm))
  catch_clusters[!is.na(terrain_features_3m_df_norm$NDVI_2017mean_1m_norm)] <- km.out.norm$cluster
  #Mar2017_terrain_3m_cropped$climate_cluster <- catch_clusters
  raster_object <- raster(extent(Mar2017_terrain_3m_cropped), resolution=res(Mar2017_terrain_3m_cropped))
  catch_clusters <- setValues(raster_object, catch_clusters)
  plot(catch_clusters)
  cluster_ID <- extract(catch_clusters, soil_0_30cm_shp)
  cluster_ID
}
#verify which is most optimum number of clusters relative to PCs


#5-var classification
soil_0_30cm_shp$class_2_PC5var <- kmeans_PCcluster(PC_results_5var, 2)
table(soil_0_30cm_shp$class_2_PC5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2_PC5var, summary)

soil_0_30cm_shp$class_3_PC5var <- kmeans_PCcluster(PC_results_5var, 3)
table(soil_0_30cm_shp$class_3_PC5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_3_PC5var, summary)

soil_0_30cm_shp$class_4_PC5var <- kmeans_PCcluster(PC_results_5var, 4)
table(soil_0_30cm_shp$class_4_PC5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_4_PC5var, summary)

soil_0_30cm_shp$class_5_PC5var <- kmeans_PCcluster(PC_results_5var, 5)
table(soil_0_30cm_shp$class_5_PC5var)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_5_PC5var, summary)

#n <- 2
#class_no <- 3
calc_strat_mean <- function(n, class_no) {
  class_proportions <- tabulate(soil_0_30cm_shp[[paste0('class_', class_no)]]) / nrow(soil_0_30cm_shp)
  results <- numeric(length=class_no)
  if(n%%class_no==0) {
    for (i in seq_along(class_proportions)) {
      results[i] <- mean(soil_0_30cm_shp$kgOrgC.m2[sample(which(soil_0_30cm_shp[[paste0('class_', class_no)]]==(1:class_no)[i]), n/class_no)]) * class_proportions[i]
    }
  }
  sum(results)
}
#final results

# strat2_results <- replicate(10000, calc_strat_mean(2, 2))
# strat3_results <- replicate(10000, calc_strat_mean(3, 3))
# strat4_results <- replicate(10000, calc_strat_mean(4, 2))
# strat4_results_v4 <- replicate(10000, calc_strat_mean(4, 4))
# strat6_results <- replicate(10000, calc_strat_mean(6, 3))
# strat6_results_v2 <- replicate(10000, calc_strat_mean(6, 2))
# strat8_results <- replicate(10000, calc_strat_mean(8, 2))
# strat8_results_v4 <- replicate(10000, calc_strat_mean(8, 4))
# strat9_results <- replicate(10000, calc_strat_mean(9, 3))
# strat10_results <- replicate(10000, calc_strat_mean(10, 2))
# 
# calculate_thresholds(strat2_results, 10000) #0.9098 0.5951 0.3140
# calculate_thresholds(strat3_results, 10000) #0.9614 0.6887 0.3835
# calculate_thresholds(strat4_results, 10000) #0.9769 0.7482 0.4233
# calculate_thresholds(strat4_results_v4, 10000) #0.9839 0.7597 0.4394
# calculate_thresholds(strat6_results, 10000) #0.9975 0.8621 0.5463
# calculate_thresholds(strat6_results_v2, 10000) #0.9955 0.8569 0.5254
# calculate_thresholds(strat8_results, 10000) #0.9992 0.9062 0.6029
# calculate_thresholds(strat8_results_v4, 10000) #0.9996 0.9127 0.6084
# calculate_thresholds(strat9_results, 10000) #0.9997 0.9404 0.6468
# calculate_thresholds(strat10_results, 10000) #0.9999 0.9457 0.6662

calc_strat_mean_v2 <- function(class_no, classes, sample_no) {
  class_proportions <- tabulate(soil_0_30cm_shp[[paste0('class_', class_no)]]) / nrow(soil_0_30cm_shp)
  results <- numeric(length=classes)
  for (i in seq_along(class_proportions)) {
      results[i] <- mean(soil_0_30cm_shp$kgOrgC.m2[sample(which(soil_0_30cm_shp[[paste0('class_', class_no)]]==(1:classes)[i]), sample_no[i])]) * class_proportions[i]
  }
  sum(results)
}


#2 var (normalized mean curv. and mean 2017 NDVI) 2 class kmeans
round(2*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat2_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(1, 1)))
calculate_thresholds(strat2_class2, 10000) #0.9000 0.5652 0.2898

round(3*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat3_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(2, 1)))
calculate_thresholds(strat3_class2, 10000) #0.9611 0.6965 0.3915

round(4*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat4_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(3, 1)))
calculate_thresholds(strat4_class2, 10000) #0.9743 0.7629 0.4414

round(5*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat5_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(3, 2)))
calculate_thresholds(strat5_class2, 10000) #0.9921 0.8191 0.5017

round(6*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat6_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(4, 2)))
calculate_thresholds(strat6_class2, 10000) #0.9958 0.8622 0.5357

round(7*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat7_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(5, 2)))
calculate_thresholds(strat7_class2, 10000) #0.9984 0.8825 0.5721

round(8*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat8_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(5, 3)))
calculate_thresholds(strat8_class2, 10000) #0.9991 0.9113 0.6060

round(9*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat9_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(6, 3)))
calculate_thresholds(strat9_class2, 10000) #0.9998 0.9307 0.6321

round(10*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat10_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(7, 3)))
calculate_thresholds(strat10_class2, 10000) # 0.9999 0.9454 0.6610

round(15*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat15_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(10, 5)))
calculate_thresholds(strat15_class2, 10000) #1.0000 0.9883 0.7788

round(20*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat20_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(14, 6)))
calculate_thresholds(strat20_class2, 10000) # 1.0000 0.9959 0.8497

round(25*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat25_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(17, 8)))
calculate_thresholds(strat25_class2, 10000) #1.0000 0.9995 0.9035

round(30*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
strat30_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(21, 9)))
calculate_thresholds(strat30_class2, 10000) #1.0000 0.9995 0.9035

###2 var (normalized mean curv. and mean 2017 NDVI) 3 class kmeans
round(3*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat3_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(1, 1, 1)))
calculate_thresholds(strat3_class3, 10000) #0.9523 0.6599 0.3648

round(4*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat4_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(1, 2, 1)))
calculate_thresholds(strat4_class3, 10000) #0.9735 0.7631 0.4479

round(5*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat5_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(2, 2, 1)))
calculate_thresholds(strat5_class3, 10000) #0.9944 0.8254 0.5159

round(6*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat6_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(2, 3, 1)))
calculate_thresholds(strat6_class3, 10000) #0.9954 0.8668 0.5489

round(7*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat7_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(3, 3, 1)))
calculate_thresholds(strat7_class3, 10000) #0.9991 0.9007 0.5860

round(8*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat8_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(3, 4, 1)))
calculate_thresholds(strat8_class3, 10000) #0.9992 0.9243 0.6305

round(9*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat9_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(4, 4, 1)))
calculate_thresholds(strat9_class3, 10000) #0.9999 0.9392 0.6510

round(10*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat10_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(4, 5, 1)))
calculate_thresholds(strat10_class3, 10000) #0.9999 0.9519 0.6784

round(12*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat12_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(5, 6, 1)))
calculate_thresholds(strat12_class3, 10000) #0.9999 0.9519 0.6784

round(15*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat15_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(6, 7, 2)))
calculate_thresholds(strat15_class3, 10000) #1.0000 0.9897 0.7916

round(18*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat18_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(7, 9, 2)))
calculate_thresholds(strat18_class3, 10000)

round(20*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat20_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(8, 10, 2)))
calculate_thresholds(strat20_class3, 10000) #1.0000 0.9972 0.8629

round(21*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat21_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(9, 10, 2)))
calculate_thresholds(strat21_class3, 10000)

round(24*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat24_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(10, 11, 3)))
calculate_thresholds(strat24_class3, 10000) #1.0000 0.9996 0.9067

round(25*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat25_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(10, 12, 3)))
calculate_thresholds(strat25_class3, 10000) #1.0000 0.9990 0.9114

round(27*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat27_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(11, 13, 3)))
calculate_thresholds(strat27_class3, 10000) #1.0000 0.9996 0.9067

round(30*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat30_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(13, 14, 3)))
calculate_thresholds(strat30_class3, 10000) #1.0000 0.9999 0.9472

###2 var (normalized mean curv. and mean 2017 NDVI) 4 class kmeans
round(4*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat4_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(1, 1, 1, 1)))
calculate_thresholds(strat4_class4, 10000) #0.9743 0.7471 0.4228

round(5*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat5_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(2, 1, 1, 1)))
calculate_thresholds(strat5_class4, 10000) #0.9832 0.8156 0.4914

round(6*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat6_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(2, 1, 1, 2)))
calculate_thresholds(strat6_class4, 10000) #0.9965 0.8640 0.5514

round(7*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat7_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(3, 1, 1, 2)))
calculate_thresholds(strat7_class4, 10000) #0.9976 0.8954 0.5773

round(8*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat8_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(3, 1, 2, 2)))
calculate_thresholds(strat8_class4, 10000) #0.9985 0.9151 0.6186

round(9*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat9_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(3, 1, 2, 3)))
calculate_thresholds(strat9_class4, 10000) #0.9999 0.9399 0.6496

round(10*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat10_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(4, 1, 2, 3)))
calculate_thresholds(strat10_class4, 10000) #1.0000 0.9482 0.6748

round(15*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat15_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(6, 1, 3, 5)))
calculate_thresholds(strat15_class4, 10000) #1.0000 0.9864 0.7798

round(20*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat20_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(8, 1, 4, 7)))
calculate_thresholds(strat20_class4, 10000) #1.0000 0.9970 0.8491

round(25*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat25_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(10, 2, 5, 8)))
calculate_thresholds(strat25_class4, 10000) #1.0000 0.9994 0.9061

round(30*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
strat30_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(12, 2, 6, 10)))
calculate_thresholds(strat30_class4, 10000) #1.0000 1.0000 0.9447

#2 class approach with 5 var
round(2*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat2_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(1, 1)))
calculate_thresholds(strat2_class2_5var, 10000) #0.8862 0.5784 0.3057

round(3*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat3_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(1, 2)))
calculate_thresholds(strat3_class2_5var, 10000) #0.9480 0.6677 0.3660

round(4*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat4_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(1, 3)))
calculate_thresholds(strat4_class2_5var, 10000) #0.9722 0.7340 0.4180

round(5*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat5_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(2, 3)))
calculate_thresholds(strat5_class2_5var, 10000) #0.9906 0.7978 0.4748

round(6*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat6_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(2, 4)))
calculate_thresholds(strat6_class2_5var, 10000) #0.9953 0.8418 0.5143

round(7*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat7_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(2, 5)))
calculate_thresholds(strat7_class2_5var, 10000) #0.9961 0.8695 0.5426

round(8*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat8_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(3, 5)))
calculate_thresholds(strat8_class2_5var, 10000) #0.9992 0.9023 0.5934

round(9*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat9_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(3, 6)))
calculate_thresholds(strat9_class2_5var, 10000) #0.9993 0.9160 0.6111

round(10*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat10_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(3, 7)))
calculate_thresholds(strat10_class2_5var, 10000) # 0.9996 0.9339 0.6320

round(15*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat15_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(5, 10)))
calculate_thresholds(strat15_class2_5var, 10000) #1.0000 0.9807 0.7521

round(20*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat20_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(7, 13)))
calculate_thresholds(strat20_class2_5var, 10000) # 1.0000 0.9944 0.8381

round(25*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat25_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(8, 17)))
calculate_thresholds(strat25_class2_5var, 10000) #1.0000 0.9990 0.8915

round(30*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
strat30_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(10, 20)))
calculate_thresholds(strat30_class2_5var, 10000)

#3 class approach with 5 var
round(3*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat3_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 1, 1)))
calculate_thresholds(strat3_class3_5var, 10000) #0.9518 0.6656 0.3650

round(4*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat4_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 2, 1)))
calculate_thresholds(strat4_class3_5var, 10000) #0.9722 0.7340 0.4180

round(5*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat5_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 2, 2)))
calculate_thresholds(strat5_class3_5var, 10000) #0.9906 0.7978 0.4748

round(6*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat6_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 3, 2)))
calculate_thresholds(strat6_class3_5var, 10000) #0.9953 0.8418 0.5143

round(7*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat7_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(2, 3, 2)))
calculate_thresholds(strat7_class3_5var, 10000) #0.9961 0.8695 0.5426

round(8*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat8_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(2, 4, 2)))
calculate_thresholds(strat8_class3_5var, 10000) #0.9992 0.9023 0.5934

round(9*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat9_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(2, 4, 3)))
calculate_thresholds(strat9_class3_5var, 10000) #0.9993 0.9160 0.6111

round(10*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat10_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(3, 4, 3)))
calculate_thresholds(strat10_class3_5var, 10000) # 0.9996 0.9339 0.6320

round(15*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat15_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(4, 6, 5)))
calculate_thresholds(strat15_class3_5var, 10000) #1.0000 0.9807 0.7521

round(20*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat20_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(5, 9, 6)))
calculate_thresholds(strat20_class3_5var, 10000) # 1.0000 0.9944 0.8381

round(25*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat25_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(6, 11, 8)))
calculate_thresholds(strat25_class3_5var, 10000) #1.0000 0.9987 0.9053

round(30*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
strat30_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(8, 13, 9)))
calculate_thresholds(strat30_class3_5var, 10000) #1.0000 0.9997 0.9406

#4 class approach with 5 var
round(4*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat4_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(1, 1, 1, 1)))
calculate_thresholds(strat4_class4_5var, 10000) # 0.9883 0.7799 0.4559

round(5*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat5_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 1, 1, 1)))
calculate_thresholds(strat5_class4_5var, 10000) #0.9952 0.8134 0.4823

round(6*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat6_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 1, 2, 1)))
calculate_thresholds(strat6_class4_5var, 10000) #0.9981 0.8459 0.5182

round(7*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat7_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 2, 2, 1)))
calculate_thresholds(strat7_class4_5var, 10000) #0.9994 0.8863 0.5669

round(8*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat8_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 2, 2, 2)))
calculate_thresholds(strat8_class4_5var, 10000) #0.9997 0.9207 0.6241

round(9*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat9_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(3, 2, 2, 2)))
calculate_thresholds(strat9_class4_5var, 10000) #0.9999 0.9407 0.6484

round(10*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat10_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(3, 2, 3, 2)))
calculate_thresholds(strat10_class4_5var, 10000) #0.9999 0.9513 0.6667

round(15*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat15_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(4, 4, 4, 3)))
calculate_thresholds(strat15_class4_5var, 10000) #1.0000 0.9884 0.7933

round(20*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat20_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(6, 5, 5, 4)))
calculate_thresholds(strat20_class4_5var, 10000) #1.0000 0.9969 0.8624

round(25*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat25_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(7, 6, 7, 5)))
calculate_thresholds(strat25_class4_5var, 10000) #1.0000 0.9997 0.9121

round(30*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
strat30_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(8, 8, 8, 6)))
calculate_thresholds(strat30_class4_5var, 10000) #1.0000 1.0000 0.9462

#class_no <- '3_NDVI'
#sample_no_total <- 6
#function to automatically calc sampling numbers
calc_strat_mean_v3 <- function(class_no, classes, sample_no_total) {
  class_proportions <- tabulate(soil_0_30cm_shp[[paste0('class_', class_no)]]) / nrow(soil_0_30cm_shp)
  if(sample_no_total == classes) {
    sample_no <- rep(1, classes)
  } else {
      sample_no <- round(ifelse(class_proportions * sample_no_total < 1, 1, class_proportions * sample_no_total), 0)
  }
  if (sum(sample_no) != sum(sample_no_total)) {
    diff <- sum(sample_no) - sum(sample_no_total)
    x <- max(sample_no[sample_no > 1] -(class_proportions*sample_no_total)[sample_no > 1])
    sample_no[which(sample_no - class_proportions * sample_no_total == x)] <- sample_no[which(sample_no - class_proportions * sample_no_total == x)] - diff
  }
  #print(sample_no)
  results <- numeric(length=classes)
  for (i in seq_along(class_proportions)) {
    results[i] <- mean(soil_0_30cm_shp$kgOrgC.m2[sample(which(soil_0_30cm_shp[[paste0('class_', class_no)]]==(1:classes)[i]), sample_no[i])]) * class_proportions[i]
  }
  sum(results)
}
#5 var predictions classified into 5 groups based on the prediction quartiles (21 in 1 to 5, respectively)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_SOC_0_30cm_5var.quintile, summary)
table(soil_0_30cm_shp$class_SOC_0_30cm_5var.quintile)
strat5_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 5))
calculate_thresholds(strat5_classSOC30_5var.quintile, 10000) #0.9972 0.8894 0.5871

strat10_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 10))
calculate_thresholds(strat10_classSOC30_5var.quintile, 10000)

strat15_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 15))
calculate_thresholds(strat15_classSOC30_5var.quintile, 10000)

strat20_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 20))
calculate_thresholds(strat20_classSOC30_5var.quintile, 10000) #20 achieves 5% accuracy #  1.0000 0.9993 0.9109

strat24_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 4, 24))
calculate_thresholds(strat24_classSOC30_5var.quintile, 10000)

strat28_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 4, 28))
calculate_thresholds(strat28_classSOC30_5var.quintile, 10000)

#5 var predictions classified into 4 groups based on the prediction quartiles (26, 26, 26, and 27 in 1 to 4, respectively)
soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile
strat2_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 2))
calculate_thresholds(strat2_classSOC30_5var.quartile, 10000)

strat4_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 4))
calculate_thresholds(strat4_classSOC30_5var.quartile, 10000)

strat6_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 6))
calculate_thresholds(strat6_classSOC30_5var.quartile, 10000)

strat8_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 8))
calculate_thresholds(strat8_classSOC30_5var.quartile, 10000)#8 achieves 95% probability of 10% accuracy

strat12_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 12))
calculate_thresholds(strat12_classSOC30_5var.quartile, 10000) 

strat16_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 16))
calculate_thresholds(strat16_classSOC30_5var.quartile, 10000)

strat20_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 20))
calculate_thresholds(strat20_classSOC30_5var.quartile, 10000) #20 achieves 5% accuracy # 1.0000 0.9995 0.9099

strat24_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 24))
calculate_thresholds(strat24_classSOC30_5var.quartile, 10000)

strat28_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 28))
calculate_thresholds(strat28_classSOC30_5var.quartile, 10000)

#strat20_classSOC30_5var.quartile_test <- replicate(10000, calc_strat_mean_v2('SOC_0_30cm_5var.quartile', 4, c(3, 7, 7, 3)))
#calculate_thresholds(strat20_classSOC30_5var.quartile_test, 10000)

#now with the terciles from the 5 var predictions
strat3_classSOC30_5var.terciles <- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 3))
calculate_thresholds(strat3_classSOC30_5var.terciles, 10000)

strat6_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 6))
calculate_thresholds(strat6_classSOC30_5var.terciles, 10000) #0.9993 0.9075 0.6027

strat9_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 9))
calculate_thresholds(strat9_classSOC30_5var.terciles, 10000)#8 achieves 95% probability of 10% accuracy #0.9998 0.9647 0.7025

strat12_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 12))
calculate_thresholds(strat12_classSOC30_5var.terciles, 10000) #1.0000 0.9854 0.7719

strat15_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 15))
calculate_thresholds(strat15_classSOC30_5var.terciles, 10000)

strat18_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 18))
calculate_thresholds(strat18_classSOC30_5var.terciles, 10000)

strat21_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 21))
calculate_thresholds(strat21_classSOC30_5var.terciles, 10000)

strat24_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 24))
calculate_thresholds(strat24_classSOC30_5var.terciles, 10000)

strat27_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 27))
calculate_thresholds(strat27_classSOC30_5var.terciles, 10000)

strat30_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 30))
calculate_thresholds(strat30_classSOC30_5var.terciles, 10000)

#lower.upper approach
strat2_classSOC30_5var.lower.upper <- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.lower.upper', 2, 2))
calculate_thresholds(strat2_classSOC30_5var.lower.upper, 10000) #

#2 class approach with NDVI
strat2_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 2))
calculate_thresholds(strat2_class2_NDVI, 10000)

strat3_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 3))
calculate_thresholds(strat3_class2_NDVI, 10000)

strat4_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 4))
calculate_thresholds(strat4_class2_NDVI, 10000) 

strat5_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 5))
calculate_thresholds(strat5_class2_NDVI, 10000) 

strat6_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 6))
calculate_thresholds(strat6_class2_NDVI, 10000) 

strat7_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 7))
calculate_thresholds(strat7_class2_NDVI, 10000) #0.9994 0.8863 0.5669

strat8_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 8))
calculate_thresholds(strat8_class2_NDVI, 10000) #0.9997 0.9207 0.6241

strat9_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 9))
calculate_thresholds(strat9_class2_NDVI, 10000) #0.9999 0.9407 0.6484

strat10_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 10))
calculate_thresholds(strat10_class2_NDVI, 10000) #0.9999 0.9513 0.6667

strat15_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 15))
calculate_thresholds(strat15_class2_NDVI, 10000) #1.0000 0.9884 0.7933

strat20_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 20))
calculate_thresholds(strat20_class2_NDVI, 10000) #1.0000 0.9969 0.8624

strat25_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 25))
calculate_thresholds(strat25_class2_NDVI, 10000) #1.0000 0.9982 0.8888

strat30_class2_NDVI <- replicate(10000, calc_strat_mean_v3('2_NDVI', 2, 30))
calculate_thresholds(strat30_class2_NDVI, 10000) #1.0000 0.9994 0.9216

#3 class NDVI approach
#2 class approach with NDVI
strat3_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 3))
calculate_thresholds(strat3_class3_NDVI, 10000)

strat4_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 4))
calculate_thresholds(strat4_class3_NDVI, 10000) 

strat5_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 5))
calculate_thresholds(strat5_class3_NDVI, 10000) 

strat6_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 6))
calculate_thresholds(strat6_class3_NDVI, 10000) 

strat7_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 7))
calculate_thresholds(strat7_class3_NDVI, 10000) 

strat8_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 8))
calculate_thresholds(strat8_class3_NDVI, 10000) 

strat9_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 9))
calculate_thresholds(strat9_class3_NDVI, 10000) 

strat10_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 10))
calculate_thresholds(strat10_class3_NDVI, 10000) 

strat15_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 15))
calculate_thresholds(strat15_class3_NDVI, 10000) 

strat20_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 20))
calculate_thresholds(strat20_class3_NDVI, 10000) 

strat25_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 25))
calculate_thresholds(strat25_class3_NDVI, 10000) 

strat30_class3_NDVI <- replicate(10000, calc_strat_mean_v3('3_NDVI', 3, 30))
calculate_thresholds(strat30_class3_NDVI, 10000) 


results_2class_2var <- data.frame(n=c(2:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat2_class2, strat3_class2, strat4_class2, strat5_class2, strat6_class2, strat7_class2, strat8_class2, strat9_class2, strat10_class2, strat15_class2, strat20_class2, strat25_class2, strat30_class2), calculate_thresholds,  iterations=10000)))
colnames(results_2class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_2class_2var

results_3class_2var <- data.frame(n=seq(3, 30, 3), do.call(rbind, lapply(list(strat3_class3, strat6_class3, strat9_class3, strat12_class3, strat15_class3, strat18_class3, strat21_class3, strat24_class3, strat27_class3, strat30_class3), calculate_thresholds,  iterations=10000)))
colnames(results_3class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_3class_2var

# results_3class_2var <- data.frame(n=c(3:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat3_class3, strat4_class3, strat5_class3, strat6_class3, strat7_class3, strat8_class3, strat9_class3, strat10_class3, strat15_class3, strat20_class3, strat25_class3, strat30_class3), calculate_thresholds,  iterations=10000)))
# colnames(results_3class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
# results_3class_2var

results_3class_MLR.5var <- data.frame(n=seq(from=3, to=30, by=3), do.call(rbind, lapply(list(strat3_classSOC30_5var.terciles, strat6_classSOC30_5var.terciles, strat9_classSOC30_5var.terciles, strat12_classSOC30_5var.quartile, strat15_classSOC30_5var.terciles, strat18_classSOC30_5var.terciles, strat21_classSOC30_5var.terciles, strat24_classSOC30_5var.terciles, strat27_classSOC30_5var.terciles, strat30_classSOC30_5var.terciles), calculate_thresholds,  iterations=10000)))
colnames(results_3class_MLR.5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_3class_MLR.5var

results_4class_MLR.5var <- data.frame(n=c(4, 8, 12, 16, 20, 24, 28), do.call(rbind, lapply(list(strat4_classSOC30_5var.quartile, strat8_classSOC30_5var.quartile, strat12_classSOC30_5var.quartile, strat16_classSOC30_5var.quartile, strat20_classSOC30_5var.quartile, strat24_classSOC30_5var.quartile, strat28_classSOC30_5var.quartile), calculate_thresholds,  iterations=10000)))
colnames(results_4class_MLR.5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_4class_MLR.5var

results_4class_2var <- data.frame(n=c(4:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat4_class4, strat5_class4, strat6_class4, strat7_class4, strat8_class4, strat9_class4, strat10_class4, strat15_class4, strat20_class4, strat25_class4, strat30_class4), calculate_thresholds,  iterations=10000)))
colnames(results_4class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_4class_2var

results_2class_5var <- data.frame(n=c(2:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat2_class2_5var, strat3_class2_5var, strat4_class2_5var, strat5_class2_5var, strat6_class2_5var, strat7_class2_5var, strat8_class2_5var, strat9_class2_5var, strat10_class2_5var, strat15_class2_5var, strat20_class2_5var, strat25_class2_5var, strat30_class2_5var), calculate_thresholds,  iterations=10000)))
colnames(results_2class_5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_2class_5var

results_3class_5var <- data.frame(n=c(3:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat3_class3_5var, strat4_class3_5var, strat5_class3_5var, strat6_class3_5var, strat7_class3_5var, strat8_class3_5var, strat9_class3_5var, strat10_class3_5var, strat15_class3_5var, strat20_class3_5var, strat25_class3_5var, strat30_class3_5var), calculate_thresholds,  iterations=10000)))
colnames(results_3class_5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_3class_5var

results_4class_5var <- data.frame(n=c(4:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat4_class4_5var, strat5_class4_5var, strat6_class4_5var, strat7_class4_5var, strat8_class4_5var, strat9_class4_5var, strat10_class4_5var, strat15_class4_5var, strat20_class4_5var, strat25_class4_5var, strat30_class4_5var), calculate_thresholds,  iterations=10000)))
colnames(results_4class_5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_4class_5var

results_random_sampling <- data.frame(n=seq(from=3, to=30, by=3), do.call(rbind, lapply(list(sample_3, sample_6, sample_9, sample_12, sample_15, sample_18, sample_21, sample_24, sample_27, sample_30), calculate_thresholds, iterations=10000)))
colnames(results_random_sampling)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_random_sampling

#results_random_sampling <- data.frame(n=c(1:10, 15, 20, 25, 30), do.call(rbind, lapply(list(sample_1, sample_2, sample_3, sample_4, sample_5, sample_6, sample_7, sample_8, sample_9, sample_10, sample_15, sample_20, sample_25, sample_30), calculate_thresholds, iterations=10000)))
#colnames(results_random_sampling)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
#results_random_sampling

#3 class unsupervised strartification using 2 vars vs. 4 class MLR.5var vs.  random
tiff(file = file.path(FiguresDir, 'random_vs_strat_Fig8.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$`prob_5%`, ylim=c(0, 1), xlab='', ylab='', type='b', col='lightgrey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$`prob_10%`, type='b', pch=3, col='lightgrey', lty=2)
lines(results_random_sampling$n, results_random_sampling$`prob_20%`, type='b', pch=16, col='lightgrey', lty=3)
lines(results_3class_2var$n, results_3class_2var$`prob_5%`, type='b', col='black', lty=1, pch=1)
lines(results_3class_2var$n, results_3class_2var$`prob_10%`, type='b', col='black', lty=2, pch=3)
lines(results_3class_2var$n, results_3class_2var$`prob_20%`, type='b', col='black', lty=3, pch=16)
lines(results_3class_MLR.5var$n, results_3class_MLR.5var$`prob_5%`, type='b', col='brown', lty=1, pch=1)
lines(results_3class_MLR.5var$n, results_3class_MLR.5var$`prob_10%`, type='b', col='brown', lty=2, pch=3)
lines(results_3class_MLR.5var$n, results_3class_MLR.5var$`prob_20%`, type='b', col='brown', lty=3, pch=16)
#points(x=2, y=calculate_thresholds(strat2_class2, 10000)[1], col='darkgrey', pch = 8)
#points(x=2, y=calculate_thresholds(strat2_class2, 10000)[2], col='darkgrey', pch = 8)
#points(x=2, y=calculate_thresholds(strat2_class2, 10000)[3], col='darkgrey', pch = 8)
mtext(text='number of samples', side=1, line=2.75)
mtext(text=paste('probability to estimate mean SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified (MLR) ± 20% ', 'stratified (k-means) ± 20%', 'random ± 20%', 'stratified (MLR) ± 10%', 'stratified (k-means) ± 10%', 'random ± 10%', 'stratified (MLR) ± 5%', 'stratified (k-means) ± 5%', 'random ± 5%'), col=c('brown', 'black', 'lightgrey', 'brown', 'black', 'lightgrey', 'brown', 'black', 'lightgrey'), lty=c(3, 3, 3, 2, 2, 2, 1, 1, 1), pch=c(16, 16, 16, 3, 3, 3, 1, 1, 1), inset = 0.05)
abline(h=0.9, lty=4, col='black')
dev.off()

#3 class 2 var unsupervised strat
tiff(file = file.path(FiguresDir, 'random_vs_strat_3class_2var.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$`prob_5%`, ylim=c(0, 1), xlab='', ylab='', type='b', col='lightgrey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$`prob_10%`, type='b', pch=3, col='lightgrey', lty=2)
lines(results_random_sampling$n, results_random_sampling$`prob_20%`, type='b', pch=16, col='lightgrey', lty=3)
lines(results_3class_2var$n, results_3class_2var$`prob_5%`, type='b', col='black', lty=1, pch=1)
lines(results_3class_2var$n, results_3class_2var$`prob_10%`, type='b', col='black', lty=2, pch=3)
lines(results_3class_2var$n, results_3class_2var$`prob_20%`, type='b', col='black', lty=3, pch=16)
# lines(results_2class_2var$n, results_2class_2var$`prob_5%`, type='b', col='darkgrey', lty=1, pch=1)
# lines(results_2class_2var$n, results_2class_2var$`prob_10%`, type='b', col='darkgrey', lty=2, pch=3)
# lines(results_2class_2var$n, results_2class_2var$`prob_20%`, type='b', col='darkgrey', lty=3, pch=16)
points(x=2, y=calculate_thresholds(strat2_class2, 10000)[1], col='darkgrey', pch = 8)
points(x=2, y=calculate_thresholds(strat2_class2, 10000)[2], col='darkgrey', pch = 8)
points(x=2, y=calculate_thresholds(strat2_class2, 10000)[3], col='darkgrey', pch = 8)
mtext(text='number of samples', side=1, line=2.75)
mtext(text=paste('probability to estimate mean SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified ± X%, k=2', 'stratified ± 20%, k=3', 'random ± 20%', 'stratified ± 10%, k=3', 'random ± 10%', 'stratified ± 5%,  k=3', 'random ± 5%'), col=c('darkgrey', 'black', 'lightgrey', 'black', 'lightgrey', 'black', 'lightgrey'), lty=c(NA, 3, 3, 2, 2, 1, 1), pch=c(8, 16, 16, 3, 3, 1, 1), inset = 0.05)
abline(h=0.9, lty=4, col='black')
dev.off()

tiff(file = file.path(FiguresDir, 'random_vs_strat_2class_2var.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$`prob_5%`, ylim=c(0, 1), xlab='', ylab='', type='b', col='grey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$`prob_10%`, type='b', pch=3, col='grey', lty=2)
lines(results_random_sampling$n, results_random_sampling$`prob_20%`, type='b', pch=16, col='grey', lty=3)
lines(results_2class_2var$n, results_2class_2var$`prob_5%`, type='b', col='black', lty=1, pch=1)
lines(results_2class_2var$n, results_2class_2var$`prob_10%`, type='b', col='black', lty=2, pch=3)
lines(results_2class_2var$n, results_2class_2var$`prob_20%`, type='b', col='black', lty=3, pch=16)
mtext(text='number of samples', side=1, line=2.75)
mtext(text=paste('probability to estimate mean SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified ±20%', 'random ±20%', 'stratified ±10%', 'random ±10%', 'stratified ±5%', 'random ±5%'), col=c('black', 'grey', 'black', 'grey', 'black', 'grey'), lty=c(3, 3, 2, 2, 1, 1), pch=c(16, 16, 3, 3, 1, 1))
abline(h=0.9, lty=4, col='black')
dev.off()

tiff(file = file.path(FiguresDir, 'random_vs_strat_2class_5var.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$`prob_5%`, ylim=c(0, 1), xlab='', ylab='', type='b', col='grey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$`prob_10%`, type='b', pch=3, col='grey', lty=2)
lines(results_random_sampling$n, results_random_sampling$`prob_20%`, type='b', pch=16, col='grey', lty=3)
lines(results_2class_5var$n, results_2class_5var$`prob_5%`, type='b', col='black', lty=1, pch=1)
lines(results_2class_5var$n, results_2class_5var$`prob_10%`, type='b', col='black', lty=2, pch=3)
lines(results_2class_5var$n, results_2class_5var$`prob_20%`, type='b', col='black', lty=3, pch=16)
mtext(text='number of samples', side=1, line=2.75)
mtext(text=paste('probability to estimate mean SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified ±20%', 'random ±20%', 'stratified ±10%', 'random ±10%', 'stratified ±5%', 'random ±5%'), col=c('black', 'grey', 'black', 'grey', 'black', 'grey'), lty=c(3, 3, 2, 2, 1, 1), pch=c(16, 16, 3, 3, 1, 1))
abline(h=0.9, lty=4, col='black')
dev.off()

tiff(file = file.path(FiguresDir, 'random_vs_strat_3class_5var.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$`prob_5%`, ylim=c(0, 1), xlab='', ylab='', type='b', col='grey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$`prob_10%`, type='b', pch=3, col='grey', lty=2)
lines(results_random_sampling$n, results_random_sampling$`prob_20%`, type='b', pch=16, col='grey', lty=3)
lines(results_3class_5var$n, results_3class_5var$`prob_5%`, type='b', col='black', lty=1, pch=1)
lines(results_3class_5var$n, results_3class_5var$`prob_10%`, type='b', col='black', lty=2, pch=3)
lines(results_3class_5var$n, results_3class_5var$`prob_20%`, type='b', col='black', lty=3, pch=16)
mtext(text='number of samples', side=1, line=2.75)
mtext(text=paste('probability to estimate mean SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified ±20%', 'random ±20%', 'stratified ±10%', 'random ±10%', 'stratified ±5%', 'random ±5%'), col=c('black', 'grey', 'black', 'grey', 'black', 'grey'), lty=c(3, 3, 2, 2, 1, 1), pch=c(16, 16, 3, 3, 1, 1))
abline(h=0.9, lty=4, col='black')
dev.off()

tiff(file = file.path(FiguresDir, 'random_vs_strat_4class_5var.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$`prob_5%`, ylim=c(0, 1), xlab='', ylab='', type='b', col='grey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$`prob_10%`, type='b', pch=3, col='grey', lty=2)
lines(results_random_sampling$n, results_random_sampling$`prob_20%`, type='b', pch=16, col='grey', lty=3)
lines(results_4class_5var$n, results_4class_5var$`prob_5%`, type='b', col='black', lty=1, pch=1)
lines(results_4class_5var$n, results_4class_5var$`prob_10%`, type='b', col='black', lty=2, pch=3)
lines(results_4class_5var$n, results_4class_5var$`prob_20%`, type='b', col='black', lty=3, pch=16)
mtext(text='number of samples', side=1, line=2.75)
mtext(text=paste('probability to estimate mean SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified ± 20%, k=4', 'random ± 20%', 'stratified ± 10%, k=4', 'random ± 10%', 'stratified ± 5%, k=4', 'random ± 5%'), col=c('black', 'grey', 'black', 'grey', 'black', 'grey'), lty=c(3, 3, 2, 2, 1, 1), pch=c(16, 16, 3, 3, 1, 1))
abline(h=0.9, lty=4, col='black')
dev.off()

#test to see if weighting of each class by aerial proportion with even allocation of sampling resources
round(15*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
strat15_even_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(6, 7, 2)))
calculate_thresholds(strat15_even_class3, 10000) #1.0000 0.9897 0.7916


#now calculate distance from forage sampling points in 2017 to soil sampling points in 2018
# forage_data <- read.csv(file.path(forageDir, 'summaries', 'forage2017_2018_summary.csv'), stringsAsFactors=FALSE)
# waypoint_forage_sp <- shapefile(file.path(spatialforageDir, 'waypoint_forage2017.shp'))
# sensor_forage_sp <- shapefile(file.path(spatialforageDir, 'sensor_forage2017.shp'))
# sensor_forage_sp <- merge(sensor_forage_sp, forage_data[,c(1, 6:9)], by='location')
# names(waypoint_forage_sp)
# names(sensor_forage_sp)
# waypoint_forage_sp$clp011618 <- NA
# waypoint_forage_sp$clp021518 <- NA
# waypoint_forage_sp$clp032218 <- NA
# waypoint_forage_sp$clp041518 <- NA
# all_forage_sp <- rbind(sensor_forage_sp, waypoint_forage_sp)
# all_forage_sp <- all_forage_sp[-which(all_forage_sp$location=='B01'),] #outlier near road
# all_forage_sp$peak_2017 <- apply(as.data.frame(all_forage_sp)[,2:5], 1, max)
# all_forage_sp$peak_2018 <- apply(as.data.frame(all_forage_sp)[,6:9], 1, max)
# hist(all_forage_sp$peak_2017)
# hist(all_forage_sp$peak_2018)
# plot(all_forage_sp, cex=all_forage_sp$peak_2017/2000, col='red', pch=2, add=TRUE)
# shapefile(all_forage_sp, file.path(spatialforageDir, 'all_pts_2017_2018.shp'))
all_forage_sp <- shapefile(file.path(spatialforageDir, 'all_pts_2017_2018.shp'))
all_forage_sp$Mar2017growth <- all_forage_sp$clp031417 - all_forage_sp$clp021517
all_forage_sp$Apr2017growth <- all_forage_sp$clp041017 - all_forage_sp$clp031417
all_forage_sp$May2017growth <- all_forage_sp$clp050117 - all_forage_sp$clp041017
all_forage_sp$Mar2018growth <- all_forage_sp$clp032218 - all_forage_sp$clp021518
all_forage_sp$Apr2018growth <- all_forage_sp$clp041518 - all_forage_sp$clp032218
all_forage_sp$meanNDVI_2017 <- extract(meanNDVI_2017, all_forage_sp, buffer=1, fun=mean)
all_forage_sp$maxNDVI_2017 <- extract(maxNDVI_2017, all_forage_sp, buffer=1, fun=mean)
all_forage_sp$NDVI_02152017 <- extract(NDVI_2017$Camatta_02152017_NDVI, all_forage_sp, buffer=1, fun=mean)
all_forage_sp$NDVI_03172017 <- extract(NDVI_2017$Camatta_03172017_NDVI, all_forage_sp, buffer=1, fun=mean)
all_forage_sp$NDVI_04062017 <- extract(NDVI_2017$Camatta_04062017_NDVI, all_forage_sp, buffer=1, fun=mean)
all_forage_sp$curvature_mean <- extract(Mar2017_terrain_3m$curvature_mean, all_forage_sp)
#all_forage_sp$clp031417
summary(lm(clp021517 ~ NDVI_02152017, data=all_forage_sp))
summary(lm(clp031417 ~ NDVI_03172017, data=all_forage_sp)) #r2=0.09
summary(lm(clp041017 ~ NDVI_04062017, data=all_forage_sp)) #r2=0.16
summary(lm(peak_2017 ~ meanNDVI_2017, data=all_forage_sp)) #r2=0.31
summary(lm(peak_2017 ~ maxNDVI_2017, data=all_forage_sp)) #r2=0.19
summary(lm(peak_2017 ~ curvature_mean, data=all_forage_sp)) #r2=0.13
summary(lm(Apr2017growth ~ curvature_mean, data = all_forage_sp)) #r2=0.18
plot(all_forage_sp$curvature_mean, all_forage_sp$Apr2017growth)
summary(lm(May2017growth ~ curvature_mean, data=all_forage_sp))

#plot soil C as interpolated map
#see labs 14 and 15 from Quant Geo for tips
#also http://rspatial.org/analysis/rst/4-interpolation.html is more refined source of information

#first look at associations
lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), function(x) plot(x, soil_0_30cm_shp$kgOrgC.m2))
rank_test <- function(x, df, y, mtd) {
  test <- cor.test(x, df[[y]], method = mtd)
  result <- data.frame(col.1=test$p.value, col.2=test$estimate)
  colnames(result) <- c(paste0(y, '.p.val.', mtd), paste0(y, if(mtd=='pearson') {'.tau.'} else {'.rho.'}, mtd))
  result
}
soil_0_30cm_shp$kgOrgC.m2_0_10cm <- soil_0_10cm_shp$kgOrgC.m2
soil_0_30cm_shp$kgOrgC.m2_10_30cm <- soil_10_30cm_shp$kgOrgC.m2
orgC_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgOrgC.m2', mtd='spearman'))
orgC0_10cm_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgOrgC.m2_0_10cm', mtd='spearman'))
orgC10_30cm_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgOrgC.m2_10_30cm', mtd='spearman'))
P_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='gP.m2', mtd='spearman'))
IC_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgIC.m2', mtd='spearman'))
clay_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='clay_wtd', mtd='spearman'))
elevation_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='elevation', mtd='spearman'))
curvmean_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='curvature_mean', mtd='spearman'))
solrad_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='annual_kwh.m2', mtd='spearman'))
slope_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='slope', mtd='spearman'))
NDVI_2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NDVI_2017mean_1m', mtd='spearman'))
NDVI_2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NDVI_2018mean_1m', mtd='spearman'))
Red_GS_2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='Red_meanGS2017', mtd='spearman'))
NIR_GS_2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NIR_meanGS2017', mtd='spearman'))
Red_GS_2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='Red_meanGS2018', mtd='spearman'))
NIR_GS_2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'gP.m2', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NIR_meanGS2018', mtd='spearman'))
predictor_corrs <- cbind(orgC_corrs, orgC0_10cm_corrs, orgC10_30cm_corrs, clay_corrs, IC_corrs, P_corrs, elevation_corrs, curvmean_corrs, solrad_corrs, slope_corrs, NDVI_2017_corrs, NDVI_2018_corrs, Red_GS_2017_corrs, Red_GS_2018_corrs, NIR_GS_2017_corrs, NIR_GS_2018_corrs)
write.csv(predictor_corrs, file.path(CarbonDir, 'correlations', 'terrain_soil_corrs_0_30cm_orgCalldepths.csv'), row.names=TRUE)

#make plots of direct associations
tiff(file = file.path(FiguresDir, 'clay_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3.5, height = 3.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$clay_wtd, soil_0_30cm_shp$kgOrgC.m2, xlab=paste('0-30 cm clay (%)'), ylab=expression(paste('0-30 cm soil organic carbon (kg', ' ', m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #ylim = c(300, 1600), xlim=c(1400, 4700), col=soil_0_30cm_shp$energy_colors,
abline(lm(kgOrgC.m2 ~ clay_wtd, data = soil_0_30cm_shp), lty=2)
text(x=18, y=5.5, labels=expression(paste(r^2, '= 0.18')))
text(x=18, y=5.1,labels=paste('p-val < 0.001'))
dev.off()
summary(lm(kgOrgC.m2 ~ clay_wtd, data = soil_0_30cm_shp)) #r2=0.18
summary(lm(kgOrgC.m2 ~ CLAY, data = soil_0_10cm_shp)) #r2=0.05
summary(lm(kgOrgC.m2 ~ CLAY, data = soil_10_30cm_shp)) #r2=0.18

tiff(file = file.path(FiguresDir, 'elevation_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$elevation, soil_0_30cm_shp$kgOrgC.m2, xlab=paste('elevation (m)'), ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #ylim = c(300, 1600), xlim=c(1400, 4700), col=soil_0_30cm_shp$energy_colors
abline(lm(kgOrgC.m2 ~ elevation, data = soil_0_30cm_shp), lty=2)
text(x=475, y=2.5, labels=expression(paste(r^2, '= 0.11')), adj=c(0,0))
text(x=475, y=2.1,labels=paste('p-val < 0.001'), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ elevation, data = soil_0_30cm_shp)) #r2=0.11
summary(lm(kgOrgC.m2 ~ elevation, data = soil_0_10cm_shp)) #NS p.val=0.13
summary(lm(kgOrgC.m2 ~ elevation, data = soil_10_30cm_shp)) #r2=0.15

tiff(file = file.path(FiguresDir, 'mean_curv_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$curvature_mean, soil_0_30cm_shp$kgOrgC.m2, xlab=paste('mean curvature'), ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ curvature_mean, data = soil_0_30cm_shp), lty=2)
text(x=-1.5, y=2.5, labels=expression(paste(r^2, '= 0.24')), adj=c(0,0))
text(x=-1.5, y=2.1,labels=paste('p-val < 0.001'),  adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ curvature_mean, data = soil_0_30cm_shp)) #r2=0.24
summary(lm(kgOrgC.m2 ~ curvature_mean, data = soil_0_10cm_shp)) #r2=0.06
summary(lm(kgOrgC.m2 ~ curvature_mean, data = soil_10_30cm_shp)) #r2=0.32

tiff(file = file.path(FiguresDir, 'solrad_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$annual_kwh.m2, soil_0_30cm_shp$kgOrgC.m2, xlab=expression(paste('clear sky radiation (kWh', ' ', yr^-1, ')')), ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_0_30cm_shp), lty=2)
text(x=1150, y=2.5, labels=expression(paste(r^2, '= 0.03')), adj=c(0,0))
text(x=1150, y=2.1,labels=paste('p-val = 0.1'), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_0_30cm_shp)) #NS; pval=0.10
summary(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_0_10cm_shp)) #NS
summary(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_10_30cm_shp)) #r2=0.09

#elevation
tiff(file = file.path(FiguresDir, 'elevation_vs_clay_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3.5, height = 3.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$elevation, soil_0_30cm_shp$kgClay.m2, xlab=paste('elevation (m)'), ylab=expression(paste('0-30 cm clay (kg', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgClay.m2 ~ elevation, data = soil_0_30cm_shp), lty=2)
text(x=480, y=60, labels=expression(paste(r^2, '= 0.26')))
text(x=480, y=50,labels=paste('p-val < 0.001'))
dev.off()
summary(lm(kgClay.m2 ~ elevation, data = soil_0_30cm_shp)) #r2=0.26
summary(lm(kgClay.m2 ~ elevation, data = soil_0_10cm_shp)) #r2=0.11
summary(lm(kgClay.m2 ~ elevation, data = soil_10_30cm_shp)) #r2=0.28

#WMPD
tiff(file = file.path(FiguresDir, 'WMPDmm_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3.5, height = 3.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$WMPD_mm, soil_0_30cm_shp$kgOrgC.m2, xlab='', ylab=expression(paste('0-30 cm soil organic carbon (kg', ' ', m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
mtext(text=paste('Weighted-mean particle diameter (mm) 0-30 cm'), side=1, line=2.5, at=0.45)
abline(lm(kgOrgC.m2 ~ WMPD_mm, data = soil_0_30cm_shp), lty=2)
text(x=0.65, y=5.5, labels=expression(paste(r^2, '= 0.13')))
text(x=0.65, y=5.1,labels=paste('p-val < 0.001'))
dev.off()
summary(lm(kgOrgC.m2 ~ WMPD_mm, data = soil_0_30cm_shp)) #r2=0.13
summary(lm(kgOrgC.m2 ~ WMPD_mm, data = soil_0_10cm_shp)) #r2=0.04
summary(lm(kgOrgC.m2 ~ WMPD_mm, data = soil_10_30cm_shp)) #r2=0.12

#slope
tiff(file = file.path(FiguresDir, 'slope_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$slope, soil_0_30cm_shp$kgOrgC.m2, xlab='slope (%)', ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ slope, data = soil_0_30cm_shp), lty=2)
text(x=16, y=3.0, labels=expression(paste(r^2, '< 0.01')), adj=c(0,0))
text(x=16, y=2.6,labels=paste('p-val = 0.6'), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ slope, data = soil_0_30cm_shp)) #NS: p-val=0.61
summary(lm(kgOrgC.m2 ~ slope, data = soil_0_10cm_shp)) #r2=0.05
summary(lm(kgOrgC.m2 ~ slope, data = soil_10_30cm_shp)) #NS: p-val=0.16

#mean NDVI 2017
tiff(file = file.path(FiguresDir, 'NDVI2017_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$NDVI_2017mean_1m, soil_0_30cm_shp$kgOrgC.m2, xlab='mean 2017 NDVI', ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ NDVI_2017mean_1m, data = soil_0_30cm_shp), lty=2)
text(x=0.42, y=5.2, labels=expression(paste(r^2, '= 0.31')), adj=c(0,0))
text(x=0.42, y=4.8,labels=paste('p-val < 0.001'), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m, data = soil_0_30cm_shp)) #r2=0.31
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m, data = soil_0_10cm_shp))
summary(lm(kgOrgC.m2 ~ NDVI_2017max_1m, data = soil_10_30cm_shp))

#2018 mean NDVI
tiff(file = file.path(FiguresDir, 'NDVI2018_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$NDVI_2018mean_1m, soil_0_30cm_shp$kgOrgC.m2, xlab='mean 2018 NDVI', ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ NDVI_2018mean_1m, data = soil_0_30cm_shp), lty=2)
text(x=0.52, y=5.2, labels=expression(paste(r^2, '< 0.01')), adj=c(0,0))
text(x=0.52, y=4.8,labels=paste('p-val = 0.65'), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_0_30cm_shp)) #r^2<0.01 p-val=0.65
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_0_10cm_shp)) #r2 < 0.01 p-val=0.34
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_10_30cm_shp)) #r2 < 0.01 p-val=0.97

#first calculate NULL cross-validated model
#simplest approach to calculate null RMSE
null <- RMSE(soil_0_30cm_shp$kgOrgC.m2, mean(soil_0_30cm_shp$kgOrgC.m2))
null #0.698
#k-fold approach to calculate null RMSE
#k <- 2
#df_pts <- soil_0_30cm_shp
#varname <- 'kgOrgC.m2'
crossval_null <- function(df_pts, varname) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    rmse[k] <- RMSE(tst[[varname]], mean(trn[[varname]]))
    predictions[kf == k] <- mean(trn[[varname]])
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'kgOrgC.m2')
#mean(orgC_0_30_rmse_null$rmse.kfold) #0.6
orgC_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'kgOrgC.m2')
#mean(orgC_0_10_rmse_null$rmse.kfold) #0.472
orgC_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'kgOrgC.m2')
#mean(orgC_10_30_rmse_null$rmse.kfold) #0.452

clay_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'clay_wtd')
#mean(clay_0_30_rmse_null$rmse.kfold) #4.7% clay
clay_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'CLAY')
#mean(clay_0_10_rmse_null$rmse.kfold) #4.1% clay
clay_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'CLAY')
#mean(clay_10_30_rmse_null$rmse.kfold) #5.4% clay

# WMPD_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_null$rmse.kfold) #0.118 mm
# WMPD_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_null$rmse.kfold) #0.110 mm
# WMPD_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_null$rmse.kfold) #0.129 mm

IC_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'kgIC.m2')
#mean(IC_0_30_rmse_null$rmse.kfold) #1.15 kg IC m2
IC_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'kgIC.m2')
#mean(IC_0_10_rmse_null$rmse.kfold) #0.365 kg IC m2
IC_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'kgIC.m2')
#mean(IC_10_30_rmse_null$rmse.kfold) #0.904 kg IC m2

soilP_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'gP.m2')
#mean(soilP_0_30_rmse_null$rmse.kfold) #0.97 g soilP m2
soilP_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'gP.m2')
#mean(soilP_0_10_rmse_null$rmse.kfold) #0.67 g soilP m2
soilP_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'gP.m2')
#mean(soilP_10_30_rmse_null$rmse.kfold) #0.479 g P soilP m2

#map organic carbon 0-30 cm
#5 var model
lm_terrain5_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
kgOrgC.m2_terrain5_0_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_0_30cm)#, filename=file.path(FiguresDir, 'kgOrgC_m2_MLRbest5var_0_30cm_FINAL.tif'), overwrite=TRUE)
quantile(kgOrgC.m2_terrain5_0_30cm)
plot(kgOrgC.m2_terrain5_0_30cm)
plot(soil_0_30cm_shp, add=TRUE)
soil_0_30cm_shp$kgOrgC.m2_lm.terrain5 <- extract(kgOrgC.m2_terrain5_0_30cm, soil_0_30cm_shp)
summary(lm(kgOrgC.m2 ~ kgOrgC.m2_lm.terrain5, data = soil_0_30cm_shp))
soil_0_30cm_shp$kgOrgC.m2_5var_residuals <- lm_terrain5_0_30cm$residuals
#2 var model
lm_terrain2_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
kgOrgC.m2_terrain2_0_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain2_0_30cm)#, filename=file.path(FiguresDir, 'kgOrgC_m2_MLRbest2var_0_30cm_FINAL.tif'), overwrite=TRUE)
quantile(kgOrgC.m2_terrain2_0_30cm)
plot(kgOrgC.m2_terrain2_0_30cm)
soil_0_30cm_shp$kgOrgC.m2_lm.terrain2 <- extract(kgOrgC.m2_terrain2_0_30cm, soil_0_30cm_shp)
summary(lm(kgOrgC.m2 ~ kgOrgC.m2_lm.terrain2, data = soil_0_30cm_shp))
soil_0_30cm_shp$kgOrgC.m2_2var_residuals <- lm_terrain2_0_30cm$residuals

#check spatial correlation of residuals from 5 var model
orgC_krig <- gstat(formula=kgOrgC.m2 ~ 1, locations = soil_0_30cm_shp)
v <- variogram(orgC_krig)
fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
plot(variogramLine(fve, 150), type='l')
points(v[,2:3], pch=20, col='red')

residuals_2var_krig <- gstat(formula=kgOrgC.m2_2var_residuals ~ 1, locations=soil_0_30cm_shp)
v <- variogram(residuals_2var_krig)
fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
plot(variogramLine(fve, 150), type='l')
points(v[,2:3], pch=20, col='red')

residuals_5var_krig <- gstat(formula=kgOrgC.m2_5var_residuals ~ 1, locations=soil_0_30cm_shp)
v <- variogram(residuals_5var_krig)
fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
plot(variogramLine(fve, 150), type='l')
points(v[,2:3], pch=20, col='red')

#do the same for 0-10 cm layer
lm_terrain5_0_10cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_10cm_shp)
kgOrgC.m2_terrain5_0_10cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_0_10cm, filename=file.path(FiguresDir, 'kgOrgC_m2_MLRbest5var_0_10cm_FINAL.tif'))
plot(kgOrgC.m2_terrain5_0_10cm)
plot(soil_0_30cm_shp, add=TRUE)
soil_0_10cm_shp$kgOrgC.m2_lm.terrain5 <- extract(kgOrgC.m2_terrain5_0_10cm, soil_0_10cm_shp)
summary(lm(kgOrgC.m2 ~ kgOrgC.m2_lm.terrain5, data = soil_0_10cm_shp))
soil_0_10cm_shp$kgOrgC.m2_5var_residuals <- lm_terrain5_0_10cm$residuals

lm_terrain2_0_10cm <- lm(kgOrgC.m2 ~ slope + NDVI_2017mean_1m, data =  soil_0_10cm_shp)
kgOrgC.m2_terrain2_0_10cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain2_0_10cm, filename=file.path(FiguresDir, 'kgOrgC_m2_MLRbest2var_0_10cm_FINAL.tif'))
plot(kgOrgC.m2_terrain2_0_10cm)
plot(soil_0_30cm_shp, add=TRUE)
soil_0_10cm_shp$kgOrgC.m2_lm.terrain2 <- extract(kgOrgC.m2_terrain2_0_10cm, soil_0_10cm_shp)
summary(lm(kgOrgC.m2 ~ kgOrgC.m2_lm.terrain2, data = soil_0_10cm_shp))
soil_0_10cm_shp$kgOrgC.m2_2var_residuals <- lm_terrain2_0_10cm$residuals

#do the same for 10-30 cm layer
lm_terrain5_10_30cm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + annual_kwh.m2 + NIR_meanGS2017 + NDVI_2017mean_1m, data =  soil_10_30cm_shp)
#kgOrgC.m2_terrain5_10_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_10_30cm, filename=file.path(FiguresDir, 'kgOrgC_m2_MLRbest5var_10_30cm_FINAL.tif'))
#plot(kgOrgC.m2_terrain5_10_30cm)
#plot(soil_0_30cm_shp, add=TRUE)
#soil_10_30cm_shp$kgOrgC.m2_lm.terrain5 <- extract(kgOrgC.m2_terrain5_10_30cm, soil_10_30cm_shp)
#summary(lm(kgOrgC.m2 ~ kgOrgC.m2_lm.terrain5, data = soil_10_30cm_shp))
soil_10_30cm_shp$kgOrgC.m2_5var_residuals <- lm_terrain5_10_30cm$residuals

lm_terrain2_10_30cm <- lm(kgOrgC.m2 ~ curvature_mean + NIR_meanGS2018, data =  soil_10_30cm_shp)
#kgOrgC.m2_terrain2_10_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain2_10_30cm, filename=file.path(FiguresDir, 'kgOrgC_m2_MLRbest2var_10_30cm_FINAL.tif'))
#plot(kgOrgC.m2_terrain2_10_30cm)
#plot(soil_0_30cm_shp, add=TRUE)
#soil_10_30cm_shp$kgOrgC.m2_lm.terrain2 <- extract(kgOrgC.m2_terrain2_10_30cm, soil_10_30cm_shp)
#summary(lm(kgOrgC.m2 ~ kgOrgC.m2_lm.terrain2, data = soil_10_30cm_shp))
soil_10_30cm_shp$kgOrgC.m2_2var_residuals <- lm_terrain2_10_30cm$residuals

#check relationship between SOC at forage sampling points and peak forage in both years
all_forage_sp$kgOrgC.m2_lm.terrain5 <- extract(kgOrgC.m2_terrain5_0_30cm, all_forage_sp)
summary(lm(peak_2017 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp)) #r2=0.21
summary(lm(peak_2018 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp)) #r2=0.22
summary(lm(Mar2017growth ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))
summary(lm(Apr2017growth ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))
summary(lm(May2017growth ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))
summary(lm(Mar2018growth ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))
summary(lm(Apr2018growth ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))

#peak 2017 forage vs. modeled SOC
tiff(file = file.path(FiguresDir, 'peak2017_vs_orgC_0_30cm_MLRbest5var.tif', sep = ''), family = 'Times New Roman', width = 3.25, height = 3.3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(all_forage_sp$kgOrgC.m2_lm.terrain5, all_forage_sp$peak_2017, xlab=expression(paste('modeled SOC (kg ', ~m^-2, ')')), ylab=expression(paste('peak forage 2017 (kg ', ~ha^-1, ')')))
abline(lm(peak_2017 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp), lty=2)
text(3.6, 1400, label=expression(paste(r^2, '= 0.21')), adj=c(0,0))
text(3.6, 1100, label='p-val = 0.007', adj=c(0,0))
dev.off()
summary(lm(peak_2017 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))

#peak 2018 vs SOC
tiff(file = file.path(FiguresDir, 'peak2018_vs_orgC_0_30cm_MLRbest5var.tif', sep = ''), family = 'Times New Roman', width = 3.25, height = 3.3, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(all_forage_sp$kgOrgC.m2_lm.terrain5, all_forage_sp$peak_2018, xlab=expression(paste('modeled SOC (kg ', ~m^-2, ')')), ylab=expression(paste('peak forage 2018 (kg ', ~ha^-1, ')')))
abline(lm(peak_2018 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp), lty=2)
text(2.75, 1100, label=expression(paste(r^2, '= 0.22')), adj=c(0,0))
text(2.75, 1000, label='p-val = 0.06', adj=c(0.0))
dev.off()
summary(lm(peak_2018 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp))

hist(soil_0_30cm_shp$kgOrgC.m2)
hist(soil_0_30cm_shp$kgOrgC.m2_lm.terrain5)
quantile(kgOrgC.m2_terrain5_0_30cm, probs=c(0.25, 0.75)) 
#25%  75% 
#3.26 4.03
quantile(kgOrgC.m2_terrain4_0_30cm, probs=c(0.33, 0.66))
quantile(soil_0_30cm_shp$kgOrgC.m2, probs=c(0.33, 0.66))
#33%  66% 
#3.38 3.84

#map organic carbon 0-10 cm
lm_terrain4_0_10cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_0_10cm_shp)
kgOrgC.m2_terrain4_0_10cm <- predict(Mar2017_terrain_3m, lm_terrain4_0_10cm)
plot(kgOrgC.m2_terrain4_0_10cm)
writeRaster(kgOrgC.m2_terrain4_0_10cm, filename = file.path(FiguresDir, 'kgOrgC_lm.terrain4_0.10cm.tif'))
all_forage_sp$kgOrgC.m2_0_10_lm.terrain4 <- extract(kgOrgC.m2_terrain4_0_10cm, all_forage_sp)
summary(lm(all_forage_sp$peak_2017 ~ all_forage_sp$kgOrgC.m2_0_10_lm.terrain4))

#map organic carbon 10-30 cm
lm_terrain4_10_30cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_10_30cm_shp)
kgOrgC.m2_terrain4_10_30cm <- predict(Mar2017_terrain_3m, lm_terrain4_10_30cm)
plot(kgOrgC.m2_terrain4_10_30cm)
writeRaster(kgOrgC.m2_terrain4_10_30cm, filename = file.path(FiguresDir, 'kgOrgC_lm.terrain4_10.30cm.tif'))
all_forage_sp$kgOrgC.m2_10_30_lm.terrain4 <- extract(kgOrgC.m2_terrain4_10_30cm, all_forage_sp)
summary(lm(all_forage_sp$peak_2017 ~ all_forage_sp$kgOrgC.m2_10_30_lm.terrain4))

#map clay 0-30 cm (as WMPD)
clay_krig <- gstat(formula=clay_wtd ~ 1, locations=soil_0_30cm_shp)
v <- variogram(clay_krig)
fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
plot(variogramLine(fve, 150), type='l', main='% clay semivariance')
points(v[,2:3], pch=20, col='red')
gs <- gstat(formula=clay_wtd ~ 1, locations=soil_0_30cm_shp, model=fve) #used model created above
claywtd_0_30cm_ordkrig <- predict(gs, as(r, 'SpatialGrid'))
#test <- predict(r, gs)
claywtd_0_30cm_ordkrig <- brick(claywtd_0_30cm_ordkrig)
writeRaster(claywtd_0_30cm_ordkrig$var1.pred, filename = file.path(FiguresDir, 'claywtd_0_30cm_ordkrig.tif'), overwrite=TRUE)
plot(claywtd_0_30cm_ordkrig$var1.pred)
soil_0_30cm_shp$clay_wtd_ordkrig <- predict(gs, soil_0_30cm_shp)$var1.pred
plot(soil_0_30cm_shp$clay_wtd_ordkrig, soil_0_30cm_shp$clay_wtd)
Mar2017_terrain_3m_cropped$clay_wtd_ordkrig <- claywtd_0_30cm_ordkrig$var1.pred
plot(Mar2017_terrain_3m_cropped$clay_wtd_ordkrig)
lm_terrain3_clay_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + clay_wtd, data =  soil_0_30cm_shp)
Mar2017_terrain_3m_cropped$kgOrgC.m2_lmterrain3clay <- predict(Mar2017_terrain_3m_cropped, lm_terrain3_clay_0_30cm)
plot(Mar2017_terrain_3m_cropped$kgOrgC.m2_lmterrain3clay)
all_forage_sp$kgOrgC.m2_lmterrain3clay <- extract(Mar2017_terrain_3m_cropped$kgOrgC.m2_lmterrain3clay, all_forage_sp)
soil_0_30cm_shp$kgOrgC.m2_lm.terrain3clay <- extract(Mar2017_terrain_3m_cropped$kgOrgC.m2_lmterrain3clay, soil_0_30cm_shp)
plot(all_forage_sp$kgOrgC.m2_lmterrain3clay, all_forage_sp$peak_2017)
summary(lm(peak_2017 ~ kgOrgC.m2_lmterrain3clay, data = all_forage_sp))
summary(lm(peak_2018 ~ kgOrgC.m2_lmterrain3clay, data = all_forage_sp))
#writeRaster(Mar2017_terrain_3m_cropped$kgOrgC.m2_lmterrain3clay, filename = file.path(FiguresDir, 'kgOrgC_lm.terrain3clay_0.30cm.tif'))
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ soil_0_30cm_shp$kgOrgC.m2_lm.terrain3clay))

crossval_lm <- function(df_pts, varname, model='~ curvature_mean + slope + annual_kwh.m2 + elevation') {
  rmse <- rep(NA, length(unique(kf)))
  predictions <- rep(NA, 105)
  for (k in 1:length(unique(kf))) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    varname_lm <- lm(as.formula(paste(varname, model)), data = trn)
    #print(summary(varname_lm))
    varname_tst_pred <- predict.lm(varname_lm, tst)
    rmse[k] <- RMSE(tst[[varname]], varname_tst_pred)
    predictions[kf == k] <- varname_tst_pred
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}

orgC_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2', model = '~  curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m') #best model identified in model selection process in both 30-fold and 20-fold tests
mean(orgC_0_30_rmse_lm$rmse.kfold) #0.478
# orgC_0_30_rmse_lm4var <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2', model = '~  curvature_mean + annual_kwh.m2 + elevation + NDVI_2017mean_1m')
# mean(orgC_0_30_rmse_lm4var$rmse.kfold) #0.497

# orgC_0_30_rmse_lm5var <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2', model = '~  curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m') # + Red_May2017 + NDVI_2017max_1m + r^2=0.41; rmse=0.500 with NDVI max -or- r^2=0.44; rmse=0.49 with Red_May2017 (much worse with Red_Nov2016); -or- r2=0.44 rmse=0.481
# mean(orgC_0_30_rmse_lm5var$rmse.kfold) #0.478

# orgC_0_30_rmse_lm6var <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2', model = '~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_meanGS2018') #r^2=0.44 with NDVI and Red_May2017; NDVI_2017max_1m
# mean(orgC_0_30_rmse_lm6var$rmse.kfold) #0.481
# orgC_0_30_rmse_lm7var <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2', model = '~ NDVI_2017mean_1m + curvature_mean + annual_kwh.m2 + slope + elevation + NIR_May2017 + Red_May2017') #r^2=0.44
# mean(orgC_0_30_rmse_lm7var$rmse.kfold)  #0.487
# orgC_0_30_rmse_lm8var <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2', model = '~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Red_May2017 + NIR_May2017') #r^2=0.44
# mean(orgC_0_30_rmse_lm8var$rmse.kfold) #0.492
# plot(orgC_0_30_rmse_lm$oob.predictions, soil_0_30cm_shp$kgOrgC.m2)
orgC_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'kgOrgC.m2', model = '~  curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m') #best model identified in model selection process
mean(orgC_0_10_rmse_lm$rmse.kfold) #0.352 kg orgC m2
plot(orgC_0_10_rmse_lm$oob.predictions, soil_0_10cm_shp$kgOrgC.m2)
orgC_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'kgOrgC.m2', '~ elevation + annual_kwh.m2 + curvature_mean + NDVI_2017mean_1m + NIR_meanGS2017')
mean(orgC_10_30_rmse_lm$rmse.kfold) #0.264 kg orgC m2; r^2=0.53 oob
plot(orgC_10_30_rmse_lm$oob.predictions, soil_10_30cm_shp$kgOrgC.m2)

clay_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'clay_wtd', model='~ annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m + NIR_meanGS2018')
mean(clay_0_30_rmse_lm$rmse.kfold) #3.3% clay; r^2=0.42
clay_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'CLAY', model="~ elevation + annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m")
mean(clay_0_10_rmse_lm$rmse.kfold) #3.2% clay; r^2=0.36
clay_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'CLAY', model='~ elevation + annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m')
mean(clay_10_30_rmse_lm$rmse.kfold) #3.9% clay; r^2=0.38

# WMPD_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_lm$rmse.kfold) #0.094 mm; r^2=0.21
# WMPD_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_lm$rmse.kfold) #0.091 mm; r^2=0.16
# WMPD_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_lm$rmse.kfold) #0.102 mm; r^2=0.22

IC_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'kgIC.m2', model='~ elevation + annual_kwh.m2 + NDVI_2018mean_1m')
mean(IC_0_30_rmse_lm$rmse.kfold) #1.08 kg IC m2; r^2=0.10
IC_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'kgIC.m2', model='~ elevation + slope + NDVI_2018mean_1m')
mean(IC_0_10_rmse_lm$rmse.kfold) #0.348 kg IC m2; r^2=0.04 (p-val=0.03)
IC_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'kgIC.m2', model='~ annual_kwh.m2 + curvature_mean + NDVI_2017mean_1m')
mean(IC_10_30_rmse_lm$rmse.kfold) #0.849 kg IC m2; r^2=0.13

# soilP_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'gP.m2')
# mean(soilP_0_30_rmse_lm$rmse.kfold) #0.981 g soilP m2; r^2=0
# soilP_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'gP.m2')
# mean(soilP_0_10_rmse_lm$rmse.kfold) #0.72 g soilP m2; r^2=0.01
# soilP_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'gP.m2')
# mean(soilP_10_30_rmse_lm$rmse.kfold) #0.486 g P soilP m2; r^2=0

#model selection approach
#on 2/6/19, commented out write.csv calls to test if 30-fold CV affects results from original 20-fold CV run
model_selection_MLR <- function(df, varname, depth, varDir) {
  if(!dir.exists(file.path(modelResults, 'MLR_model_selection', varDir))) {
    dir.create(file.path(modelResults, 'MLR_model_selection', varDir))
  }
  models_to_test <- expand.grid(elevation=c(TRUE, FALSE), slope=c(TRUE, FALSE), annual_kwh.m2=c(TRUE, FALSE), curvature_mean=c(TRUE, FALSE), NDVI_2017mean_1m=c(TRUE, FALSE), NDVI_2018mean_1m=c(TRUE, FALSE), NIR_meanGS2017=c(TRUE, FALSE), NIR_meanGS2018=c(TRUE, FALSE))
  models_to_test <- models_to_test[1:(nrow(models_to_test)-1), ]
  model_selection_results <- list(length=nrow(models_to_test))
  #print(length(model_selection_results[1]))
  for (i in 1:nrow(models_to_test)) {
    print(i)
    model_to_test <- paste(colnames(models_to_test)[unlist(models_to_test[i,])], collapse = ' + ')
    model_selection_results[[i]] <- crossval_lm(df, varname, model = paste('~ ', model_to_test))
    print(mean(model_selection_results[[i]]$rmse.kfold))
  }
  mean_RMSEs <- unlist(lapply(model_selection_results, function(x) mean(x$rmse.kfold)))
  mean_RMSEs_all <- do.call(cbind, lapply(model_selection_results, function(x) x$rmse.kfold))
  oob_predictions_all <- do.call(cbind, lapply(model_selection_results, function(x) x$oob.predictions))
  #make summary
  df_model <- apply(models_to_test, 1, sum)
  meanRMSEs <- apply(mean_RMSEs_all, 2, mean)
  oob_r2s <- apply(oob_predictions_all, 2, function(x) summary(lm(df[[varname]] ~ x))$r.squared)
  summary <- data.frame(model=apply(models_to_test, 1, function(x) paste(colnames(models_to_test)[x], collapse = ' + ')), meanRMSE=meanRMSEs, oob_r.squared=oob_r2s, df_model=df_model)
  summary_by_model_df <- split(summary, summary$df_model)
  best_models <- unlist(lapply(summary_by_model_df, function(x) as.character(x$model[which.min(x$meanRMSE)])))
  best_rmses <- unlist(lapply(summary_by_model_df, function(x) min(x$meanRMSE)))
  best_oobs_r2 <- unlist(lapply(summary_by_model_df, function(x) x$oob_r.squared[which.min(x$meanRMSE)]))
  final_summary <- data.frame(model_df=1:8, model_name=best_models, meanRMSE_20foldCV=best_rmses, OOB_r2=best_oobs_r2)
  #write.csv(mean_RMSEs_all, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_CV_RMSEs.csv')), row.names = FALSE)
  #write.csv(oob_predictions_all, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_oob_pred.csv')), row.names = FALSE)
  #write.csv(models_to_test, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_model_test_grid.csv')), row.names = FALSE)
  #write.csv(final_summary, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_BEST_models.csv')), row.names = FALSE)
  list(RMSEs=mean_RMSEs_all, OOBs=oob_predictions_all, test_grid=models_to_test, best_models=final_summary)
}
#SOC content
orgC_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'kgOrgC.m2', '0_30cm', 'SOC')
orgC_0_30_MLR8var$best_models
orgC_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'kgOrgC.m2', '0_10cm', 'SOC')
orgC_0_10_MLR8var$best_models
orgC_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'kgOrgC.m2', '10_30cm', 'SOC')
orgC_10_30_MLR8var$best_models

#SOC percent
orgCper_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'orgC.percent', '0_30cm', 'SOC_percent')
orgCper_0_30_MLR8var$best_models
orgCper_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'orgC.percent', '0_10cm', 'SOC_percent')
orgCper_0_10_MLR8var$best_models
orgCper_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'orgC.percent', '10_30cm', 'SOC_percent')
orgCper_10_30_MLR8var$best_models

#SIC
SIC_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'kgIC.m2', '0_30cm', 'SIC')
SIC_0_30_MLR8var$best_models
SIC_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'kgIC.m2', '0_10cm', 'SIC')
SIC_0_10_MLR8var$best_models
SIC_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'kgIC.m2', '10_30cm', 'SIC')
SIC_10_30_MLR8var$best_models

#TC
TC_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'kgTC.m2', '0_30cm', 'TC')
TC_0_30_MLR8var$best_models
TC_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'kgTC.m2', '0_10cm', 'TC')
TC_0_10_MLR8var$best_models
TC_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'kgTC.m2', '10_30cm', 'TC')
TC_10_30_MLR8var$best_models

#clay
clay_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'clay_wtd', '0_30cm', 'clay')
clay_0_30_MLR8var$best_models
clay_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'CLAY', '0_10cm', 'clay')
clay_0_10_MLR8var$best_models
clay_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'CLAY', '10_30cm', 'clay')
clay_10_30_MLR8var$best_models

#Olsen P
P_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'gP.m2', '0_30cm', 'P')
P_0_30_MLR8var$best_models
P_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'gP.m2', '0_10cm', 'P')
P_0_10_MLR8var$best_models
P_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'gP.m2', '10_30cm', 'P')
P_10_30_MLR8var$best_models

#inverse distance weighted model for 0-30 cm organic carbon
library(gstat)
orgC_nn <- gstat(formula=kgOrgC.m2~1, locations=soil_0_30cm_shp, maxdist = 35, set=list(idp = 0)) #nearest neighbors approach:  idp=0
orgC_interpolated_nn <- interpolate(r, orgC_nn)
orgC_idw <- gstat(formula = kgOrgC.m2 ~ 1, locations = soil_0_30cm_shp, maxdist = 45, set=list(idp = 2.5))
orgC_interpolated_idw <- interpolate(r, orgC_idw)
summary(orgC_interpolated_nn)
summary(orgC_interpolated_idw)
plot(orgC_interpolated_nn)
plot(orgC_interpolated_idw)
soil_0_30cm_shp$orgC_est_idw <- extract(orgC_interpolated_idw, soil_0_30cm_shp) 
soil_0_30cm_shp$orgC_est_nn <- extract(orgC_interpolated_nn, soil_0_30cm_shp)
plot(soil_0_30cm_shp$orgC_est_idw, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$orgC_est_nn, soil_0_30cm_shp$kgOrgC.m2)
plot(all_forage_sp, cex=all_forage_sp$peak_2017/2000, col='red', pch=2, add=TRUE)
all_forage_sp$orgC_est_idw <- extract(orgC_interpolated_idw, all_forage_sp)
all_forage_sp$orgC_est_nn <- extract(orgC_interpolated_nn, all_forage_sp)
plot(all_forage_sp$orgC_est_idw, all_forage_sp$peak_2017)
plot(all_forage_sp$orgC_est_nn, all_forage_sp$peak_2017) #col=ifelse(all_forage_sp$location %in% as.character(1:16), 'red', 'black'))
summary(lm(peak_2017 ~ orgC_est_idw, data = all_forage_sp))
summary(lm(peak_2017 ~ orgC_est_nn, data = all_forage_sp))


#x <- c(3, 20, 0.2)
#test <- soil_0_30cm_shp[i,]
#train <- soil_0_30cm_shp[-i,]
f1 <- function(x, test, train) {
  #nmx <- x[1] if so desired
  mxdist <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula=kgOrgC.m2~1, locations=train, maxdist = mxdist, set=list(idp=idp)) #nmax=nmx if so desired
  p <- predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(test$kgOrgC.m2, p)
}
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
opt_results <- data.frame(maxdist=rep(NA, 20), idp=rep(NA, 20)) #if so desired nmax=rep(NA, 20)
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  opt <- optim(c(45, 1.8), f1, test=tst, train=trn)
  opt_results[k,] <- opt$par
}
opt_results
apply(opt_results, 2, mean)

#cross-validate idw model
predictions <- rep(NA, 105)
predictions[kf==k] <- p$var1.pred
print(summary(lm(df_pts[[varname]] ~ predictions)))
list(rmse.kfold=rmse, oob.predictions=predictions)
crossval_idw <- function(df_pts, varname, idp=1.7, maxdist=40) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    gs <- gstat(formula=as.formula(paste(varname, '~1')), locations=trn, maxdist = maxdist, set=list(idp=idp))
    p <- predict(gs, tst)
    rmse[k] <- RMSE(tst[[varname]], p$var1.pred)
    predictions[kf==k] <- p$var1.pred
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'kgOrgC.m2')
mean(orgC_0_30_rmse_idw$rmse.kfold) #0.5989389; r^2=0.19
orgC_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_idw$rmse.kfold) #0.3868148; r^2=0.08
orgC_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_idw$rmse.kfold) #0.3314461; r^2=0.29

clay_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_idw$rmse.kfold) #2.9% clay; r^2=0.63
clay_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_idw$rmse.kfold) #2.8% clay; r^2=0.51
clay_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_idw$rmse.kfold) #3.4% clay; r^2=0.58

# WMPD_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_idw$rmse.kfold) #0.055 mm; r^2=0.73
# WMPD_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_idw$rmse.kfold) #0.057 mm; r^2=0.65
# WMPD_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_idw$rmse.kfold) #0.064 mm; r^2=0.68

IC_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_idw$rmse.kfold) #1.05 kg IC m2; r^2=0.17
IC_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_idw$rmse.kfold) #0.345 kg IC m2; r^2=0.09
IC_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_idw$rmse.kfold) #0.835 kg IC m2; r^2=0.14

# soilP_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'gP.m2')
# mean(soilP_0_30_rmse_idw$rmse.kfold) #1.04 g soilP m2; r^2=0
# soilP_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'gP.m2')
# mean(soilP_0_10_rmse_idw$rmse.kfold) #0.768 g soilP m2; r^2=0
# soilP_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'gP.m2')
# mean(soilP_10_30_rmse_idw$rmse.kfold) #0.530 g P soilP m2; r^2=0.01

#cross validate nearest neighbors model
crossval_nn <- function(df_pts, varname, maxdist=25, idp=0) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    gs <- gstat(formula=as.formula(paste(varname, '~1')), locations=trn, maxdist = maxdist, set=list(idp=idp))
    p <- predict(gs, tst)
    rmse[k] <- RMSE(tst[[varname]], p$var1.pred)
    predictions[kf==k] <- p$var1.pred
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'kgOrgC.m2')
mean(orgC_0_30_rmse_nn$rmse.kfold) #0.6152075; r^2: 0.17
orgC_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_nn$rmse.kfold) #0.4063699; r^2: 0.08
orgC_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_nn$rmse.kfold) #0.3396059; r^2: 0.26

clay_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_nn$rmse.kfold) #2.80241% clay; r^2:0.64
clay_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_nn$rmse.kfold) #2.94755% clay; r^2:0.48
clay_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_nn$rmse.kfold) #3.312241% clay; r^2:0.60

# WMPD_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_nn$rmse.kfold) #0.05174892 mm; r^2:0.76
# WMPD_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_nn$rmse.kfold) #0.05679778 mm; r^2:0.65
# WMPD_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_nn$rmse.kfold) #0.06007388 mm; r^2:0.72

IC_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_nn$rmse.kfold) #0.981 kg IC m2; r^2:0.28
IC_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_nn$rmse.kfold) #0.325 kg IC m2; r^2:0.21
IC_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_nn$rmse.kfold) #0.802 kg IC m2; r^2:0.22

# soilP_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'gP.m2')
# mean(soilP_0_30_rmse_nn$rmse.kfold) # 1.13 g soilP m2; r^2:0
# soilP_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'gP.m2')
# mean(soilP_0_10_rmse_nn$rmse.kfold) #0.802 g soilP m2; r^2:0
# soilP_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'gP.m2')
# mean(soilP_10_30_rmse_nn$rmse.kfold) #0.549 g P soilP m2; r^2:0.02

#ordinary kriging
#borrowed code from http://rspatial.org/analysis/rst/4-interpolation.html
orgC_krig <- gstat(formula=kgOrgC.m2~1, locations=soil_0_30cm_shp)
v <- variogram(orgC_krig)
v
plot(v)
fve <- fit.variogram(v, vgm(psill=0.3, model="Sph", range=50, nugget=0.2))
fve
#model     psill    range
#1   Nug 0.1879249  0.00000
#2   Sph 0.3257305 63.60272
plot(variogramLine(fve, 100), type='l', ylim=c(0,0.7))
points(v[,2:3], pch=20, col='red')
krig_model <- gstat(formula=kgOrgC.m2~1, locations=soil_0_30cm_shp, model=fve)
# predicted values
orgC_krigged <- predict(krig_model, as(r, 'SpatialGrid'))
## [using ordinary kriging]
spplot(orgC_krigged)
orgC_krigged <- brick(orgC_krigged)
plot(orgC_krigged$var1.pred)
all_forage_sp$orgC_est_krig <- extract(orgC_krigged$var1.pred, all_forage_sp)
plot(all_forage_sp$orgC_est_krig, all_forage_sp$peak_2017)
abline(lm(peak_2017 ~ orgC_est_krig, data = all_forage_sp), lty=2)
summary(lm(peak_2017 ~ orgC_est_krig, data = all_forage_sp))
plot(all_forage_sp$orgC_est_krig, all_forage_sp$peak_2018)
abline(lm(peak_2018 ~ orgC_est_krig, data = all_forage_sp), lty=2)
summary(lm(peak_2018 ~ orgC_est_krig, data = all_forage_sp)) #r2=0.17; p-val=0.12
soil_0_30cm_shp$orgC_est_krig <- extract(orgC_krigged$var1.pred, soil_0_30cm_shp)
summary(lm(kgOrgC.m2 ~ orgC_est_krig, data = soil_0_30cm_shp)) #r2=0.87
plot(soil_0_30cm_shp$orgC_est_krig, soil_0_30cm_shp$kgOrgC.m2)

#cross-validate ordinary krigged
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
tapply(kf, kf, length)

crossval_ordkrig <- function(df_pts, varname, psill=0.3, model='Sph', range=50, nugget=0.2) {
  rmse <- rep(NA, length(unique(kf)))
  predictions <- rep(NA, 105)
  for (k in 1:length(unique(kf))) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    varname_krig <- gstat(formula=as.formula(paste(varname, '~1')), locations=trn)
    v <- variogram(varname_krig)
    #test <- autofitVariogram(formula=as.formula(paste(varname, '~1')), input_data = trn, verbose = FALSE, kappa=c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
    #v <- variogram(varname_krig, boundaries=test$exp_var$dist + 10)
    fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
    #fve <- autofitVariogram(formula=as.formula(paste(varname, '~1')), input_data = trn, verbose = FALSE, kappa=c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
    #fve <- fit.variogram(v, vgm(psill=fve$var_model$psill[2], model=as.character(fve$var_model$model)[2], range=fve$var_model$range[2], nugget = fve$var_model$psill[1], kappa = fve$var_model$kappa[2]))
    #v <- variogram(varname_krig)
    #fve <- fit.variogram(v, vgm(psill=psill, model=model, range=range, nugget=nugget))
    #print(fve)
    plot(variogramLine(fve, 150), type='l', main=paste0(k, '-fold plot'))
    points(v[,2:3], pch=20, col='red')
    gs <- gstat(formula=as.formula(paste(varname, '~1')), locations=trn, model=fve) #used model created above
    p <- predict(gs, tst)
    rmse[k] <- RMSE(tst[[varname]], p$var1.pred)
    predictions[kf==k] <- p$var1.pred
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'kgOrgC.m2')
mean(orgC_0_30_rmse_ordkrig$rmse.kfold) #0.611 exp model is better than Sph; r^2=0.15
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ orgC_0_30_rmse_ordkrig$oob.predictions))
orgC_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_ordkrig$rmse.kfold) #0.391; r^2=0.04
orgC_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_ordkrig$rmse.kfold) #0.339; r^2=0.26

clay_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_ordkrig$rmse.kfold) #2.79% clay; r^2=0.65
clay_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_ordkrig$rmse.kfold) #2.87% clay; r2=0.50
clay_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_ordkrig$rmse.kfold) #3.42% clay; r2=0.57

# WMPD_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_ordkrig$rmse.kfold) #0.05555203 mm; r^2:0.72
# WMPD_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_ordkrig$rmse.kfold) #0.05769833 mm; r^2:0.64
# WMPD_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_ordkrig$rmse.kfold) #0.0621748 mm; r^2:0.70

IC_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_ordkrig$rmse.kfold) #1.01 kg IC m2; r2=0.21
IC_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_ordkrig$rmse.kfold) #0.330 kg IC m2; r2=0.16
IC_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_ordkrig$rmse.kfold) #0.852 kg IC m2; r2=0.07

# soilP_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'gP.m2')
# mean(soilP_0_30_rmse_ordkrig$rmse.kfold) # 1.07 g soilP m2; r2=0
# soilP_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'gP.m2')
# mean(soilP_0_10_rmse_ordkrig$rmse.kfold) #0.811 g soilP m2; r2=0
# soilP_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'gP.m2')
# mean(soilP_10_30_rmse_ordkrig$rmse.kfold) # 0.535 g P soilP m2; r2=0

#regression krigging model 0-30 cm
library(relaimpo)
library(car)
library(gstat)
library(spdep)
autocorr_test_soil <- function(df_shp, varname, nsim) {
  set.seed(19801976)
  #then, make an inverse distance weighted matrix
  idw <- 1/pointDistance(df_shp, latlon=FALSE)  #equivalent to 1/as.matrix(dist(coordinates(forage_data_sp))), see GEO200CN lab 14
  diag(idw) <- 0 #set Inf back to zero
  idw_list <- mat2listw(idw)
  result <- moran.mc(df_shp[[varname]], idw_list, nsim = nsim)
  print(result)
  result
  results <- cbind(result$statistic, result$p.value)
  results <- as.data.frame(results)
  colnames(results) <- c('Moran I statistic', 'p_value')
  results$n_pts <- nrow(df_shp)
  results$varname <- varname
  results
}
names(soil_0_10cm_shp)
#Monte-Carlo simulation of Moran I
soil_0_30_autocorr <- do.call(rbind, lapply(c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'clay_wtd', 'sand_wtd', 'silt_wtd', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'kgOrgC.m2_2var_residuals', 'kgOrgC.m2_5var_residuals'), function(x) autocorr_test_soil(soil_0_30cm_shp, varname = x, nsim = 999)))
soil_0_30_autocorr
write.csv(soil_0_30_autocorr, file.path(modelResults, 'autocorrelation', 'soil_0_30cm_autocorrelation_2_12_19.csv'), row.names = FALSE)

names(soil_0_10cm_shp)
soil_0_10_autocorr <- do.call(rbind, lapply(c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'kgOrgC.m2_2var_residuals', 'kgOrgC.m2_5var_residuals'), function(x) autocorr_test_soil(soil_0_10cm_shp, varname = x, nsim = 999)))
soil_0_10_autocorr
write.csv(soil_0_10_autocorr, file.path(modelResults, 'autocorrelation', 'soil_0_10cm_autocorrelation.csv'), row.names = FALSE)

soil_10_30_autocorr <- do.call(rbind, lapply(c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'kgOrgC.m2_2var_residuals', 'kgOrgC.m2_5var_residuals'), function(x) autocorr_test_soil(soil_10_30cm_shp, varname = x, nsim = 999)))
soil_10_30_autocorr
write.csv(soil_10_30_autocorr, file.path(modelResults, 'autocorrelation', 'soil_10_30cm_autocorrelation.csv'), row.names = FALSE)

#regression kriging approach
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + coords.x1 + coords.x2, data=soil_0_30cm_shp))
orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + annual_kwh.m2 + NDVI_2017mean_1m, data=soil_0_30cm_shp)
summary(orgC_lm)
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + annual_kwh.m2 + WMPD_mm, data=soil_0_30cm_shp))
vif(orgC_lm)
orgC_terrain_pred <- predict(Mar2017_terrain_3m, orgC_lm, fun=predict)
plot(orgC_terrain_pred)
orgC_terrain_pred <- resample(orgC_terrain_pred, r, method='bilinear')
plot(orgC_terrain_pred)
soil_0_30cm_shp$org_lm_residuals <- orgC_lm$residuals
#check autocorrelation in residuals
soil_0_30cm_shp$orgC_lm_predictions <- orgC_lm$fitted.values
summary(lm(kgOrgC.m2 ~ orgC_lm_predictions, data = soil_0_30cm_shp))
plot(soil_0_30cm_shp$orgC_lm_predictions, soil_0_30cm_shp$kgOrgC.m2)
orgC_reg_krig <- gstat(formula=org_lm_residuals~1, locations =  soil_0_30cm_shp)
v <- variogram(orgC_reg_krig) #not specifying distance here
test <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Lin', 'Bes', 'Log', 'Ste')), fit.kappa = seq(0.3,30,0.05))
as.character(test$model)[2]
test

library(sp)
library(automap)
test <- autofitVariogram(formula=org_lm_residuals~1, input_data = soil_0_30cm_shp, verbose = FALSE)
test$exp_var
test$var_model
plot(test)
v <- variogram(orgC_reg_krig, boundaries=test$exp_var$dist + 10)
v
mean(v$dist[2:length(v$dist)] - v$dist[1:(length(v$dist)-1)])
show.vgms()

fve <- test
fve <- fit.variogram(v, vgm(psill=fve$var_model$psill[2], model=as.character(fve$var_model$model)[2], range=fve$var_model$range[2], nugget = fve$var_model$psill[1], kappa = fve$var_model$kappa[2]), fit.kappa = seq(0.3,30,0.05))
fve

#fve <- fit.variogram(v, vgm(psill=0.3, model="Exp", range=50, nugget = 0.2)) #works for width = 21
#fve <- fit.variogram(v, vgm(psill=2.3, model="Sph", range=2900, nugget = 0.2))
#finally, convergence!
#model      psill    range
#1   Nug 0.23304615   0.0000
#2   Sph 0.07379567 168.1028
plot(variogramLine(fve, 200), type='l', ylim=c(0,0.7))
points(v[,2:3], pch=20, col='red')
regkrig_model <- gstat(formula=org_lm_residuals ~ 1, locations = soil_0_30cm_shp, model=fve)
# predicted values
orgC_res_regkrig <- predict(regkrig_model, as(r, 'SpatialGrid'))
## [using ordinary kriging]
spplot(orgC_res_regkrig)
orgC_res_regkrig <- brick(orgC_res_regkrig)
plot(orgC_res_regkrig$var1.pred)
soil_0_30cm_shp$orgC_regkrig_correction <- extract(orgC_res_regkrig$var1.pred, soil_0_30cm_shp) 
plot(soil_0_30cm_shp$WMPD_mm, soil_0_30cm_shp$orgC_regkrig_correction)
summary(lm(soil_0_30cm_shp$orgC_regkrig_correction ~ soil_0_30cm_shp$WMPD_mm))
plot(soil_0_30cm_shp$orgC_lm_predictions, soil_0_30cm_shp$org_lm_residuals)
plot(orgC_terrain_pred)
orgC_regkrig_est <- orgC_res_regkrig$var1.pred + orgC_terrain_pred
plot(orgC_regkrig_est)
soil_0_30cm_shp$orgC_est_regkrig <- extract(orgC_regkrig_est, soil_0_30cm_shp)
plot(soil_0_30cm_shp$kgClay.m2, soil_0_30cm_shp$orgC_est_regkrig)
plot(soil_0_30cm_shp$orgC_est_regkrig, soil_0_30cm_shp$kgOrgC.m2)
summary(lm(kgOrgC.m2 ~ orgC_est_regkrig, data=soil_0_30cm_shp)) #r^2=0.5; r^2=0.99 when using other manually selected variogram object that results in more exact interpolation
all_forage_sp$orgC_est_regkrig <- extract(orgC_regkrig_est, all_forage_sp)
plot(all_forage_sp$orgC_est_regkrig, all_forage_sp$peak_2017)
abline(lm(peak_2017 ~ orgC_est_regkrig, data = all_forage_sp), lty=2)
summary(lm(peak_2017 ~ orgC_est_regkrig, data = all_forage_sp)) #r^2=0.23;p-val=0.006
plot(all_forage_sp$orgC_est_regkrig, all_forage_sp$peak_2018)
abline(lm(peak_2018 ~ orgC_est_regkrig, data = all_forage_sp), lty=2)
summary(lm(peak_2018 ~ orgC_est_regkrig, data = all_forage_sp)) #r^2=0.17;p-val=0.11

#alternative set-up is to uncomment lines labeled (a), (b), (c), and (d) and comment lines labeled (1) and (2)
crossval_regkrig <- function(df_pts, varname, model='~ curvature_mean + elevation + slope + annual_kwh.m2') {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    varname_lm <- lm(as.formula(paste(varname, model)), data=trn)
    trn$residuals <- varname_lm$residuals
    # terrain_pred <- predict(Mar2017_terrain_3m, varname_lm, fun=predict)
    # tst_pred <- extract(terrain_pred, tst)
    tst_pred <- predict.lm(varname_lm, tst)
    reg_krig <- gstat(formula=residuals ~ 1, locations = trn)
    test <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE) # (a)
    # v <- variogram(reg_krig) #(1)
    v <- variogram(reg_krig, boundaries=test$exp_var$dist + 10) # (b)
    fve <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE, kappa=c(0.05, seq(0.1,5,0.1), seq(5, 250, 5))) #(c)
    fve <- fit.variogram(v, vgm(psill=fve$var_model$psill[2], model=as.character(fve$var_model$model)[2], range=fve$var_model$range[2], nugget = fve$var_model$psill[1], kappa = fve$var_model$kappa[2])) #(d)
    #fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5))) #(2)
    #fve <- fit.variogram(v, vgm(psill=0.1, model="Sph", range=150, nugget = 0.2))
    regkrig_model <- gstat(formula=residuals ~ 1, locations = trn, model=fve)
    p_res_correction <- predict(regkrig_model, tst)
    p <- p_res_correction$var1.pred + tst_pred
    rmse[k] <- RMSE(tst[[varname]], p)
    predictions[kf == k] <- p
    print(fve)
    plot(variogramLine(fve, 150), type='l', main=paste0(k, '-fold plot'))
    points(v[,2:3], pch=20, col='red')
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}

orgC_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'kgOrgC.m2', model='~  curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m')
plot(orgC_0_30_rmse_regkrig$oob.predictions, soil_0_30cm_shp$kgOrgC.m2)
abline(0,1,lty=2)
mean(orgC_0_30_rmse_regkrig$rmse.kfold) #rmse: 0.477; r^2:0.44 (same as MLR)
orgC_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'kgOrgC.m2', model = '~ elevation + slope + annual_kwh.m2 + curvature_mean + NDVI_2017mean_1m')
mean(orgC_0_10_rmse_regkrig$rmse.kfold) #0.357; r^2=0.16
orgC_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'kgOrgC.m2', model='~ elevation + annual_kwh.m2 + curvature_mean + NDVI_2017mean_1m + NIR_meanGS2017')
mean(orgC_10_30_rmse_regkrig$rmse.kfold) #0.268; r^2=0.52
#all orgC regression krigging are the exact same as MLR results


clay_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'clay_wtd', model='~ annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m + NIR_meanGS2018')
mean(clay_0_30_rmse_regkrig$rmse.kfold) #rmse: 3.22; r^2:0.51
clay_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'CLAY', model = "~ elevation + annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m")
mean(clay_0_10_rmse_regkrig$rmse.kfold) #3.07% clay; r2=???
clay_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'CLAY', model= '~ elevation + annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m')
mean(clay_10_30_rmse_regkrig$rmse.kfold) #3.75% clay; r2=0.49
#all clay regression krigging performed worse than ordinary krigging but better than MLR

# WMPD_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_regkrig$rmse.kfold) #0.05617986 mm; r^2=0.71
# WMPD_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_regkrig$rmse.kfold) #0.06184034 mm; r^2=0.59
# WMPD_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_regkrig$rmse.kfold) #0.06532887; r^2=0.66

IC_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'kgIC.m2', model='~ elevation + annual_kwh.m2 + NDVI_2018mean_1m')
mean(IC_0_30_rmse_regkrig$rmse.kfold) #1.03 kg IC m2; r2=0.15
IC_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'kgIC.m2', model='~ elevation + slope + NDVI_2018mean_1m')
mean(IC_0_10_rmse_regkrig$rmse.kfold) #0.320 kg IC m2; r2=0.15
IC_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'kgIC.m2', model='~ annual_kwh.m2 + curvature_mean + NDVI_2017mean_1m')
mean(IC_10_30_rmse_regkrig$rmse.kfold) #0.843 kg IC m2; r2=0.12

# soilP_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'gP.m2')
# mean(soilP_0_30_rmse_regkrig$rmse.kfold) #rmse:1.13; r^2:0
# soilP_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'gP.m2')
# mean(soilP_0_10_rmse_regkrig$rmse.kfold) # 0.825 g soilP m2; r^2:0
# soilP_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'gP.m2')
# mean(soilP_10_30_rmse_regkrig$rmse.kfold) # 0.487 g P soilP m2; r^2=0

#random forest test
#best current lm model: kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp
library(randomForest)
tuneRF(x=as.data.frame(soil_0_30cm_shp)[,c('curvature_mean', 'elevation', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m')], soil_0_30cm_shp$kgOrgC.m2, ntreeTry = 400, stepFactor = 1, improve = 0.02)
RF_kgOrgC_0_30cm <- randomForest(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m + NDVI_2018mean_1m + elevation+ annual_kwh.m2, data = soil_0_30cm_shp, mtry=1) #Mean of squared residuals: 0.3649699
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ RF_kgOrgC_0_30cm$predicted))
RF_kgOrgC_0_30cm$importance

new_grid <- expand.grid(mtry = 1:10)
train_result <- train(x=as.data.frame(soil_0_30cm_shp)[,c('curvature_mean', 'elevation', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m')], y=soil_0_30cm_shp$kgOrgC.m2, method = 'rf', tuneGrid = new_grid)
train_result$results


mean(RF_kgOrgC_0_30cm$mse)
hist(RF_kgOrgC_0_30cm$predicted)
RF_kgOrgC_0_30cm_8var <- randomForest(kgOrgC.m2 ~ NDVI_2017max_1m + curvature_mean + annual_kwh.m2 + slope + elevation + Red_Nov2016 + Red_May2017 + NIR_May2017, data = soil_0_30cm_shp, mtry=2)
hist(RF_kgOrgC_0_30cm_8var$predicted)
hist(RF_kgOrgC_0_30cm$predicted)
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ RF_kgOrgC_0_30cm$predicted)) #r2=0.24
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ RF_kgOrgC_0_30cm_8var$predicted)) #r2=0.31
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ predict(RF_kgOrgC_0_30cm, soil_0_30cm_shp))) #r2=0.90
plot(predict(RF_kgOrgC_0_30cm, soil_0_30cm_shp), soil_0_30cm_shp$kgOrgC.m2)
RF_kgOrgC_0_30cm_clay <- randomForest(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data = soil_0_30cm_shp, mtry=1)
hist(RF_kgOrgC_0_30cm_clay$predicted)
kgOrgC.m2_RF.terrain4_0_30cm <- predict(RF_kgOrgC_0_30cm, Mar2017_terrain_3m)
forage_terrain_energy$kgOrgC.m2_RF.terrain4 <- predict(RF_kgOrgC_0_30cm, forage_terrain_energy)
summary(lm(peak2017 ~ kgOrgC.m2_RF.terrain4, data = forage_terrain_energy)) #Multiple R-squared:  0.1488,	Adjusted R-squared:  0.088 
#F-statistic: 2.447 on 1 and 14 DF,  p-value: 0.14

#cross validate RF
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
tapply(kf, kf, length)
rmse <- rep(NA, 20)
#k <- 1
predictions <- rep(NA, 105)
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  mtry_tuned <- tuneRF(x=as.data.frame(soil_0_30cm_shp)[,c('curvature_mean', 'slope', 'annual_kwh.m2', 'elevation', 'NDVI_2017mean_1m')], soil_0_30cm_shp$kgOrgC.m2, ntreeTry = 100, stepFactor = 1, improve = 0.02)[1]
  RF_kgOrgC_0_30cm <- randomForest(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + NDVI_2017mean_1m, data = trn, mtry=mtry_tuned, ntree=100)
  #orgC_terrain_pred <- predict(Mar2017_terrain_3m, RF_kgOrgC_0_30cm, fun=predict)
  #orgC_tst_pred <- extract(orgC_terrain_pred, tst)
  orgC_tst_pred <- predict(RF_kgOrgC_0_30cm, tst)
  rmse[k] <- RMSE(tst$kgOrgC.m2, orgC_tst_pred)
  predictions[kf==k] <- orgC_tst_pred
}
plot(RF_kgOrgC_0_30cm)
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ predictions)) #r2=0.36
mean(rmse) #0.535269 with clay as predictor
randomForest(CLAY ~ curvature_mean + slope + annual_kwh.m2 + elevation, data = soilC_0_10cm_df, mtry=1)

#random forest CV function[modify so that model selection algorithm can be applied]
# df_pts <- soil_0_30cm_shp
# varname <- 'kgOrgC.m2'
# model <- "~ elevation + slope + annual_kwh.m2 + curvature_mean + NDVI_2017mean_1m + NDVI_2018mean_1m + NIR_meanGS2017 + NIR_meanGS2018"
crossval_RF <- function(df_pts, varname, ntree=75, model='~ curvature_mean + slope + annual_kwh.m2 + elevation') { #varnames=c("curvature_mean", "slope", "annual_kwh.m2", "elevation")
  rmse <- rep(NA, 20) 
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    model_vector <- gsub('[+]', '', model)
    model_vector <- gsub('~ ', '', model_vector)
    model_vector <- unlist(strsplit(model_vector, '  '))
    mtry_tuned <- tuneRF(x=as.data.frame(df_pts)[ ,model_vector], df_pts[[varname]], ntreeTry = 100, stepFactor = 1, improve = 0.02, trace = FALSE)[1]
    #print(mtry_tuned)
    RF_varname <- randomForest(as.formula(paste(varname, model)), data = trn, mtry=mtry_tuned, ntree=ntree)
    varname_tst_pred <- predict(RF_varname, tst)
    rmse[k] <- RMSE(tst[[varname]], varname_tst_pred)
    predictions[kf == k] <- varname_tst_pred
    #plot(RF_varname)
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'kgOrgC.m2', model='~ elevation + curvature_mean + NDVI_2017mean_1m + NDVI_2018mean_1m')
plot(orgC_0_30_rmse_RF$oob.predictions, soil_0_30cm_shp$kgOrgC.m2)
abline(0,1,lty=2)
mean(orgC_0_30_rmse_RF$rmse.kfold) #rmse: 0.514; r^2:0.35

orgC_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'kgOrgC.m2', model='~ curvature_mean + NDVI_2017mean_1m')
mean(orgC_0_10_rmse_RF$rmse.kfold) #0.356; r^2=0.19
orgC_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'kgOrgC.m2', model = '~ curvature_mean + NDVI_2018mean_1m + NIR_meanGS2017 + NIR_meanGS2018')
mean(orgC_10_30_rmse_RF$rmse.kfold) #0.295; r^2=0.43

clay_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'clay_wtd', model= '~ elevation + slope + NIR_meanGS2017 + NIR_meanGS2018')
mean(clay_0_30_rmse_RF$rmse.kfold) #rmse: 3.36; r^2:0.47
clay_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'CLAY', model= '~ slope + annual_kwh.m2 + curvature_mean + NDVI_2018mean_1m + NIR_meanGS2017')
mean(clay_0_10_rmse_RF$rmse.kfold) #3.23% clay; r2=0.36
clay_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'CLAY', model='~ elevation + slope + NIR_meanGS2017 + NIR_meanGS2018')
mean(clay_10_30_rmse_RF$rmse.kfold) #3.79% clay; r2=0.47

# WMPD_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'WMPD_mm')
# mean(WMPD_0_30_rmse_RF$rmse.kfold) #0.05617986 mm; r^2=0.71
# WMPD_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'WMPD_mm')
# mean(WMPD_0_10_rmse_RF$rmse.kfold) #0.06184034 mm; r^2=0.59
# WMPD_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'WMPD_mm')
# mean(WMPD_10_30_rmse_RF$rmse.kfold) #0.06532887; r^2=0.66

IC_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'kgIC.m2', model='~ elevation + slope + curvature_mean + NDVI_2018mean_1m + NIR_meanGS2018')
mean(IC_0_30_rmse_RF$rmse.kfold) #1.04 kg IC m2; r2=0.21
IC_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'kgIC.m2', model='~ elevation + annual_kwh.m2 + NDVI_2017mean_1m + NDVI_2018mean_1m')
mean(IC_0_10_rmse_RF$rmse.kfold) #0.341 kg IC m2; r2=0.09
IC_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'kgIC.m2', model='~ elevation + annual_kwh.m2 + curvature_mean + NDVI_2018mean_1m + NIR_meanGS2018')
mean(IC_10_30_rmse_RF$rmse.kfold) #0.81 kg IC m2; r2=0.21

# soilP_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'gP.m2')
# mean(soilP_0_30_rmse_RF$rmse.kfold) #rmse:1.01; r^2:0
# soilP_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'gP.m2')
# mean(soilP_0_10_rmse_RF$rmse.kfold) #0.74 g soilP m2; r^2:0
# soilP_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'gP.m2')
# mean(soilP_10_30_rmse_RF$rmse.kfold) #0.49 g P soilP m2; r^2=0

#model selection process for RF
# df <- soil_0_30cm_shp
# varname <- 'kgOrgC.m2'
# depth <- '0_30cm'
# varDir <- 'SOC'
model_selection_RF <- function(df, varname, depth, varDir) {
  if(!dir.exists(file.path(modelResults, 'RF_model_selection', varDir))) {
    dir.create(file.path(modelResults, 'RF_model_selection', varDir))
  }
  models_to_test <- expand.grid(elevation=c(TRUE, FALSE), slope=c(TRUE, FALSE), annual_kwh.m2=c(TRUE, FALSE), curvature_mean=c(TRUE, FALSE), NDVI_2017mean_1m=c(TRUE, FALSE), NDVI_2018mean_1m=c(TRUE, FALSE), NIR_meanGS2017=c(TRUE, FALSE), NIR_meanGS2018=c(TRUE, FALSE))
  models_to_test <- models_to_test[1:(nrow(models_to_test)-1), ]
  df_sum <- apply(models_to_test, 1, sum)
  models_to_test <- models_to_test[df_sum>1,]
  model_selection_results <- list(length=nrow(models_to_test))
  #print(length(model_selection_results[1]))
  for (i in 1:nrow(models_to_test)) {
    print(i)
    model_to_test <- paste(colnames(models_to_test)[unlist(models_to_test[i,])], collapse = ' + ')
    model_to_test <- paste('~', model_to_test)
    model_selection_results[[i]] <- crossval_RF(df_pts = df, varname = varname, model = model_to_test)
    print(mean(model_selection_results[[i]]$rmse.kfold))
  }
  mean_RMSEs <- unlist(lapply(model_selection_results, function(x) mean(x$rmse.kfold)))
  mean_RMSEs_all <- do.call(cbind, lapply(model_selection_results, function(x) x$rmse.kfold))
  oob_predictions_all <- do.call(cbind, lapply(model_selection_results, function(x) x$oob.predictions))
  #make summary
  df_model <- apply(models_to_test, 1, sum)
  meanRMSEs <- apply(mean_RMSEs_all, 2, mean)
  oob_r2s <- apply(oob_predictions_all, 2, function(x) summary(lm(df[[varname]] ~ x))$r.squared)
  summary <- data.frame(model=apply(models_to_test, 1, function(x) paste(colnames(models_to_test)[x], collapse = ' + ')), meanRMSE=meanRMSEs, oob_r.squared=oob_r2s, df_model=df_model)
  summary_by_model_df <- split(summary, summary$df_model)
  best_models <- unlist(lapply(summary_by_model_df, function(x) as.character(x$model[which.min(x$meanRMSE)])))
  best_rmses <- unlist(lapply(summary_by_model_df, function(x) min(x$meanRMSE)))
  best_oobs_r2 <- unlist(lapply(summary_by_model_df, function(x) x$oob_r.squared[which.min(x$meanRMSE)]))
  final_summary <- data.frame(model_df=2:8, model_name=best_models, meanRMSE_20foldCV=best_rmses, OOB_r2=best_oobs_r2) #because 1 var models were dropped above
  write.csv(mean_RMSEs_all, file.path(modelResults, 'RF_model_selection', varDir, paste(varname, '_', depth, '_RF_8var_selection_CV_RMSEs.csv')), row.names = FALSE)
  write.csv(oob_predictions_all, file.path(modelResults, 'RF_model_selection', varDir, paste(varname, '_', depth, '_RF_8var_selection_oob_pred.csv')), row.names = FALSE)
  write.csv(models_to_test, file.path(modelResults, 'RF_model_selection', varDir, paste(varname, '_', depth, '_RF_8var_selection_model_test_grid.csv')), row.names = FALSE)
  write.csv(final_summary, file.path(modelResults, 'RF_model_selection', varDir, paste(varname, '_', depth, '_RF_8var_selection_BEST_models.csv')), row.names = FALSE)
  list(RMSEs=mean_RMSEs_all, OOBs=oob_predictions_all, test_grid=models_to_test, best_models=final_summary)
}
orgC_0_30_RF8var <- model_selection_RF(soil_0_30cm_shp, 'kgOrgC.m2', '0_30cm', 'SOC')
orgC_0_30_RF8var$best_models
orgC_0_10_RF8var <- model_selection_RF(soil_0_10cm_shp, 'kgOrgC.m2', '0_10cm', 'SOC')
orgC_0_10_RF8var$best_models
orgC_10_30_RF8var <- model_selection_RF(soil_10_30cm_shp, 'kgOrgC.m2', '10_30cm', 'SOC')
orgC_10_30_RF8var$best_models

clay_0_30_RF8var <- model_selection_RF(soil_0_30cm_shp, 'clay_wtd', '0_30cm', 'Clay')
clay_0_30_RF8var$best_models
clay_0_10_RF8var <- model_selection_RF(soil_0_10cm_shp, 'CLAY', '0_10cm', 'Clay')
clay_0_10_RF8var$best_models
clay_10_30_RF8var <- model_selection_RF(soil_10_30cm_shp, 'CLAY', '10_30cm', 'Clay')
clay_10_30_RF8var$best_models

P_0_30_RF8var <- model_selection_RF(soil_0_30cm_shp, 'gP.m2', '0_30cm', 'OlsenP')
P_0_30_RF8var$best_models
P_0_10_RF8var <- model_selection_RF(soil_0_10cm_shp, 'gP.m2', '0_10cm', 'OlsenP')
P_0_10_RF8var$best_models
P_10_30_RF8var <- model_selection_RF(soil_10_30cm_shp, 'gP.m2', '10_30cm', 'OlsenP')
P_10_30_RF8var$best_models

IC_0_30_RF8var <- model_selection_RF(soil_0_30cm_shp, 'kgIC.m2', '0_30cm', 'SIC')
IC_0_30_RF8var$best_models
IC_0_10_RF8var <- model_selection_RF(soil_0_10cm_shp, 'kgIC.m2', '0_10cm', 'SIC')
IC_0_10_RF8var$best_models
IC_10_30_RF8var <- model_selection_RF(soil_10_30cm_shp, 'kgIC.m2', '10_30cm', 'SIC')
IC_10_30_RF8var$best_models

#export model testing results to csvs
#soil organic carbon results
orgC_0_30_model_comparison_RMSEs <- data.frame(null=orgC_0_30_rmse_null$rmse.kfold, idw=orgC_0_30_rmse_idw$rmse.kfold, nn=orgC_0_30_rmse_nn$rmse.kfold, ordkrig=orgC_0_30_rmse_ordkrig$rmse.kfold, regkrig=orgC_0_30_rmse_regkrig$rmse.kfold, MLR_terrain=orgC_0_30_rmse_lm$rmse.kfold,  RF_terrain=orgC_0_30_rmse_RF$rmse.kfold)
write.csv(orgC_0_30_model_comparison_RMSEs, file.path(modelResults, 'model_comparison_v2', 'orgC_0_30cm_model_comps_RMSEs.csv'), row.names=TRUE)#row.names will be kfold

orgC_0_10_model_comparison_RMSEs <- data.frame(null=orgC_0_10_rmse_null$rmse.kfold, idw=orgC_0_10_rmse_idw$rmse.kfold, nn=orgC_0_10_rmse_nn$rmse.kfold, ordkrig=orgC_0_10_rmse_ordkrig$rmse.kfold, regkrig=orgC_0_10_rmse_regkrig$rmse.kfold, MLR_terrain=orgC_0_10_rmse_lm$rmse.kfold, RF_terrain=orgC_0_10_rmse_RF$rmse.kfold)
write.csv(orgC_0_10_model_comparison_RMSEs, file.path(modelResults, 'model_comparison_v2', 'orgC_0_10cm_model_comps_RMSEs.csv'), row.names=TRUE)

orgC_10_30_model_comparison_RMSEs <- data.frame(null=orgC_10_30_rmse_null$rmse.kfold, idw=orgC_10_30_rmse_idw$rmse.kfold, nn=orgC_10_30_rmse_nn$rmse.kfold, ordkrig=orgC_10_30_rmse_ordkrig$rmse.kfold, regkrig=orgC_10_30_rmse_regkrig$rmse.kfold, MLR_terrain=orgC_10_30_rmse_lm$rmse.kfold, RF_terrain=orgC_10_30_rmse_RF$rmse.kfold)
write.csv(orgC_10_30_model_comparison_RMSEs, file.path(modelResults, 'model_comparison_v2', 'orgC_10_30cm_model_comps_RMSEs.csv'), row.names=TRUE)

orgC_0_30_model_comparison_OOBs <- data.frame(null=orgC_0_30_rmse_null$oob.predictions, idw=orgC_0_30_rmse_idw$oob.predictions, nn=orgC_0_30_rmse_nn$oob.predictions, ordkrig=orgC_0_30_rmse_ordkrig$oob.predictions, regkrig=orgC_0_30_rmse_regkrig$oob.predictions, MLR_terrain=orgC_0_30_rmse_lm$oob.predictions, RF_terrain=orgC_0_30_rmse_RF$oob.predictions)
write.csv(orgC_0_30_model_comparison_OOBs, file.path(modelResults, 'model_comparison_v2', 'orgC_0_30cm_model_comps_OOBs.csv'), row.names=TRUE)#row.names will be kfold

orgC_0_10_model_comparison_OOBs <- data.frame(null=orgC_0_10_rmse_null$oob.predictions, idw=orgC_0_10_rmse_idw$oob.predictions, nn=orgC_0_10_rmse_nn$oob.predictions, ordkrig=orgC_0_10_rmse_ordkrig$oob.predictions, regkrig=orgC_0_10_rmse_regkrig$oob.predictions, lm_4terrain=orgC_0_10_rmse_lm$oob.predictions, MLR_terrain_clay=orgC_0_10_rmse_lm_clay$oob.predictions, RF_4terrain=orgC_0_10_rmse_RF$oob.predictions)
write.csv(orgC_0_10_model_comparison_OOBs, file.path(modelResults, 'orgC_0_10cm_model_comps_OOBs.csv'), row.names=TRUE)

orgC_10_30_model_comparison_OOBs <- data.frame(null=orgC_10_30_rmse_null$oob.predictions, idw=orgC_10_30_rmse_idw$oob.predictions, nn=orgC_10_30_rmse_nn$oob.predictions, ordkrig=orgC_10_30_rmse_ordkrig$oob.predictions, regkrig=orgC_10_30_rmse_regkrig$oob.predictions, lm_4terrain=orgC_10_30_rmse_lm$oob.predictions, lm_3terrain_clay=orgC_10_30_rmse_lm_clay$oob.predictions, RF_4terrain=orgC_10_30_rmse_RF$oob.predictions)
write.csv(orgC_10_30_model_comparison_OOBs, file.path(modelResults, 'orgC_10_30cm_model_comps_OOBs.csv'), row.names=TRUE)

#WMPD results
WMPD_0_30_model_comparison_RMSEs <- data.frame(null=WMPD_0_30_rmse_null$rmse.kfold, idw=WMPD_0_30_rmse_idw$rmse.kfold, nn=WMPD_0_30_rmse_nn$rmse.kfold, ordkrig=WMPD_0_30_rmse_ordkrig$rmse.kfold, regkrig=WMPD_0_30_rmse_regkrig$rmse.kfold, lm_4terrain=WMPD_0_30_rmse_lm$rmse.kfold, RF_4terrain=WMPD_0_30_rmse_RF$rmse.kfold)
write.csv(WMPD_0_30_model_comparison_RMSEs, file.path(modelResults, 'WMPD_0_30cm_model_comps_RMSEs.csv'), row.names=TRUE)#row.names will be kfold

WMPD_0_10_model_comparison_RMSEs <- data.frame(null=WMPD_0_10_rmse_null$rmse.kfold, idw=WMPD_0_10_rmse_idw$rmse.kfold, nn=WMPD_0_10_rmse_nn$rmse.kfold, ordkrig=WMPD_0_10_rmse_ordkrig$rmse.kfold, regkrig=WMPD_0_10_rmse_regkrig$rmse.kfold, lm_4terrain=WMPD_0_10_rmse_lm$rmse.kfold, RF_4terrain=WMPD_0_10_rmse_RF$rmse.kfold)
write.csv(WMPD_0_10_model_comparison_RMSEs, file.path(modelResults, 'WMPD_0_10cm_model_comps_RMSEs.csv'), row.names=TRUE)

WMPD_10_30_model_comparison_RMSEs <- data.frame(null=WMPD_10_30_rmse_null$rmse.kfold, idw=WMPD_10_30_rmse_idw$rmse.kfold, nn=WMPD_10_30_rmse_nn$rmse.kfold, ordkrig=WMPD_10_30_rmse_ordkrig$rmse.kfold, regkrig=WMPD_10_30_rmse_regkrig$rmse.kfold, lm_4terrain=WMPD_10_30_rmse_lm$rmse.kfold, RF_4terrain=WMPD_10_30_rmse_RF$rmse.kfold)
write.csv(WMPD_10_30_model_comparison_RMSEs, file.path(modelResults, 'WMPD_10_30cm_model_comps_RMSEs.csv'), row.names=TRUE)

WMPD_0_30_model_comparison_OOBs <- data.frame(null=WMPD_0_30_rmse_null$oob.predictions, idw=WMPD_0_30_rmse_idw$oob.predictions, nn=WMPD_0_30_rmse_nn$oob.predictions, ordkrig=WMPD_0_30_rmse_ordkrig$oob.predictions, regkrig=WMPD_0_30_rmse_regkrig$oob.predictions, lm_4terrain=WMPD_0_30_rmse_lm$oob.predictions, RF_4terrain=WMPD_0_30_rmse_RF$oob.predictions)
write.csv(WMPD_0_30_model_comparison_OOBs, file.path(modelResults, 'WMPD_0_30cm_model_comps_OOBs.csv'), row.names=TRUE)#row.names will be kfold

WMPD_0_10_model_comparison_OOBs <- data.frame(null=WMPD_0_10_rmse_null$oob.predictions, idw=WMPD_0_10_rmse_idw$oob.predictions, nn=WMPD_0_10_rmse_nn$oob.predictions, ordkrig=WMPD_0_10_rmse_ordkrig$oob.predictions, regkrig=WMPD_0_10_rmse_regkrig$oob.predictions, lm_4terrain=WMPD_0_10_rmse_lm$oob.predictions, RF_4terrain=WMPD_0_10_rmse_RF$oob.predictions)
write.csv(WMPD_0_10_model_comparison_OOBs, file.path(modelResults, 'WMPD_0_10cm_model_comps_OOBs.csv'), row.names=TRUE)

WMPD_10_30_model_comparison_OOBs <- data.frame(null=WMPD_10_30_rmse_null$oob.predictions, idw=WMPD_10_30_rmse_idw$oob.predictions, nn=WMPD_10_30_rmse_nn$oob.predictions, ordkrig=WMPD_10_30_rmse_ordkrig$oob.predictions, regkrig=WMPD_10_30_rmse_regkrig$oob.predictions, lm_4terrain=WMPD_10_30_rmse_lm$oob.predictions, RF_4terrain=WMPD_10_30_rmse_RF$oob.predictions)
write.csv(WMPD_10_30_model_comparison_OOBs, file.path(modelResults, 'WMPD_10_30cm_model_comps_OOBs.csv'), row.names=TRUE)

#soil inorganic carbon results
IC_0_30_model_comparison_RMSEs <- data.frame(null=IC_0_30_rmse_null$rmse.kfold, idw=IC_0_30_rmse_idw$rmse.kfold, nn=IC_0_30_rmse_nn$rmse.kfold, ordkrig=IC_0_30_rmse_ordkrig$rmse.kfold, regkrig=IC_0_30_rmse_regkrig$rmse.kfold, lm_4terrain=IC_0_30_rmse_lm$rmse.kfold, RF_4terrain=IC_0_30_rmse_RF$rmse.kfold)
write.csv(IC_0_30_model_comparison_RMSEs, file.path(modelResults, 'IC_0_30cm_model_comps_RMSEs.csv'), row.names=TRUE)#row.names will be kfold

IC_0_10_model_comparison_RMSEs <- data.frame(null=IC_0_10_rmse_null$rmse.kfold, idw=IC_0_10_rmse_idw$rmse.kfold, nn=IC_0_10_rmse_nn$rmse.kfold, ordkrig=IC_0_10_rmse_ordkrig$rmse.kfold, regkrig=IC_0_10_rmse_regkrig$rmse.kfold, lm_4terrain=IC_0_10_rmse_lm$rmse.kfold, RF_4terrain=IC_0_10_rmse_RF$rmse.kfold)
write.csv(IC_0_10_model_comparison_RMSEs, file.path(modelResults, 'IC_0_10cm_model_comps_RMSEs.csv'), row.names=TRUE)

IC_10_30_model_comparison_RMSEs <- data.frame(null=IC_10_30_rmse_null$rmse.kfold, idw=IC_10_30_rmse_idw$rmse.kfold, nn=IC_10_30_rmse_nn$rmse.kfold, ordkrig=IC_10_30_rmse_ordkrig$rmse.kfold, regkrig=IC_10_30_rmse_regkrig$rmse.kfold, lm_4terrain=IC_10_30_rmse_lm$rmse.kfold, RF_4terrain=IC_10_30_rmse_RF$rmse.kfold)
write.csv(IC_10_30_model_comparison_RMSEs, file.path(modelResults, 'IC_10_30cm_model_comps_RMSEs.csv'), row.names=TRUE)

IC_0_30_model_comparison_OOBs <- data.frame(null=IC_0_30_rmse_null$oob.predictions, idw=IC_0_30_rmse_idw$oob.predictions, nn=IC_0_30_rmse_nn$oob.predictions, ordkrig=IC_0_30_rmse_ordkrig$oob.predictions, regkrig=IC_0_30_rmse_regkrig$oob.predictions, lm_4terrain=IC_0_30_rmse_lm$oob.predictions, RF_4terrain=IC_0_30_rmse_RF$oob.predictions)
write.csv(IC_0_30_model_comparison_OOBs, file.path(modelResults, 'IC_0_30cm_model_comps_OOBs.csv'), row.names=TRUE)#row.names will be kfold

IC_0_10_model_comparison_OOBs <- data.frame(null=IC_0_10_rmse_null$oob.predictions, idw=IC_0_10_rmse_idw$oob.predictions, nn=IC_0_10_rmse_nn$oob.predictions, ordkrig=IC_0_10_rmse_ordkrig$oob.predictions, regkrig=IC_0_10_rmse_regkrig$oob.predictions, lm_4terrain=IC_0_10_rmse_lm$oob.predictions, RF_4terrain=IC_0_10_rmse_RF$oob.predictions)
write.csv(IC_0_10_model_comparison_OOBs, file.path(modelResults, 'IC_0_10cm_model_comps_OOBs.csv'), row.names=TRUE)

IC_10_30_model_comparison_OOBs <- data.frame(null=IC_10_30_rmse_null$oob.predictions, idw=IC_10_30_rmse_idw$oob.predictions, nn=IC_10_30_rmse_nn$oob.predictions, ordkrig=IC_10_30_rmse_ordkrig$oob.predictions, regkrig=IC_10_30_rmse_regkrig$oob.predictions, lm_4terrain=IC_10_30_rmse_lm$oob.predictions, RF_4terrain=IC_10_30_rmse_RF$oob.predictions)
write.csv(IC_10_30_model_comparison_OOBs, file.path(modelResults, 'IC_10_30cm_model_comps_OOBs.csv'), row.names=TRUE)

#soil P results
soilP_0_30_model_comparison_RMSEs <- data.frame(null=soilP_0_30_rmse_null$rmse.kfold, idw=soilP_0_30_rmse_idw$rmse.kfold, nn=soilP_0_30_rmse_nn$rmse.kfold, ordkrig=soilP_0_30_rmse_ordkrig$rmse.kfold, regkrig=soilP_0_30_rmse_regkrig$rmse.kfold, lm_4terrain=soilP_0_30_rmse_lm$rmse.kfold, RF_4terrain=soilP_0_30_rmse_RF$rmse.kfold)
write.csv(soilP_0_30_model_comparison_RMSEs, file.path(modelResults, 'soilP_0_30cm_model_comps_RMSEs.csv'), row.names=TRUE)#row.names will be kfold

soilP_0_10_model_comparison_RMSEs <- data.frame(null=soilP_0_10_rmse_null$rmse.kfold, idw=soilP_0_10_rmse_idw$rmse.kfold, nn=soilP_0_10_rmse_nn$rmse.kfold, ordkrig=soilP_0_10_rmse_ordkrig$rmse.kfold, regkrig=soilP_0_10_rmse_regkrig$rmse.kfold, lm_4terrain=soilP_0_10_rmse_lm$rmse.kfold, RF_4terrain=soilP_0_10_rmse_RF$rmse.kfold)
write.csv(soilP_0_10_model_comparison_RMSEs, file.path(modelResults, 'soilP_0_10cm_model_comps_RMSEs.csv'), row.names=TRUE)

soilP_10_30_model_comparison_RMSEs <- data.frame(null=soilP_10_30_rmse_null$rmse.kfold, idw=soilP_10_30_rmse_idw$rmse.kfold, nn=soilP_10_30_rmse_nn$rmse.kfold, ordkrig=soilP_10_30_rmse_ordkrig$rmse.kfold, regkrig=soilP_10_30_rmse_regkrig$rmse.kfold, lm_4terrain=soilP_10_30_rmse_lm$rmse.kfold, RF_4terrain=soilP_10_30_rmse_RF$rmse.kfold)
write.csv(soilP_10_30_model_comparison_RMSEs, file.path(modelResults, 'soilP_10_30cm_model_comps_RMSEs.csv'), row.names=TRUE)

soilP_0_30_model_comparison_OOBs <- data.frame(null=soilP_0_30_rmse_null$oob.predictions, idw=soilP_0_30_rmse_idw$oob.predictions, nn=soilP_0_30_rmse_nn$oob.predictions, ordkrig=soilP_0_30_rmse_ordkrig$oob.predictions, regkrig=soilP_0_30_rmse_regkrig$oob.predictions, lm_4terrain=soilP_0_30_rmse_lm$oob.predictions, RF_4terrain=soilP_0_30_rmse_RF$oob.predictions)
write.csv(soilP_0_30_model_comparison_OOBs, file.path(modelResults, 'soilP_0_30cm_model_comps_OOBs.csv'), row.names=TRUE)#row.names will be kfold

soilP_0_10_model_comparison_OOBs <- data.frame(null=soilP_0_10_rmse_null$oob.predictions, idw=soilP_0_10_rmse_idw$oob.predictions, nn=soilP_0_10_rmse_nn$oob.predictions, ordkrig=soilP_0_10_rmse_ordkrig$oob.predictions, regkrig=soilP_0_10_rmse_regkrig$oob.predictions, lm_4terrain=soilP_0_10_rmse_lm$oob.predictions, RF_4terrain=soilP_0_10_rmse_RF$oob.predictions)
write.csv(soilP_0_10_model_comparison_OOBs, file.path(modelResults, 'soilP_0_10cm_model_comps_OOBs.csv'), row.names=TRUE)

soilP_10_30_model_comparison_OOBs <- data.frame(null=soilP_10_30_rmse_null$oob.predictions, idw=soilP_10_30_rmse_idw$oob.predictions, nn=soilP_10_30_rmse_nn$oob.predictions, ordkrig=soilP_10_30_rmse_ordkrig$oob.predictions, regkrig=soilP_10_30_rmse_regkrig$oob.predictions, lm_4terrain=soilP_10_30_rmse_lm$oob.predictions, RF_4terrain=soilP_10_30_rmse_RF$oob.predictions)
write.csv(soilP_10_30_model_comparison_OOBs, file.path(modelResults, 'soilP_10_30cm_model_comps_OOBs.csv'), row.names=TRUE)

#regression krigging function development for maps
varname <- 'clay_wtd'
df_pts <- soil_0_30cm_shp
width <- 21
psill <- 85
model <- "Sph"
range <- 75
nugget <- 200
regressionKrig <- function(varname, df_pts, width=21, psill=85, model="Sph", range=75, nugget = 200) {
  varname_lm <- lm(as.formula(paste(varname, '~ curvature_mean + elevation + slope + annual_kwh.m2')), data=df_pts)
  print(summary(varname_lm)) #r^2=0.30
  terrain_pred <- predict(Mar2017_terrain_3m, varname_lm, fun=predict)
  terrain_pred <- resample(terrain_pred, r, method='bilinear')
  df_pts[[paste0(varname, 'lm_residuals')]] <- varname_lm$residuals
#check autocorrelation in residuals
  autocorrtest <- autocorr_test_soil(df_pts, varname, nsim = 999) #p<0.001
  df_pts[[paste0(varname, 'lm_predictions')]] <- varname_lm$fitted.values
  varname_reg_krig <- gstat(formula=as.formula(paste(paste0(varname, 'lm_residuals'), '~1')), locations =  df_pts)
  v <- variogram(varname_reg_krig, width=width)
  plot(v)
  fve <- fit.variogram(v, vgm(psill=85, model="Sph", range=75, nugget = 200))
  print(fve)
  plot(variogramLine(fve, 200), type='l', ylim=c(0,350))
  points(v[,2:3], pch=20, col='red')
  regkrig_model <- gstat(formula=as.formula(paste(paste0(varname, 'lm_residuals'), '~1')), locations = df_pts, model=fve)
# predicted values
  varname_res_regkrig <- predict(regkrig_model, as(r, 'SpatialGrid'))
  varname_res_regkrig <- brick(varname_res_regkrig)
  plot(varname_res_regkrig$var1.pred)
  plot(df_pts[[paste0(varname, 'lm_residuals')]], df_pts[[varname]])
  varname_regkrig_est <- varname_res_regkrig$var1.pred + terrain_pred
  plot(varname_regkrig_est)
  df_pts[[paste0(varname, '_est_regkrig')]] <- extract(varname_regkrig_est, df_pts)
  plot(df_pts[[paste0(varname, '_est_regkrig')]], df_pts[[varname]])
  summary(lm(as.formula(paste(varname, '~', paste0(varname, '_est_regkrig'))), data=df_pts)) #r^2=0.91!
  all_forage_sp[[paste0(varname, '_est_regkrig')]] <- extract(varname_regkrig_est, all_forage_sp)
#plot(all_forage_sp[[paste0(varname, '_est_regkrig')]], all_forage_sp$peak_2017)
#abline(lm(peak_2017 ~ clay_est_regkrig, data = all_forage_sp), lty=2)
  summary(lm(as.formula(paste('peak_2017 ~', paste0(varname, '_est_regkrig'))), data = all_forage_sp)) #r^2<0.01
#plot(all_forage_sp$clay_est_regkrig, all_forage_sp$peak_2018)
#abline(lm(peak_2018 ~ clay_est_regkrig, data = all_forage_sp), lty=2)
  summary(lm(as.formula(paste('peak_2018 ~', paste0(varname, '_est_regkrig'))), data = all_forage_sp))
}

#exploratory models
colnames(soil_0_30cm_df)
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + annual_kwh.m2_norm, data = soil_0_30cm_df)) #only r2=0.26 with interaction offering no improvement
summary(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df))
summary(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m + annual_kwh.m2, data = soil_0_30cm_df))

#best 0-30 cm model
summary(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df)) #r2=0.5
vif(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df))
vif(lm(kgOrgC.m2 ~ curvature_mean_norm + annual_kwh.m2_norm + slope_norm + elevation_norm + NDVI_2017mean_1m_norm, data = soil_0_30cm_df))
plot(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df))


summary(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df[-c(2,82),])) #r2=0.62

#5var predicted vs. observed plot
tiff(file = file.path(FiguresDir, 'predictedSOC_vs_observedSOC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4, 4, 1, 1))
plot(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values, soil_0_30cm_shp$kgOrgC.m2, xlab='', ylab='', pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
#points(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values[c(2, 5, 82)], soil_0_30cm_shp$kgOrgC.m2[c(2, 5, 82)], col='red', pch=19)
abline(0, 1, lty=2)
mtext(text=expression(paste('predicted 0-30 cm SOC (kg', ~m^-2, ')')), side=1, line=2.75)
mtext(text=expression(paste('observed 0-30 cm SOC (kg ', ~m^-2, ')')), side=2, line=2.75)
text(x=2, y=4.8, labels=expression(paste(r^2, '= 0.50, all data')), adj=c(0,0))
text(x=2, y=5.2, labels=paste('5-predictor model'), adj=c(0,0))
#text(x=2, y=4.8,labels=expression(paste(r^2, '= 0.62, red pts removed')), adj=c(0,0))
dev.off()


summary(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)) #r2=0.41
summary(lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm, data = soil_0_30cm_df))
plot(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df))
summary(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df[-c(2,41,82),])) #r2=0.51
tiff(file = file.path(FiguresDir, 'predictedSOC_vs_observedSOC_0_30cm_2var.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=150)
par(mar=c(4, 4, 1, 1))
plot(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values, soil_0_30cm_shp$kgOrgC.m2, xlab='', ylab='', pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
#points(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values[c(2, 5, 82)], soil_0_30cm_shp$kgOrgC.m2[c(2, 5, 82)], col='red', pch=19)
abline(0, 1, lty=2)
mtext(text=expression(paste('predicted 0-30 cm SOC (kg', ~m^-2, ')')), side=1, line=2.75)
mtext(text=expression(paste('observed 0-30 cm SOC (kg ', ~m^-2, ')')), side=2, line=2.75)
#text(x=2, y=5.2, labels=expression(paste(r^2, '= 0.50, all data')), adj=c(0,0))
text(x=2.25, y=5.2, labels=paste('2-predictor model'), adj=c(0,0))
text(x=2.25, y=4.8, labels=expression(paste(r^2, '= 0.41, all data')), adj=c(0,0))
dev.off()