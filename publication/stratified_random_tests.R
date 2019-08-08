#version 1 of this is embedded in the behemoth, soilC_spatial.R
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
FiguresDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/Figures'
ResultsDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/analysis/stratified random tests'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
library(raster)
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
#define functions
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
calc_rand_mean <- function(n) {
  mean(soil_0_30cm_shp$kgOrgC.m2[sample(1:105, n)])
}
iterate_rand_mean <- function(iterations, n) {
  replicate(iterations, calc_rand_mean(n))
}
kmeans_cluster <- function(classes, vars, writetofile=FALSE) {
  km.out.norm <- kmeans(na.omit(terrain_features_3m_df_norm[,vars]), classes) #verified this omits all rows where any var is NA
  catch_clusters <- rep(NA, nrow(terrain_features_3m_df_norm))
  catch_clusters[!is.na(terrain_features_3m_df_norm$NDVI_2017mean_1m_norm)] <- km.out.norm$cluster
  #Mar2017_terrain_3m_cropped$climate_cluster <- catch_clusters
  raster_object <- raster(extent(Mar2017_terrain_3m_cropped), resolution=res(Mar2017_terrain_3m_cropped), crs=crs(Mar2017_terrain_3m_cropped))
  catch_clusters <- setValues(raster_object, catch_clusters)
  if(writetofile) {
    writeRaster(catch_clusters, filename = file.path(FiguresDir, paste0('cluster', classes, '_', length(vars), 'vars.tif')))
  }
  #plot(catch_clusters)
  cluster_ID <- extract(catch_clusters, soil_0_30cm_shp)
  cluster_ID
}
#function to calculate mean from stratified random sampling
# calc_strat_mean_v2 <- function(class_no, classes, sample_no) {
#   class_proportions <- tabulate(soil_0_30cm_shp[[paste0('class_', class_no)]]) / nrow(soil_0_30cm_shp)
#   results <- numeric(length=classes)
#   for (i in seq_along(class_proportions)) {
#     results[i] <- mean(soil_0_30cm_shp$kgOrgC.m2[sample(which(soil_0_30cm_shp[[paste0('class_', class_no)]]==(1:classes)[i]), sample_no[i])]) * class_proportions[i]
#   }
#   sum(results)
# }
#revised function
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

#start with random sampling approach
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
#these above results are synthesized at the end of the script

#read-in terrain properties
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan', '3m_filtered'), full.names = TRUE))
names(Mar2017_terrain_3m)
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
Mar2017_terrain_3m$annual_kwh.m2 <- solrad_raster

NDVI_stack <- stack(list.files(NDVIDir, full.names = TRUE))
names(NDVI_stack)
#cellStats(NDVI_stack, 'mean')
NDVI_2017 <- NDVI_stack[[c(1,3,5,7,9)]]
NDVI_2018 <- NDVI_stack[[c(2,4,6,8)]]
meanNDVI_2017 <- calc(NDVI_2017, fun=mean)#, filename=file.path(FiguresDir, 'meanNDVI2017.tif'))
meanNDVI_2018 <- calc(NDVI_2018, fun=mean)#, filename=file.path(FiguresDir, 'meanNDVI2018.tif'))

#create grid for spatial predictions and upscaling NDVI & reflectance data to same grid
r <- raster(Mar2017_terrain_3m)
#adjuster <- 30
e <- extent(meanNDVI_2017)
#e <- extent(c(xmin(soil_0_30cm_shp)- adjuster, xmax(soil_0_30cm_shp)+adjuster, ymin(soil_0_30cm_shp) - adjuster, ymax(soil_0_30cm_shp) + adjuster))
r <- crop(r, e)

#now do stratified random hypothetical exercise
#set up for kmeans unsupervised classification
Mar2017_terrain_3m_cropped <- crop(Mar2017_terrain_3m, r)
#Mar2017_terrain_3m_cropped$maxNDVI_2017 <- resample(maxNDVI_2017, Mar2017_terrain_3m_cropped)
Mar2017_terrain_3m_cropped$NDVI_2017mean_1m <- resample(meanNDVI_2017, Mar2017_terrain_3m_cropped)
terrain_features_3m_df <- as.data.frame(Mar2017_terrain_3m_cropped)
dim(terrain_features_3m_df)
terrain_features_3m_df <- terrain_features_3m_df[,colnames(terrain_features_3m_df) %in% c('curvature_mean', 'elevation', 'slope', 'annual_kwh.m2', 'NDVI_2017mean_1m')] #'slope', 'annual_kwh.m2'; 3 classes with these works pretty well:'curvature_mean', 'elevation', 'meanNDVI_2017'
dim(terrain_features_3m_df)
terrain_features_3m_df_norm <- as.data.frame(lapply(terrain_features_3m_df, function(x) {(x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE)}))
colnames(terrain_features_3m_df_norm) <- c('curvature_mean_norm', 'elevation_norm', 'slope_norm', 'annual_kwh.m2_norm', 'NDVI_2017mean_1m_norm')

#2-var unsupervised classification

soil_0_30cm_shp$class_2 <- kmeans_cluster(2, c('NDVI_2017mean_1m_norm', 'curvature_mean_norm'))
table(soil_0_30cm_shp$class_2)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2, summary)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2, mean)
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_2, sd)
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

#5-var unsupervised classification
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


#DANGER, order of classes changes randomly, so third argument to calc_strat_mean_v2 needs to be manually checked each time this script is run
#fixed with implementaton of v3 function version
#2 var (normalized mean curv. and mean 2017 NDVI) 2 class kmeans
#round(2*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat2_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(1, 1)))
strat2_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 2))

#round(3*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat3_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(2, 1)))
strat3_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 3))

#round(4*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat4_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(3, 1)))
strat4_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 4))

#round(5*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat5_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(3, 2)))
strat5_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 5))

#round(6*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat6_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(4, 2)))
strat6_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 6))

#round(7*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat7_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(5, 2)))
strat7_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 7))

#round(8*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat8_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(5, 3)))
strat8_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 8))

#round(9*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat9_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(6, 3)))
strat9_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 9))

#round(10*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat10_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(7, 3)))
strat10_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 10))

#round(15*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat15_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(10, 5)))
strat15_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 15))

#round(20*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat20_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(14, 6)))
strat20_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 20))

#round(25*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat25_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(17, 8)))
strat25_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 25))

#round(30*(tabulate(soil_0_30cm_shp$class_2) /105), 1)
#strat30_class2 <- replicate(10000, calc_strat_mean_v2(2, 2, c(21, 9)))
strat30_class2 <- replicate(10000, calc_strat_mean_v3('2', 2, 30))


###2 var (normalized mean curv. and mean 2017 NDVI) 3 class kmeans
#round(3*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat3_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(1, 1, 1)))
strat3_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 3))

#round(4*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat4_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(1, 2, 1)))
strat4_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 4))

#round(5*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat5_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(2, 2, 1)))
strat5_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 5))

#round(6*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat6_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(2, 3, 1)))
strat6_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 6))

#round(7*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat7_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(3, 3, 1)))
strat7_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 7))

#round(8*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat8_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(3, 4, 1)))
strat8_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 8))

#round(9*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat9_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(4, 4, 1)))
strat9_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 9))

#round(10*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat10_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(4, 5, 1)))
strat10_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 10))

#round(12*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat12_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(5, 6, 1)))
strat12_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 12))

#round(15*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat15_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(6, 7, 2)))
strat15_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 15))

#round(18*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat18_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(7, 9, 2)))
strat18_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 18))

#round(20*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat20_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(8, 10, 2)))
strat20_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 20))

#round(21*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat21_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(9, 10, 2)))
strat21_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 21))

#round(24*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat24_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(10, 11, 3)))
strat24_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 24))

#round(25*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat25_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(10, 12, 3)))
strat25_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 25))

#round(27*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat27_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(11, 13, 3)))
strat27_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 27))

#round(30*(tabulate(soil_0_30cm_shp$class_3) /105), 1)
#strat30_class3 <- replicate(10000, calc_strat_mean_v2(3, 3, c(13, 14, 3)))
strat30_class3 <- replicate(10000, calc_strat_mean_v3('3', 3, 30))

###2 var (normalized mean curv. and mean 2017 NDVI) 4 class kmeans
#round(4*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat4_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(1, 1, 1, 1)))
strat4_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 4))

#round(5*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat5_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(2, 1, 1, 1)))
strat5_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 5))

#round(6*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat6_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(2, 1, 1, 2)))
strat6_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 6))

#round(7*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat7_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(3, 1, 1, 2)))
strat7_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 7))

#round(8*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat8_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(3, 1, 2, 2)))
strat8_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 8))

#round(9*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat9_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(3, 1, 2, 3)))
strat9_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 9))

#round(10*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat10_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(4, 1, 2, 3)))
strat10_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 10))

#round(15*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat15_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(6, 1, 3, 5)))
strat15_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 15))

#round(20*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat20_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(8, 1, 4, 7)))
strat20_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 20))

#round(25*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat25_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(10, 2, 5, 8)))
strat25_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 25))

#round(30*(tabulate(soil_0_30cm_shp$class_4) /105), 1)
#strat30_class4 <- replicate(10000, calc_strat_mean_v2(4, 4, c(12, 2, 6, 10)))
strat30_class4 <- replicate(10000, calc_strat_mean_v3('4', 4, 30))

#2 class approach with 5 var
#round(2*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat2_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(1, 1)))
strat2_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 2))

#round(3*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat3_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(1, 2)))
strat3_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 3))

#round(4*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat4_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(1, 3)))
strat4_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 4))

#round(5*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat5_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(2, 3)))
strat5_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 5))

#round(6*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat6_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(2, 4)))
strat6_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 6))

#round(7*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat7_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(2, 5)))
strat7_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 7))

#round(8*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat8_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(3, 5)))
strat8_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 8))

#round(9*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat9_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(3, 6)))
strat9_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 9))

#round(10*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat10_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(3, 7)))
strat10_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 10))

#round(15*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat15_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(5, 10)))
strat15_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 15))

#round(20*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat20_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(7, 13)))
strat20_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 20))

#round(25*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat25_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(8, 17)))
strat25_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 25))

#round(30*(tabulate(soil_0_30cm_shp$class_2_5var) /105), 1)
#strat30_class2_5var <- replicate(10000, calc_strat_mean_v2('2_5var', 2, c(10, 20)))
strat30_class2_5var <- replicate(10000, calc_strat_mean_v3('2_5var', 2, 30))

#3 class approach with 5 var
#round(3*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat3_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 1, 1)))
strat3_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 3))

#round(4*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat4_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 2, 1)))
strat4_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 4))

#round(5*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat5_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 2, 2)))
strat5_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 5))

#round(6*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat6_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(1, 3, 2)))
strat6_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 6))

#round(7*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat7_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(2, 3, 2)))
strat7_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 7))

#round(8*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat8_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(2, 4, 2)))
strat8_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 8))

#round(9*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat9_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(2, 4, 3)))
strat9_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 9))

#round(10*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat10_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(3, 4, 3)))
strat10_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 10))

#round(15*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat15_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(4, 6, 5)))
strat15_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 15))

#round(20*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat20_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(5, 9, 6)))
strat20_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 20))

#round(25*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat25_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(6, 11, 8)))
strat25_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 25))

#round(30*(tabulate(soil_0_30cm_shp$class_3_5var) /105), 1)
#strat30_class3_5var <- replicate(10000, calc_strat_mean_v2('3_5var', 3, c(8, 13, 9)))
strat30_class3_5var <- replicate(10000, calc_strat_mean_v3('3_5var', 3, 30))

#4 class approach with 5 var
#round(4*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat4_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(1, 1, 1, 1)))
strat4_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 4))

#round(5*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat5_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 1, 1, 1)))
strat5_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 5))

#round(6*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat6_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 1, 2, 1)))
strat6_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 6))

#round(7*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat7_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 2, 2, 1)))
strat7_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 7))

#round(8*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat8_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(2, 2, 2, 2)))
strat8_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 8))

#round(9*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat9_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(3, 2, 2, 2)))
strat9_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 9))

#round(10*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat10_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(3, 2, 3, 2)))
strat10_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 10))

#round(15*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat15_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(4, 4, 4, 3)))
strat15_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 15))

#round(20*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat20_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(6, 5, 5, 4)))
strat20_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 20))

#round(25*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat25_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(7, 6, 7, 5)))
strat25_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 25))

#round(30*(tabulate(soil_0_30cm_shp$class_4_5var) /105), 1)
#strat30_class4_5var <- replicate(10000, calc_strat_mean_v2('4_5var', 4, c(8, 8, 8, 6)))
strat30_class4_5var <- replicate(10000, calc_strat_mean_v3('4_5var', 4, 30))

#5 var predictions classified into 5 groups based on the prediction quartiles (21 in 1 to 5, respectively)
soil_0_30cm_shp$SOC_0_30cm_5var.prediction <- lm(kgOrgC.m2 ~ curvature_mean_norm + NDVI_2017mean_1m_norm + elevation_norm + annual_kwh.m2_norm + slope_norm, data=soil_0_30cm_shp)$fitted.values
summary(soil_0_30cm_shp$SOC_0_30cm_5var.prediction)
SOCprediction_quintiles <- quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, c(0.2, 0.4, 0.6, 0.8))
soil_0_30cm_shp$class_SOC_0_30cm_5var.quintile <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quintiles[1], 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quintiles[2], 2, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quintiles[3], 3, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quintiles[4], 4, 5))))
tapply(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$class_SOC_0_30cm_5var.quintile, summary)
table(soil_0_30cm_shp$class_SOC_0_30cm_5var.quintile) #evenly split
strat5_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 5))

strat10_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 10))

strat15_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 15))

strat20_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 20))

strat24_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 24))

strat28_classSOC30_5var.quintile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quintile', 5, 28))

#5 var predictions classified into 4 groups based on the prediction quartiles (26, 26, 26, and 27 in 1 to 4, respectively)
SOCprediction_quartiles <- quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, c(0.25, 0.5, 0.75))
soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quartiles[1], 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quartiles[2], 2, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quartiles[3], 3, 4)))
table(soil_0_30cm_shp$class_SOC_0_30cm_5var.quartile)

strat4_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 4))

strat6_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 6))

strat8_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 8))

strat12_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 12))

strat16_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 16))

strat20_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 20))

strat24_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 24))

strat28_classSOC30_5var.quartile<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.quartile', 4, 28))

#now with the terciles from the 5 var predictions
SOCprediction_quartiles <- quantile(soil_0_30cm_shp$SOC_0_30cm_5var.prediction, c(0.333, 0.666))
soil_0_30cm_shp$class_SOC_0_30cm_5var.terciles <- ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quartiles[1], 1, ifelse(soil_0_30cm_shp$SOC_0_30cm_5var.prediction < SOCprediction_quartiles[2], 2, 3))
table(soil_0_30cm_shp$class_SOC_0_30cm_5var.terciles)

strat3_classSOC30_5var.terciles <- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 3))

strat6_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 6))

strat9_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 9))

strat12_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 12))

strat15_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 15))

strat18_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 18))

strat21_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 21))

strat24_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 24))

strat27_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 27))

strat30_classSOC30_5var.terciles<- replicate(10000, calc_strat_mean_v3('SOC_0_30cm_5var.terciles', 3, 30))

#synthesis of results
results_2class_2var <- data.frame(n=c(2:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat2_class2, strat3_class2, strat4_class2, strat5_class2, strat6_class2, strat7_class2, strat8_class2, strat9_class2, strat10_class2, strat15_class2, strat20_class2, strat25_class2, strat30_class2), calculate_thresholds,  iterations=10000)))
colnames(results_2class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_2class_2var
write.csv(results_2class_2var, file.path(ResultsDir, 'NDVI_curv_two_class.csv'), row.names = FALSE)

results_3class_2var <- data.frame(n=seq(3, 30, 3), do.call(rbind, lapply(list(strat3_class3, strat6_class3, strat9_class3, strat12_class3, strat15_class3, strat18_class3, strat21_class3, strat24_class3, strat27_class3, strat30_class3), calculate_thresholds,  iterations=10000)))
colnames(results_3class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_3class_2var
write.csv(results_3class_2var, file.path(ResultsDir, 'NDVI_curv_three_class.csv'), row.names = FALSE)
results_3class_2var <- read.csv(file.path(ResultsDir, 'NDVI_curv_three_class.csv'), stringsAsFactors = FALSE)

results_3class_MLR.5var <- data.frame(n=seq(from=3, to=30, by=3), do.call(rbind, lapply(list(strat3_classSOC30_5var.terciles, strat6_classSOC30_5var.terciles, strat9_classSOC30_5var.terciles, strat12_classSOC30_5var.quartile, strat15_classSOC30_5var.terciles, strat18_classSOC30_5var.terciles, strat21_classSOC30_5var.terciles, strat24_classSOC30_5var.terciles, strat27_classSOC30_5var.terciles, strat30_classSOC30_5var.terciles), calculate_thresholds,  iterations=10000)))
colnames(results_3class_MLR.5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_3class_MLR.5var
write.csv(results_3class_MLR.5var, file.path(ResultsDir, 'MLRbest5var_three_class.csv'), row.names = FALSE)
results_3class_MLR.5var <- read.csv(file.path(ResultsDir, 'MLRbest5var_three_class.csv'), stringsAsFactors = FALSE)

results_4class_MLR.5var <- data.frame(n=c(4, 8, 12, 16, 20, 24, 28), do.call(rbind, lapply(list(strat4_classSOC30_5var.quartile, strat8_classSOC30_5var.quartile, strat12_classSOC30_5var.quartile, strat16_classSOC30_5var.quartile, strat20_classSOC30_5var.quartile, strat24_classSOC30_5var.quartile, strat28_classSOC30_5var.quartile), calculate_thresholds,  iterations=10000)))
colnames(results_4class_MLR.5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_4class_MLR.5var
write.csv(results_4class_MLR.5var, file.path(ResultsDir, 'MLRbest5var_four_class.csv'), row.names = FALSE)

results_4class_2var <- data.frame(n=c(4:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat4_class4, strat5_class4, strat6_class4, strat7_class4, strat8_class4, strat9_class4, strat10_class4, strat15_class4, strat20_class4, strat25_class4, strat30_class4), calculate_thresholds,  iterations=10000)))
colnames(results_4class_2var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_4class_2var
write.csv(results_4class_2var, file.path(ResultsDir, 'NDVI_curv_four_class.csv'), row.names = FALSE)

results_2class_5var <- data.frame(n=c(2:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat2_class2_5var, strat3_class2_5var, strat4_class2_5var, strat5_class2_5var, strat6_class2_5var, strat7_class2_5var, strat8_class2_5var, strat9_class2_5var, strat10_class2_5var, strat15_class2_5var, strat20_class2_5var, strat25_class2_5var, strat30_class2_5var), calculate_thresholds,  iterations=10000)))
colnames(results_2class_5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_2class_5var
write.csv(results_2class_5var, file.path(ResultsDir, 'FiveVar_two_class.csv'), row.names = FALSE)

results_3class_5var <- data.frame(n=c(3:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat3_class3_5var, strat4_class3_5var, strat5_class3_5var, strat6_class3_5var, strat7_class3_5var, strat8_class3_5var, strat9_class3_5var, strat10_class3_5var, strat15_class3_5var, strat20_class3_5var, strat25_class3_5var, strat30_class3_5var), calculate_thresholds,  iterations=10000)))
colnames(results_3class_5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_3class_5var
write.csv(results_3class_5var, file.path(ResultsDir, 'FiveVar_three_class.csv'), row.names = FALSE)

results_4class_5var <- data.frame(n=c(4:10, 15, 20, 25, 30), do.call(rbind, lapply(list(strat4_class4_5var, strat5_class4_5var, strat6_class4_5var, strat7_class4_5var, strat8_class4_5var, strat9_class4_5var, strat10_class4_5var, strat15_class4_5var, strat20_class4_5var, strat25_class4_5var, strat30_class4_5var), calculate_thresholds,  iterations=10000)))
colnames(results_4class_5var)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_4class_5var
write.csv(results_4class_5var, file.path(ResultsDir, 'FiveVar_four_class.csv'), row.names = FALSE)

results_random_sampling <- data.frame(n=seq(from=3, to=30, by=3), do.call(rbind, lapply(list(sample_3, sample_6, sample_9, sample_12, sample_15, sample_18, sample_21, sample_24, sample_27, sample_30), calculate_thresholds, iterations=10000)))
colnames(results_random_sampling)[2:4] <- c('prob_20%', 'prob_10%', 'prob_5%')
results_random_sampling
write.csv(results_random_sampling, file.path(ResultsDir, 'random_sampling.csv'), row.names = FALSE)
results_random_sampling <- read.csv(file.path(ResultsDir, 'random_sampling.csv'), stringsAsFactors = FALSE)


#plot the comparison (Fig 7)
tiff(file = file.path(FiguresDir, 'Fig7.tif', sep = ''), family = 'Times New Roman', width = 4.5, height = 4.5, pointsize = 11, units = 'in', res=800, compression='lzw')
par(mar=c(4.5, 4.5, 1, 1))
plot(results_random_sampling$n, results_random_sampling$prob_5., ylim=c(0, 1), xlab='', ylab='', type='b', col='lightgrey', lty=1, pch=1)
lines(results_random_sampling$n, results_random_sampling$prob_10., type='b', pch=3, col='lightgrey', lty=2)
lines(results_random_sampling$n, results_random_sampling$prob_20., type='b', pch=16, col='lightgrey', lty=3)
lines(results_3class_2var$n, results_3class_2var$prob_5., type='b', col='black', lty=1, pch=1)
lines(results_3class_2var$n, results_3class_2var$prob_10., type='b', col='black', lty=2, pch=3)
lines(results_3class_2var$n, results_3class_2var$prob_20., type='b', col='black', lty=3, pch=16)
lines(results_3class_MLR.5var$n, results_3class_MLR.5var$prob_5., type='b', col='brown', lty=1, pch=1)
lines(results_3class_MLR.5var$n, results_3class_MLR.5var$prob_10., type='b', col='brown', lty=2, pch=3)
lines(results_3class_MLR.5var$n, results_3class_MLR.5var$prob_20., type='b', col='brown', lty=3, pch=16)
#points(x=2, y=calculate_thresholds(strat2_class2, 10000)[1], col='darkgrey', pch = 8)
#points(x=2, y=calculate_thresholds(strat2_class2, 10000)[2], col='darkgrey', pch = 8)
#points(x=2, y=calculate_thresholds(strat2_class2, 10000)[3], col='darkgrey', pch = 8)
mtext(text='Number of samples', side=1, line=2.75)
mtext(text=paste('Probability to estimate catchment SOC with accuracy of', "\u00B1", 'X%'), side=2, line=2.75)
legend('bottomright', legend = c('stratified (MLR) ± 20% ', 'stratified (k-means) ± 20%', 'random ± 20%', 'stratified (MLR) ± 10%', 'stratified (k-means) ± 10%', 'random ± 10%', 'stratified (MLR) ± 5%', 'stratified (k-means) ± 5%', 'random ± 5%'), col=c('brown', 'black', 'lightgrey', 'brown', 'black', 'lightgrey', 'brown', 'black', 'lightgrey'), lty=c(3, 3, 3, 2, 2, 2, 1, 1, 1), pch=c(16, 16, 16, 3, 3, 3, 1, 1, 1), inset = 0.05)
abline(h=0.9, lty=4, col='black')
dev.off()