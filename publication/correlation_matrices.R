#version 1 of this is embedded in the behemoth, soilC_spatial.R
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
modelResults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/analysis/model_results'
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
library(raster)
#correlation test function
rank_test <- function(x, df, y, mtd) {
  test <- cor.test(x, df[[y]], method = mtd)
  result <- data.frame(col.1=test$p.value, col.2=test$estimate)
  colnames(result) <- c(paste0(y, '.p.val.', mtd), paste0(y, if(mtd=='pearson') {'.cor.'} else {'.rho.'}, mtd))
  result
}
#read-in intermediate results for 0-30 cm, 0-10 cm, and 10-30 cm
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_0_10cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_10_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

#read-in NDVI data
NDVI_stack <- stack(list.files(NDVIDir, full.names = TRUE))
names(NDVI_stack)
# cellStats(NDVI_stack, 'mean')
NDVI_2017 <- NDVI_stack[[c(1,3,5,7,9)]]
NDVI_2018 <- NDVI_stack[[c(2,4,6,8)]]


#for Helen: add monthly NDVI to each data.frame
soil_0_30cm_shp$NDVI_Jan2017mean_1m <- extract(NDVI_2017$Camatta_01162017_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Feb2017mean_1m <- extract(NDVI_2017$Camatta_02152017_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Mar2017mean_1m <- extract(NDVI_2017$Camatta_03172017_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Apr2017mean_1m <- extract(NDVI_2017$Camatta_04062017_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_May2017mean_1m <- extract(NDVI_2017$Camatta_04302017_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Jan2018mean_1m <- extract(NDVI_2018$Camatta_01162018_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Feb2018mean_1m <- extract(NDVI_2018$Camatta_02152018_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Mar2018mean_1m <- extract(NDVI_2018$Camatta_03192018_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
soil_0_30cm_shp$NDVI_Apr2018mean_1m <- extract(NDVI_2018$Camatta_04142018_NDVI, soil_0_30cm_shp, buffer=1, fun=mean)
lapply(list(mean2017='NDVI_2017mean_1m', Jan2017='NDVI_Jan2017mean_1m', Feb2017='NDVI_Feb2017mean_1m', Mar2017='NDVI_Mar2017mean_1m', Apr2017='NDVI_Apr2017mean_1m', May2017='NDVI_May2017mean_1m'), function(x) {
  summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ soil_0_30cm_shp[[x]]))
}
)
lapply(list(mean2018='NDVI_2018mean_1m', Jan2018='NDVI_Jan2018mean_1m', Feb2018='NDVI_Feb2018mean_1m', Mar2018='NDVI_Mar2018mean_1m', Apr2018='NDVI_Apr2018mean_1m'), function(x) {
  summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ soil_0_30cm_shp[[x]]))
}
)

lapply(list(mean2018='NDVI_2018mean_1m', Jan2018='NDVI_Jan2018mean_1m', Feb2018='NDVI_Feb2018mean_1m', Mar2018='NDVI_Mar2018mean_1m', Apr2018='NDVI_Apr2018mean_1m'), function(x) {
  hist(soil_0_30cm_shp[[x]])
}
)
lapply(list(mean2018='NDVI_2018mean_1m', Jan2018='NDVI_Jan2018mean_1m', Feb2018='NDVI_Feb2018mean_1m', Mar2018='NDVI_Mar2018mean_1m', Apr2018='NDVI_Apr2018mean_1m'), function(x) {
  summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ soil_0_30cm_shp[[x]]))
}
)

#for Helen: NDVI cross correlation matrix
names(soil_0_30cm_shp)
mtd_corr <- 'spearman'
NDVIJan2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Jan2017mean_1m', mtd=mtd_corr))
NDVIFeb2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Feb2017mean_1m', mtd=mtd_corr))
NDVIMar2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Mar2017mean_1m', mtd=mtd_corr))
NDVIApr2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Apr2017mean_1m', mtd=mtd_corr))
NDVIApr2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Apr2017mean_1m', mtd=mtd_corr))
NDVIMay2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_May2017mean_1m', mtd=mtd_corr))
NDVIJan2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Jan2018mean_1m', mtd=mtd_corr))
NDVIFeb2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Feb2018mean_1m', mtd=mtd_corr))
NDVIMar2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Mar2018mean_1m', mtd=mtd_corr))
NDVIApr2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[ ,c("NDVI_Jan2017mean_1m", "NDVI_Feb2017mean_1m", "NDVI_Mar2017mean_1m", "NDVI_Apr2017mean_1m", "NDVI_May2017mean_1m", "NDVI_Jan2018mean_1m", "NDVI_Feb2018mean_1m", "NDVI_Mar2018mean_1m", "NDVI_Apr2018mean_1m")]), rank_test, df=soil_0_30cm_shp, y='NDVI_Apr2018mean_1m', mtd=mtd_corr))
#bind 'em
NDVI_corrs <- cbind(NDVIJan2017_corrs, NDVIFeb2017_corrs, NDVIMar2017_corrs, NDVIApr2017_corrs, NDVIMay2017_corrs, NDVIJan2018_corrs, NDVIFeb2018_corrs, NDVIMar2018_corrs, NDVIApr2018_corrs)
write.csv(NDVI_corrs, file.path(modelResults, 'correlations', 'NDVI_pearson_corr_matrix.csv'), row.names=TRUE)

#now soil properties
soil_0_30cm_shp$kgOrgC.m2_0_10cm <- soil_0_10cm_shp$kgOrgC.m2
soil_0_30cm_shp$kgOrgC.m2_10_30cm <- soil_10_30cm_shp$kgOrgC.m2
mtd_corr <- 'spearman'
orgC_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgOrgC.m2', mtd=mtd_corr))
orgC0_10cm_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgOrgC.m2_0_10cm', mtd=mtd_corr))
orgC10_30cm_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgOrgC.m2_10_30cm', mtd=mtd_corr))
BD_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='bulk_density_g_cm3', mtd=mtd_corr))
IC_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='kgIC.m2', mtd=mtd_corr))
clay_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='clay_wtd', mtd=mtd_corr))
elevation_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='elevation', mtd=mtd_corr))
curvmean_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='curvature_mean', mtd=mtd_corr))
solrad_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='annual_kwh.m2', mtd=mtd_corr))
slope_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='slope', mtd=mtd_corr))
NDVI_2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NDVI_2017mean_1m', mtd=mtd_corr))
NDVI_2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NDVI_2018mean_1m', mtd=mtd_corr))
Red_GS_2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='Red_meanGS2017', mtd=mtd_corr))
NIR_GS_2017_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NIR_meanGS2017', mtd=mtd_corr))
Red_GS_2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='Red_meanGS2018', mtd=mtd_corr))
NIR_GS_2018_corrs <- do.call(rbind, lapply(as.data.frame(soil_0_30cm_shp[,c('kgOrgC.m2', 'kgOrgC.m2_0_10cm', 'kgOrgC.m2_10_30cm', 'clay_wtd', 'kgIC.m2', 'bulk_density_g_cm3', 'elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m', 'Red_meanGS2017', 'Red_meanGS2018', 'NIR_meanGS2017', 'NIR_meanGS2018')]), rank_test, df=soil_0_30cm_shp, y='NIR_meanGS2018', mtd=mtd_corr))
predictor_corrs <- cbind(orgC_corrs, orgC0_10cm_corrs, orgC10_30cm_corrs, clay_corrs, IC_corrs, BD_corrs, elevation_corrs, curvmean_corrs, solrad_corrs, slope_corrs, NDVI_2017_corrs, NDVI_2018_corrs, Red_GS_2017_corrs, Red_GS_2018_corrs, NIR_GS_2017_corrs, NIR_GS_2018_corrs)
write.csv(predictor_corrs, file.path(modelResults, 'terrain_soil_corrs_0_30cm_orgCalldepths.csv'), row.names=TRUE) #last version replaced P corrs with bulk density
