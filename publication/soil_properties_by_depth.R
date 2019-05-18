soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
AnalysisDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/analysis/'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
modelResults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/analysis/model_results'
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

#test mean difference in carbonates by depth
t.test(x=soil_0_10cm_df$inorgC.percent, y=soil_10_30cm_df$inorgC.percent, paired = TRUE)
t.test(x=soil_0_10cm_df$orgC.percent, y=soil_10_30cm_df$orgC.percent, paired = TRUE)
t.test(x=soil_0_10cm_df$totC.percent, y=soil_10_30cm_df$totC.percent, paired = TRUE)
t.test(x=soil_0_10cm_df$totN.percent, y=soil_10_30cm_df$totN.percent, paired = TRUE)
t.test(x=soil_0_10cm_df$CLAY, y=soil_10_30cm_df$CLAY, paired = TRUE)
t.test(x=soil_0_10cm_df$SILT, y=soil_10_30cm_df$SILT, paired = TRUE)
t.test(x=soil_0_10cm_df$SAND, y=soil_10_30cm_df$SAND, paired = TRUE)
t.test(x=soil_0_10cm_df$bulk_density_g_cm3, y=soil_10_30cm_df$bulk_density_g_cm3, paired = TRUE)
t.test(x=soil_0_10cm_df$PH, y=soil_10_30cm_df$PH, paired = TRUE)
t.test(x=soil_0_10cm_df$HCO3_P, y=soil_10_30cm_df$HCO3_P, paired = TRUE)

#create summary stats for 105 pts
summary_stats_0_30cm <- as.data.frame(lapply(as.data.frame(soil_0_30cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'clay_wtd', 'sand_wtd', 'silt_wtd', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'elevation', 'slope', 'annual_kwh.m2', 'curvature_mean', 'NDVI_2017mean_1m', 'NDVI_2018mean_1m')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_0_30cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_0_30cm <- summary_stats_0_30cm[,c(ncol(summary_stats_0_30cm), 1:(ncol(summary_stats_0_30cm)-1))]
summary_stats_0_30cm
write.csv(summary_stats_0_30cm, file.path(AnalysisDir, 'soil property summaries', 'summary_stats_0_30cm.csv'), row.names = FALSE)

summary_stats_0_10cm <- as.data.frame(lapply(as.data.frame(soil_0_10cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_0_10cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_0_10cm <- summary_stats_0_10cm[,c(ncol(summary_stats_0_10cm), 1:(ncol(summary_stats_0_10cm)-1))]
summary_stats_0_10cm
write.csv(summary_stats_0_10cm, file.path(AnalysisDir, 'soil property summaries', 'summary_stats_0_10cm.csv'), row.names = FALSE)

summary_stats_10_30cm <- as.data.frame(lapply(as.data.frame(soil_10_30cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_10_30cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_10_30cm <- summary_stats_10_30cm[,c(ncol(summary_stats_10_30cm), 1:(ncol(summary_stats_10_30cm)-1))]
summary_stats_10_30cm
write.csv(summary_stats_10_30cm, file.path(AnalysisDir, 'soil property summaries', 'summary_stats_10_30cm.csv'), row.names = FALSE)


#concentration summary
summary_stats_0_10cm <- as.data.frame(lapply(as.data.frame(soil_0_10cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'inorgC.percent', 'totC.percent', 'totN.percent', 'CLAY', 'SAND', 'SILT', 'HCO3_P', 'bulk_density_g_cm3', 'PH')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_0_10cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_0_10cm <- summary_stats_0_10cm[,c(ncol(summary_stats_0_10cm), 1:(ncol(summary_stats_0_10cm)-1))]
write.csv(summary_stats_0_10cm, file.path(AnalysisDir, 'soil property summaries', 'summary_stats_0_10cm_conc.csv'), row.names = FALSE)

summary_stats_10_30cm <- as.data.frame(lapply(as.data.frame(soil_10_30cm_shp)[,c('orgC.percent', 'kgOrgC.m2', 'inorgC.percent', 'totC.percent', 'totN.percent', 'CLAY', 'SAND', 'SILT', 'HCO3_P', 'bulk_density_g_cm3', 'PH')], function(x) c(as.numeric(summary(x)), sd(x))))
summary_stats_10_30cm$stat <- c('min', 'Q1', 'median', 'mean', 'Q3', 'max', 'sd')
summary_stats_10_30cm <- summary_stats_10_30cm[,c(ncol(summary_stats_10_30cm), 1:(ncol(summary_stats_10_30cm)-1))]
write.csv(summary_stats_10_30cm, file.path(AnalysisDir, 'soil property summaries', 'summary_stats_10_30cm_conc.csv'), row.names = FALSE)

#soil property spatial autocorrelation
#residuals for 0-30 cm
#elev residuals
lm_elev_0_30cmSOC <- lm(kgOrgC.m2 ~ elevation, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_elev_residuals <- lm_elev_0_30cmSOC$residuals
#annual solrad residuals
lm_solrad_0_30cmSOC <- lm(kgOrgC.m2 ~ annual_kwh.m2, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_solrad_residuals <- lm_solrad_0_30cmSOC$residuals
#slope residuals
lm_slope_0_30cmSOC <- lm(kgOrgC.m2 ~ slope, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_slope_residuals <- lm_slope_0_30cmSOC$residuals
#curv residuals
lm_curv_0_30cmSOC <- lm(kgOrgC.m2 ~ curvature_mean, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_curv_residuals <- lm_curv_0_30cmSOC$residuals
#NDVI 2017 residuals
lm_NDVI2017_0_30cmSOC <- lm(kgOrgC.m2 ~ soil_0_30cm_shp$NDVI_2017mean_1m, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_NDVI2017_residuals <- lm_NDVI2017_0_30cmSOC$residuals
#NDVI 2018 residuals
lm_NDVI2018_0_30cmSOC <- lm(kgOrgC.m2 ~ soil_0_30cm_shp$NDVI_2018mean_1m, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_NDVI2018_residuals <- lm_NDVI2018_0_30cmSOC$residuals
#2 var model residuals
lm_terrain2_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_2var_residuals <- lm_terrain2_0_30cm$residuals
#5 var model residuals
lm_terrain5_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2_5var_residuals <- lm_terrain5_0_30cm$residuals

#residuals for 0-10 cm
lm_terrain5_0_10cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_10cm_shp)
soil_0_10cm_shp$kgOrgC.m2_5var_residuals <- lm_terrain5_0_10cm$residuals
lm_terrain2_0_10cm <- lm(kgOrgC.m2 ~ slope + NDVI_2017mean_1m, data =  soil_0_10cm_shp)
soil_0_10cm_shp$kgOrgC.m2_lm.terrain2 <- extract(kgOrgC.m2_terrain2_0_10cm, soil_0_10cm_shp)
soil_0_10cm_shp$kgOrgC.m2_2var_residuals <- lm_terrain2_0_10cm$residuals
#residuals for 10-30 cm
lm_terrain5_10_30cm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + annual_kwh.m2 + NIR_meanGS2017 + NDVI_2017mean_1m, data =  soil_10_30cm_shp)
soil_10_30cm_shp$kgOrgC.m2_5var_residuals <- lm_terrain5_10_30cm$residuals
lm_terrain2_10_30cm <- lm(kgOrgC.m2 ~ curvature_mean + NIR_meanGS2018, data =  soil_10_30cm_shp)
soil_10_30cm_shp$kgOrgC.m2_2var_residuals <- lm_terrain2_10_30cm$residuals
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
soil_0_30_autocorr <- do.call(rbind, lapply(c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'clay_wtd', 'sand_wtd', 'silt_wtd', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'kgOrgC.m2_2var_residuals', 'kgOrgC.m2_5var_residuals', 'kgOrgC.m2_elev_residuals', 'kgOrgC.m2_curv_residuals', 'kgOrgC.m2_slope_residuals', 'kgOrgC.m2_solrad_residuals', 'kgOrgC.m2_NDVI2017_residuals', 'kgOrgC.m2_NDVI2018_residuals'), function(x) autocorr_test_soil(soil_0_30cm_shp, varname = x, nsim = 999)))
soil_0_30_autocorr
write.csv(soil_0_30_autocorr, file.path(modelResults, 'autocorrelation', 'soil_0_30cm_autocorrelation.csv'), row.names = FALSE)

names(soil_0_10cm_shp)
soil_0_10_autocorr <- do.call(rbind, lapply(c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'kgOrgC.m2_2var_residuals', 'kgOrgC.m2_5var_residuals'), function(x) autocorr_test_soil(soil_0_10cm_shp, varname = x, nsim = 999)))
soil_0_10_autocorr
write.csv(soil_0_10_autocorr, file.path(modelResults, 'autocorrelation', 'soil_0_10cm_autocorrelation.csv'), row.names = FALSE)

soil_10_30_autocorr <- do.call(rbind, lapply(c('orgC.percent', 'kgOrgC.m2', 'kgIC.m2', 'kgTC.m2', 'kgTN.m2', 'CLAY', 'SAND', 'SILT', 'gP.m2', 'bulk_density_g_cm3', 'PH', 'kgOrgC.m2_2var_residuals', 'kgOrgC.m2_5var_residuals'), function(x) autocorr_test_soil(soil_10_30cm_shp, varname = x, nsim = 999)))
soil_10_30_autocorr
write.csv(soil_10_30_autocorr, file.path(modelResults, 'autocorrelation', 'soil_10_30cm_autocorrelation.csv'), row.names = FALSE)