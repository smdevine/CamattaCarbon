##TO-DO

#spatial analysis
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
list.files(file.path(soilCresults, 'shapefiles'))
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
library(raster)
library(gstat)
library(spdep)
library(automap)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
plot(soil_0_30cm_shp, cex=soil_0_30cm_shp$kgOrgC.m2/2, pch=20)
sd(soil_0_30cm_shp$kgOrgC.m2)/mean(soil_0_30cm_shp$kgOrgC.m2) #CV is 19.5%

#0-10 dataset (modified orgC, TN, clay, IC, and P colnames to match naming conventions for 0-30)
soil_0_10cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
plot(soil_0_10cm_shp, cex=soil_0_10cm_shp$kgOrgC.m2/1.5, pch=20)

#10-30 dataset (modified orgC, TN, clay, IC, and P colnames to match naming conventions for 0-30)
soil_10_30cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
plot(soil_10_30cm_shp, cex=soil_10_30cm_shp$kgOrgC.m2/1.5, pch=20)

#create grid for spatial predictions
r <- raster(xmn=(xmin(soil_0_30cm_shp)-10), xmx=(xmax(soil_0_30cm_shp)+10), ymn=(ymin(soil_0_30cm_shp) - 10), ymx=(ymax(soil_0_30cm_shp) + 10), resolution=3, crs=crs(soil_0_30cm_shp))

#read-in terrain properties
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan'), full.names = TRUE))
names(Mar2017_terrain_3m)
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
plot(solrad_raster)
plot(soil_0_30cm_shp, pch=1, add=TRUE)
Mar2017_terrain_3m$annual_kwh.m2 <- solrad_raster
plot(Mar2017_terrain_3m$curvature_mean)
plot(soil_0_30cm_shp, pch=1, cex=soil_0_30cm_shp$kgOrgC.m2/2, add=TRUE)
plot(Mar2017_terrain_3m$elevation)

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


# distance_matrix <- as.data.frame(pointDistance(coordinates(all_forage_sp)[,1:2], coordinates(soil_0_30cm_shp)[,1:2], lonlat = FALSE))
# colnames(distance_matrix) <- paste('pt_', soil_0_30cm_shp$point_no)
# distance_matrix <- cbind(clip_plot=all_forage_sp$location, distance_matrix)
# class(distance_matrix)
# dim(distance_matrix)
# head(distance_matrix)
# write.csv(distance_matrix, file.path(soilCresults, 'distance_matrix', 'distance_matrix_clip_plots2017.csv'), row.names = FALSE)
# distance_matrix <- read.csv(file.path(soilCresults, 'distance_matrix', 'distance_matrix_clip_plots2017.csv'), stringsAsFactors = FALSE)
# pts_prox <- data.frame(pts_less_than=apply(distance_matrix[,2:ncol(distance_matrix)], 1, function(x) sum(x < 20)))
# pts_prox$location <- distance_matrix$clip_plot
# pts_prox[pts_prox == 0,]

#plot soil C as interpolated map
#see labs 14 and 15 from Quant Geo for tips
#also http://rspatial.org/analysis/rst/4-interpolation.html is more refined source of information

#first look at associations
lapply(as.data.frame(soil_0_30cm_shp[,c('elevation', 'curvature_mean', 'annual_kwh.m2', 'slope', 'clay_wtd')]), function(x) plot(x, soil_0_30cm_shp$kgOrgC.m2))
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd + silt_wtd, data =  soil_0_30cm_shp))

#first calculate NULL cross-validated model
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
#simplest approach to calculate null RMSE
null <- RMSE(soil_0_30cm_shp$kgOrgC.m2, mean(soil_0_30cm_shp$kgOrgC.m2))
null #0.706
#k-fold approach to calculate null RMSE
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
k <- 2
df_pts <- soil_0_30cm_shp
varname <- 'kgOrgC.m2'
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
mean(orgC_0_30_rmse_null$rmse.kfold) #0.671
orgC_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_null$rmse.kfold) #0.472
orgC_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_null$rmse.kfold) #0.452

clay_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_null) #5.2% clay
clay_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_null) #4.5% clay
clay_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_null) #5.9% clay

WMPD_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_null) #0.118 mm
WMPD_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_null) #0.110 mm
WMPD_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_null) #0.129 mm

IC_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_null) #0.876 kg IC m2
IC_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_null) #0.427 kg IC m2
IC_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_null) #0.551 kg IC m2

soilP_0_30_rmse_null <- crossval_null(soil_0_30cm_shp, 'gP.m2')
mean(soilP_0_30_rmse_null) #2.98 g soilP m2
soilP_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'gP.m2')
mean(soilP_0_10_rmse_null) #2.20 g soilP m2
soilP_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'gP.m2')
mean(soilP_10_30_rmse_null) #1.12 g P soilP m2

#multiple linear regression CV test
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + clay_wtd, data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_0_10cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_10_30cm_shp))
crossval_lm <- function(df_pts, varname) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    varname_lm <- lm(as.formula(paste(varname, '~ curvature_mean + slope + annual_kwh.m2 + elevation')), data = trn)
    #print(summary(varname_lm))
    varname_tst_pred <- predict.lm(varname_lm, tst)
    rmse[k] <- RMSE(tst[[varname]], varname_tst_pred)
    predictions[kf == k] <- varname_tst_pred
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'kgOrgC.m2')
mean(orgC_0_30_rmse_lm$rmse.kfold) #0.533 kg orgC m2; r^2=0.36 oob
orgC_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_lm$rmse.kfold) #0.363 kg orgC m2; r^2=0.14 oob
orgC_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_lm$rmse.kfold) #0.291 kg orgC m2; r^2=0.45 oob

clay_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_lm$rmse.kfold) #4.0% clay; r^2=0.24
clay_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_lm$rmse.kfold) #3.5% clay; r^2=0.23
clay_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_lm$rmse.kfold) #4.7% clay; r^2=0.21

WMPD_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_lm$rmse.kfold) #0.094 mm; r^2=0.21
WMPD_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_lm$rmse.kfold) #0.091 mm; r^2=0.16
WMPD_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_lm$rmse.kfold) #0.102 mm; r^2=0.22

IC_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_lm$rmse.kfold) #0.703 kg IC m2; r^2=0.05
IC_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_lm$rmse.kfold) #0.367 kg IC m2; r^2=0
IC_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_lm$rmse.kfold) #0.432 kg IC m2; r^2=0.09

soilP_0_30_rmse_lm <- crossval_lm(soil_0_30cm_shp, 'gP.m2')
mean(soilP_0_30_rmse_lm$rmse.kfold) #2.18 g soilP m2; r^2=0.08
soilP_0_10_rmse_lm <- crossval_lm(soil_0_10cm_shp, 'gP.m2')
mean(soilP_0_10_rmse_lm$rmse.kfold) #1.66 g soilP m2; r^2=0.04
soilP_10_30_rmse_lm <- crossval_lm(soil_10_30cm_shp, 'gP.m2')
mean(soilP_10_30_rmse_lm$rmse.kfold) #0.83 g P soilP m2; r^2=0.19

crossval_lm_clay <- function(df_pts, varname) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    varname_lm <- lm(as.formula(paste(varname, '~ curvature_mean + slope + annual_kwh.m2 + WMPD_mm')), data = trn)
    print(summary(varname_lm))
    varname_tst_pred <- predict.lm(varname_lm, tst)
    rmse[k] <- RMSE(tst[[varname]], varname_tst_pred)
    predictions[kf == k] <- varname_tst_pred
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_lm_clay <- crossval_lm_clay(soil_0_30cm_shp, 'kgOrgC.m2') #r^2=0.41
orgC_0_10_rmse_lm_clay <- crossval_lm_clay(soil_0_10cm_shp, 'kgOrgC.m2') #r^2=0.17 
orgC_10_30_rmse_lm_clay <- crossval_lm_clay(soil_10_30cm_shp, 'kgOrgC.m2') #r^2=0.46
mean(orgC_0_30_rmse_lm_clay$rmse.kfold) #0.5033317; r^2=0.41
mean(orgC_10_30_rmse_lm_clay$rmse.kfold)
t.test(x=orgC_0_30_rmse_lm_clay$rmse.kfold, y=orgC_0_30_rmse_null$rmse.kfold, alternative = 'less', paired = TRUE) #t = -3.9311, df = 19, p-value = 0.0004484
t.test(x=orgC_0_30_rmse_lm_clay$rmse.kfold, y=orgC_0_30_rmse_lm$rmse.kfold, alternative = 'less', paired = TRUE) #t = -1.5417, df = 19, p-value = 0.06981

#multiple linear regression CV test all vars
# summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + curvature_profile  + TCI, data = soil_0_30cm_shp))
# rmse <- rep(NA, 20)
# for (k in 1:20) {
#   tst <- soil_0_30cm_shp[kf == k, ]
#   trn <- soil_0_30cm_shp[kf != k, ]
#   kgOrgC_0_30cm_lm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + curvature_profile  + TCI, data = trn)
#   orgC_tst_pred <- predict.lm(kgOrgC_0_30cm_lm, tst)
#   rmse[k] <- RMSE(tst$kgOrgC.m2, orgC_tst_pred)
# }
# rmse
# mean(rmse) #0.5418224

# #test with clay as predictor
# summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data = soil_0_30cm_shp)) #r^2=0.48
# summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd + kgIC.m2, data = soil_0_30cm_shp)) #r^2=0.48
# rmse <- rep(NA, 20)
# predictions <- rep(NA, 105)
# for (k in 1:20) {
#   tst <- soil_0_30cm_shp[kf == k, ]
#   trn <- soil_0_30cm_shp[kf != k, ]
#   kgOrgC_0_30cm_lm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data = trn)
#   print(summary(kgOrgC_0_30cm_lm))
#   orgC_tst_pred <- predict.lm(kgOrgC_0_30cm_lm, tst)
#   rmse[k] <- RMSE(tst$kgOrgC.m2, orgC_tst_pred)
#   predictions[kf==k] <- orgC_tst_pred
# }
# rmse
# mean(rmse) #0.4994609 only very slight improvement if inorganic carbon or elevation included
# summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ predictions)) #r2=0.43

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

WMPD_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_idw$rmse.kfold) #0.055 mm; r^2=0.73
WMPD_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_idw$rmse.kfold) #0.057 mm; r^2=0.65
WMPD_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_idw$rmse.kfold) #0.064 mm; r^2=0.68

IC_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_idw$rmse.kfold) #0.659 kg IC m2; r^2=0.17
IC_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_idw$rmse.kfold) #0.345 kg IC m2; r^2=0.09
IC_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_idw$rmse.kfold) #0.418 kg IC m2; r^2=0.14

soilP_0_30_rmse_idw <- crossval_idw(soil_0_30cm_shp, 'gP.m2')
mean(soilP_0_30_rmse_idw$rmse.kfold) #2.03 g soilP m2; r^2=0.19
soilP_0_10_rmse_idw <- crossval_idw(soil_0_10cm_shp, 'gP.m2')
mean(soilP_0_10_rmse_idw$rmse.kfold) #1.51 g soilP m2; r^2=0.16
soilP_10_30_rmse_idw <- crossval_idw(soil_10_30cm_shp, 'gP.m2')
mean(soilP_10_30_rmse_idw$rmse.kfold) #0.83 g P soilP m2; r^2=0.17

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

WMPD_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_nn$rmse.kfold) #0.05174892 mm; r^2:0.76
WMPD_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_nn$rmse.kfold) #0.05679778 mm; r^2:0.65
WMPD_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_nn$rmse.kfold) #0.06007388 mm; r^2:0.72

IC_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_nn$rmse.kfold) #0.6083737 kg IC m2; r^2:0.31
IC_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_nn$rmse.kfold) #0.3248757 kg IC m2; r^2:0.21
IC_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_nn$rmse.kfold) #0.4009364 kg IC m2; r^2:0.22

soilP_0_30_rmse_nn <- crossval_nn(soil_0_30cm_shp, 'gP.m2')
mean(soilP_0_30_rmse_nn$rmse.kfold) # 1.91 g soilP m2; r^2:0.29
soilP_0_10_rmse_nn <- crossval_nn(soil_0_10cm_shp, 'gP.m2')
mean(soilP_0_10_rmse_nn$rmse.kfold) #1.43 g soilP m2; r^2:0.26
soilP_10_30_rmse_nn <- crossval_nn(soil_10_30cm_shp, 'gP.m2')
mean(soilP_10_30_rmse_nn$rmse.kfold) #0.79 g P soilP m2; r^2:0.29

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
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
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
mean(clay_0_30_rmse_ordkrig$rmse.kfold) #2.84% clay; r^2=0.63
clay_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_ordkrig$rmse.kfold) #2.88% clay; r2=0.50
clay_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_ordkrig$rmse.kfold) #3.42% clay; r2=0.57

WMPD_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_ordkrig$rmse.kfold) #0.05555203 mm; r^2:0.72
WMPD_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_ordkrig$rmse.kfold) #0.05769833 mm; r^2:0.64
WMPD_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_ordkrig$rmse.kfold) #0.0621748 mm; r^2:0.70

IC_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_ordkrig$rmse.kfold) #0.622131 kg IC m2; r2=0.24
IC_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_ordkrig$rmse.kfold) #0.3295604 kg IC m2; r2=0.16
IC_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_ordkrig$rmse.kfold) #0.4256804 kg IC m2; r2=0.07

soilP_0_30_rmse_ordkrig <- crossval_ordkrig(soil_0_30cm_shp, 'gP.m2', psill = 7, nugget = 0, model = 'Exp', range=25)
mean(soilP_0_30_rmse_ordkrig$rmse.kfold) # 1.96 g soilP m2
soilP_0_10_rmse_ordkrig <- crossval_ordkrig(soil_0_10cm_shp, 'gP.m2', psill = 7, nugget = 0, model = 'Exp', range=25)
mean(soilP_0_10_rmse_ordkrig) #1.47 g soilP m2
soilP_10_30_rmse_ordkrig <- crossval_ordkrig(soil_10_30cm_shp, 'gP.m2', psill = 1, nugget = 0, model = 'Exp', range=25)
mean(soilP_10_30_rmse_ordkrig) #0.79 g P soilP m2

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
soilC_0_30_autocorr <- autocorr_test_soil(soil_0_30cm_shp, 'kgOrgC.m2', nsim = 999)
soilC_0_30_autocorr
#Monte-Carlo simulation of Moran I
soil_0_30_autocorr <- do.call(rbind, lapply(c('kgOrgC.m2', 'kgTN.m2', 'kgClay.m2', 'WMPD_mm', 'sand_wtd', 'silt_wtd', 'clay_wtd', 'kgIC.m2', 'gP.m2'), function(x) autocorr_test_soil(soil_0_30cm_shp, varname = x, nsim = 999)))
soil_0_30_autocorr

names(soil_0_10cm_shp)
soil_0_10_autocorr <- do.call(rbind, lapply(c('kgOrgC.m2', 'kgTN.m2', 'kgClay.m2', 'WMPD_mm', 'SAND', 'SILT', 'CLAY', 'kgIC.m2', 'gP.m2'), function(x) autocorr_test_soil(soil_0_10cm_shp, varname = x, nsim = 999)))
soil_0_10_autocorr

soil_10_30_autocorr <- do.call(rbind, lapply(c('kgOrgC.m2', 'kgTN.m2', 'kgClay.m2', 'WMPD_mm', 'SAND', 'SILT', 'CLAY', 'kgIC.m2', 'gP.m2'), function(x) autocorr_test_soil(soil_10_30cm_shp, varname = x, nsim = 999)))
soil_10_30_autocorr

#regression kriging approach
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + coords.x1 + coords.x2, data=soil_0_30cm_shp))
orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + annual_kwh.m2, data=soil_0_30cm_shp)
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

#cross-validate orgC regression krigging for 0-30 cm dataset
#turned this into flexible function for multiple depths and varnames
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
table(kf)
rmse <- rep(NA, 20)
predictions <- rep(NA, 105)
#kappa list example from autofitVariogram function: c(0.05, seq(0.2, 2, 0.1), 5, 10)
#k <- 1
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + annual_kwh.m2, data=trn)
  trn$residuals <- orgC_lm$residuals
  #orgC_terrain_pred <- predict(Mar2017_terrain_3m, orgC_lm, fun=predict)
  #orgC_tst_pred <- extract(orgC_terrain_pred, tst)
  #print(summary(varname_lm))
  orgC_tst_pred <- predict.lm(orgC_lm, tst)
  soil_0_30cm_shp$kgOrgC.m2.oob.predictions[kf == k] <- orgC_tst_pred
  orgC_reg_krig <- gstat(formula=residuals ~ 1, locations = trn)
  v <- variogram(orgC_reg_krig)
  test <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE)
  v <- variogram(orgC_reg_krig, boundaries=test$exp_var$dist + 10)
  #fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Lin', 'Bes', 'Log', 'Ste')), fit.kappa = seq(0.3,30,0.05))
  fve <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE, kappa=c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
  fve <- fit.variogram(v, vgm(psill=fve$var_model$psill[2], model=as.character(fve$var_model$model)[2], range=fve$var_model$range[2], nugget = fve$var_model$psill[1], kappa = fve$var_model$kappa[2]))
  #fve <- fit.variogram(v, vgm(psill=0.3, model="Exp", range=50, nugget = 0.2, kappa = 10))
  #fve <- fit.variogram(v, vgm(psill=0.3, model="Sph", range=50, nugget = 0.2))
  regkrig_model <- gstat(formula=residuals ~ 1, locations = trn, model=fve)
  p_res_correction <- predict(regkrig_model, tst)
  p <- p_res_correction$var1.pred + orgC_tst_pred
  rmse[k] <- RMSE(tst$kgOrgC.m2, p)
  predictions[kf==k] <- p
  print(fve)
  plot(variogramLine(fve, 200), type='l', ylim=c(0,0.7), main=paste0(k, '-fold plot'))
  points(v[,2:3], pch=20, col='red')
}
plot(predictions, soil_0_30cm_shp$kgOrgC.m2)
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ predictions)) #OOB prediction is r^2=0.36
rmse #0.4939221 0.5600666 1.0633800 0.6527512 0.4721887 0.3774705 0.7145091 0.4467885 0.7634694 0.4930032 0.4282195 0.4026534 0.3921021 0.6156166 0.4385753 0.8660794 0.7585025 0.5583349 0.6681401 0.2705451
mean(rmse) #0.5098548 using v <- variogram(orgC_reg_krig, width = 22) and fve <- fit.variogram(v, vgm(psill=0.3, model="Exp", range=50, nugget = 0.2)), best so far except for multiple linear regression
#0.525458 using autofit function
#0.5266747 using fit.variogram function that considers multiple models
#turn this into function
#0.523046 using 
crossval_regkrig <- function(df_pts, varname) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    varname_lm <- lm(as.formula(paste(varname, '~ curvature_mean + elevation + slope + annual_kwh.m2')), data=trn)
    trn$residuals <- varname_lm$residuals
    # terrain_pred <- predict(Mar2017_terrain_3m, varname_lm, fun=predict)
    # tst_pred <- extract(terrain_pred, tst)
    tst_pred <- predict.lm(varname_lm, tst)
    reg_krig <- gstat(formula=residuals ~ 1, locations = trn)
    test <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE)
    #v <- variogram(reg_krig)
    v <- variogram(reg_krig, boundaries=test$exp_var$dist + 10)
    fve <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE, kappa=c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
    fve <- fit.variogram(v, vgm(psill=fve$var_model$psill[2], model=as.character(fve$var_model$model)[2], range=fve$var_model$range[2], nugget = fve$var_model$psill[1], kappa = fve$var_model$kappa[2]))
    #fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Ste')), fit.kappa = c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
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

orgC_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'kgOrgC.m2')
plot(orgC_0_30_rmse_regkrig$oob.predictions, soil_0_30cm_shp$kgOrgC.m2)
abline(0,1,lty=2)
mean(orgC_0_30_rmse_regkrig$rmse.kfold) #rmse: 0.5278; r^2:0.38
orgC_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_regkrig$rmse.kfold) #0.363; r^2=0.16
orgC_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_regkrig$rmse.kfold) #0.284; r^2=0.48

clay_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_regkrig$rmse.kfold) #rmse: 2.86; r^2:0.62
clay_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_regkrig$rmse.kfold) #2.83% clay; r2=0.49
clay_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_regkrig$rmse.kfold) #3.44% clay; r2=0.57

WMPD_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_regkrig$rmse.kfold) #0.05617986 mm; r^2=0.71
WMPD_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_regkrig$rmse.kfold) #0.06184034 mm; r^2=0.59
WMPD_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_regkrig$rmse.kfold) #0.06532887; r^2=0.66

IC_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_regkrig$rmse.kfold) #0.6593583 kg IC m2; r2=0.14
IC_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_regkrig$rmse.kfold) #0.3390994 kg IC m2; r2=0.10
IC_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_regkrig$rmse.kfold) #0.4339998 kg IC m2; r2=0.07

soilP_0_30_rmse_regkrig <- crossval_regkrig(soil_0_30cm_shp, 'gP.m2')
mean(soilP_0_30_rmse_regkrig$rmse.kfold) #rmse:1.9197; r^2:0.22
soilP_0_10_rmse_regkrig <- crossval_regkrig(soil_0_10cm_shp, 'gP.m2')
mean(soilP_0_10_rmse_regkrig$rmse.kfold) #1.52 g soilP m2; r^2:0.17
soilP_10_30_rmse_regkrig <- crossval_regkrig(soil_10_30cm_shp, 'gP.m2')
mean(soilP_10_30_rmse_regkrig$rmse.kfold) #0.81 g P soilP m2; r^2=0.20

#export model testing results to csvs
#soil organic carbon results
orgC_0_30_model_comparison_RMSEs <- data.frame(null=orgC_0_30_rmse_null$rmse.kfold, idw=orgC_0_30_rmse_idw$rmse.kfold, nn=orgC_0_30_rmse_nn$rmse.kfold, ordkrig=orgC_0_30_rmse_ordkrig$rmse.kfold, regkrig=orgC_0_30_rmse_regkrig$rmse.kfold, lm_4terrain=orgC_0_30_rmse_lm$rmse.kfold, lm_3terrain_clay=orgC_0_30_rmse_lm_clay$rmse.kfold, RF_4terrain=orgC_0_30_rmse_RF$rmse.kfold)
write.csv(orgC_0_30_model_comparison_RMSEs, file.path(modelResults, 'orgC_0_30cm_model_comps_RMSEs.csv'), row.names=TRUE)#row.names will be kfold

orgC_0_10_model_comparison_RMSEs <- data.frame(null=orgC_0_10_rmse_null$rmse.kfold, idw=orgC_0_10_rmse_idw$rmse.kfold, nn=orgC_0_10_rmse_nn$rmse.kfold, ordkrig=orgC_0_10_rmse_ordkrig$rmse.kfold, regkrig=orgC_0_10_rmse_regkrig$rmse.kfold, lm_4terrain=orgC_0_10_rmse_lm$rmse.kfold, lm_3terrain_clay=orgC_0_10_rmse_lm_clay$rmse.kfold, RF_4terrain=orgC_0_10_rmse_RF$rmse.kfold)
write.csv(orgC_0_10_model_comparison_RMSEs, file.path(modelResults, 'orgC_0_10cm_model_comps_RMSEs.csv'), row.names=TRUE)

orgC_10_30_model_comparison_RMSEs <- data.frame(null=orgC_10_30_rmse_null$rmse.kfold, idw=orgC_10_30_rmse_idw$rmse.kfold, nn=orgC_10_30_rmse_nn$rmse.kfold, ordkrig=orgC_10_30_rmse_ordkrig$rmse.kfold, regkrig=orgC_10_30_rmse_regkrig$rmse.kfold, lm_4terrain=orgC_10_30_rmse_lm$rmse.kfold, lm_3terrain_clay=orgC_10_30_rmse_lm_clay$rmse.kfold, RF_4terrain=orgC_10_30_rmse_RF$rmse.kfold)
write.csv(orgC_10_30_model_comparison_RMSEs, file.path(modelResults, 'orgC_10_30cm_model_comps_RMSEs.csv'), row.names=TRUE)

orgC_0_30_model_comparison_OOBs <- data.frame(null=orgC_0_30_rmse_null$oob.predictions, idw=orgC_0_30_rmse_idw$oob.predictions, nn=orgC_0_30_rmse_nn$oob.predictions, ordkrig=orgC_0_30_rmse_ordkrig$oob.predictions, regkrig=orgC_0_30_rmse_regkrig$oob.predictions, lm_4terrain=orgC_0_30_rmse_lm$oob.predictions, lm_3terrain_clay=orgC_0_30_rmse_lm_clay$oob.predictions, RF_4terrain=orgC_0_30_rmse_RF$oob.predictions)
write.csv(orgC_0_30_model_comparison_OOBs, file.path(modelResults, 'orgC_0_30cm_model_comps_OOBs.csv'), row.names=TRUE)#row.names will be kfold

orgC_0_10_model_comparison_OOBs <- data.frame(null=orgC_0_10_rmse_null$oob.predictions, idw=orgC_0_10_rmse_idw$oob.predictions, nn=orgC_0_10_rmse_nn$oob.predictions, ordkrig=orgC_0_10_rmse_ordkrig$oob.predictions, regkrig=orgC_0_10_rmse_regkrig$oob.predictions, lm_4terrain=orgC_0_10_rmse_lm$oob.predictions, lm_3terrain_clay=orgC_0_10_rmse_lm_clay$oob.predictions, RF_4terrain=orgC_0_10_rmse_RF$oob.predictions)
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

#random forest test
library(randomForest)
tuneRF(x=as.data.frame(soil_0_30cm_shp)[,c('curvature_mean', 'slope', 'annual_kwh.m2', 'elevation')], soil_0_30cm_shp$kgOrgC.m2, ntreeTry = 100, stepFactor = 1, improve = 0.02)
RF_kgOrgC_0_30cm <- randomForest(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data = soil_0_30cm_shp, mtry=1) #Mean of squared residuals: 0.3649699
summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ predict(RF_kgOrgC_0_30cm, soil_0_30cm_shp))) #r2=0.90
plot(predict(RF_kgOrgC_0_30cm, soil_0_30cm_shp), soil_0_30cm_shp$kgOrgC.m2)
RF_kgOrgC_0_30cm_clay <- randomForest(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data = soil_0_30cm_shp, mtry=1)

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
  RF_kgOrgC_0_30cm <- randomForest(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data = trn, mtry=2, ntree=100)
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

#random forest CV function
crossval_RF <- function(df_pts, varname, mtry=1, ntree=75) {
  rmse <- rep(NA, 20)
  predictions <- rep(NA, 105)
  for (k in 1:20) {
    tst <- df_pts[kf == k, ]
    trn <- df_pts[kf != k, ]
    RF_varname <- randomForest(as.formula(paste(varname, '~ curvature_mean + slope + annual_kwh.m2 + elevation')), data = trn, mtry=mtry, ntree=ntree)
    varname_tst_pred <- predict(RF_varname, tst)
    rmse[k] <- RMSE(tst[[varname]], varname_tst_pred)
    predictions[kf == k] <- varname_tst_pred
    plot(RF_varname)
  }
  print(summary(lm(df_pts[[varname]] ~ predictions)))
  list(rmse.kfold=rmse, oob.predictions=predictions)
}
orgC_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'kgOrgC.m2')
plot(orgC_0_30_rmse_RF$oob.predictions, soil_0_30cm_shp$kgOrgC.m2)
abline(0,1,lty=2)
mean(orgC_0_30_rmse_RF$rmse.kfold) #rmse: 0.5278; r^2:0.38
orgC_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'kgOrgC.m2')
mean(orgC_0_10_rmse_RF$rmse.kfold) #0.363; r^2=0.16
orgC_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'kgOrgC.m2')
mean(orgC_10_30_rmse_RF$rmse.kfold) #0.284; r^2=0.48

clay_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'clay_wtd')
mean(clay_0_30_rmse_RF$rmse.kfold) #rmse: 2.86; r^2:0.62
clay_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'CLAY')
mean(clay_0_10_rmse_RF$rmse.kfold) #2.83% clay; r2=0.49
clay_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'CLAY')
mean(clay_10_30_rmse_RF$rmse.kfold) #3.44% clay; r2=0.57

WMPD_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'WMPD_mm')
mean(WMPD_0_30_rmse_RF$rmse.kfold) #0.05617986 mm; r^2=0.71
WMPD_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'WMPD_mm')
mean(WMPD_0_10_rmse_RF$rmse.kfold) #0.06184034 mm; r^2=0.59
WMPD_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'WMPD_mm')
mean(WMPD_10_30_rmse_RF$rmse.kfold) #0.06532887; r^2=0.66

IC_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'kgIC.m2')
mean(IC_0_30_rmse_RF$rmse.kfold) #0.6593583 kg IC m2; r2=0.14
IC_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'kgIC.m2')
mean(IC_0_10_rmse_RF$rmse.kfold) #0.3390994 kg IC m2; r2=0.10
IC_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'kgIC.m2')
mean(IC_10_30_rmse_RF$rmse.kfold) #0.4339998 kg IC m2; r2=0.07

soilP_0_30_rmse_RF <- crossval_RF(soil_0_30cm_shp, 'gP.m2')
mean(soilP_0_30_rmse_RF$rmse.kfold) #rmse:1.9197; r^2:0.22
soilP_0_10_rmse_RF <- crossval_RF(soil_0_10cm_shp, 'gP.m2')
mean(soilP_0_10_rmse_RF$rmse.kfold) #1.52 g soilP m2; r^2:0.17
soilP_10_30_rmse_RF <- crossval_RF(soil_10_30cm_shp, 'gP.m2')
mean(soilP_10_30_rmse_RF$rmse.kfold) #0.81 g P soilP m2; r^2=0.20
