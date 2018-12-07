#spatial analysis
mainDir <- 'C:/Users/smdevine/Desktop/rangeland project'
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data/soil C'
soilDataDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
#0_30 dataset
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'shapefiles', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
library(raster)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))
plot(soil_0_30cm_shp, cex=soil_0_30cm_shp$kgOrgC.m2/2, pch=20)
#0-10 dataset

#10-30 dataset

#now calculate distance from forage sampling points in 2017 to soil sampling points in 2018
waypoint_forage_sp <- shapefile(file.path(forageDir, 'waypoint_forage2017.shp'))
sensor_forage_sp <- shapefile(file.path(forageDir, 'sensor_forage2017.shp'))
names(waypoint_forage_sp)
names(sensor_forage_sp)
all_forage_sp <- rbind(sensor_forage_sp, waypoint_forage_sp)
all_forage_sp <- all_forage_sp[-which(all_forage_sp$location=='B01'),] #outlier near road
all_forage_sp$peak_2017 <- apply(as.data.frame(all_forage_sp)[,2:5], 1, max)
plot(all_forage_sp, cex=all_forage_sp$peak_2017/2000, col='red', pch=2, add=TRUE)
distance_matrix <- as.data.frame(pointDistance(coordinates(all_forage_sp)[,1:2], coordinates(soil_0_30cm_shp)[,1:2], lonlat = FALSE))
colnames(distance_matrix) <- paste('pt_', soil_0_30cm_shp$point_no)
distance_matrix <- cbind(clip_plot=all_forage_sp$location, distance_matrix)
class(distance_matrix)
dim(distance_matrix)
head(distance_matrix)
write.csv(distance_matrix, file.path(soilCresults, 'distance_matrix', 'distance_matrix_clip_plots2017.csv'), row.names = FALSE)
distance_matrix <- read.csv(file.path(soilCresults, 'distance_matrix', 'distance_matrix_clip_plots2017.csv'), stringsAsFactors = FALSE)
pts_prox <- data.frame(pts_less_than=apply(distance_matrix[,2:ncol(distance_matrix)], 1, function(x) sum(x < 20)))
pts_prox$location <- distance_matrix$clip_plot
pts_prox[pts_prox == 0,]

#plot soil C as interpolated map
#see labs 14 and 15 from Quant Geo for tips
#also http://rspatial.org/analysis/rst/4-interpolation.html is more refined source of information
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
rmse <- rep(NA, 20)
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  rmse[k] <- RMSE(trn$kgOrgC.m2, mean(tst$kgOrgC.m2))
}
rmse
mean(rmse) #0.7752293
r <- raster(xmn=(xmin(soil_0_30cm_shp)-20), xmx=(xmax(soil_0_30cm_shp)+20), ymn=(ymin(soil_0_30cm_shp) - 20), ymx=(ymax(soil_0_30cm_shp) + 20), resolution=3, crs=crs(soil_0_30cm_shp))
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
rmse <- rep(NA, 20)
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  gs <- gstat(formula=kgOrgC.m2~1, locations=trn, maxdist = 40, set=list(idp=1.7))
  p <- predict(gs, tst)
  rmse[k] <- RMSE(tst$kgOrgC.m2, p$var1.pred)
}
rmse
mean(rmse) #0.6244324


#cross validate nearest neighbors model
for (k in 1:7) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  gs <- gstat(formula=kgOrgC.m2~1, locations=trn, maxdist = 25, set=list(idp=0))
  p <- predict(gs, tst)
  rmse[k] <- RMSE(tst$kgOrgC.m2, p$var1.pred)
}
rmse
mean(rmse)
# null <- RMSE(soil_0_30cm_shp$kgOrgC.m2, mean(soil_0_30cm_shp$kgOrgC.m2))
# null #0.706


#ordinary kriging
#borrowed code from http://rspatial.org/analysis/rst/4-interpolation.html
orgC_krig <- gstat(formula=kgOrgC.m2~1, locations=soil_0_30cm_shp)
v <- variogram(orgC_krig, width=14)
v
head(v)
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

#cross-validate ordinary krigged
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
tapply(kf, kf, length)
rmse <- rep(NA, 20)
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  orgC_krig <- gstat(formula=kgOrgC.m2~1, locations=trn)
  v <- variogram(orgC_krig, width=13)
  fve <- fit.variogram(v, vgm(psill=0.4, model="Sph", range=50, nugget=0.1))
  print(fve)
  plot(variogramLine(fve, 100), type='l', ylim=c(0,0.7), main=paste0(k, '-fold plot'))
  points(v[,2:3], pch=20, col='red')
  gs <- gstat(formula=kgOrgC.m2~1, locations=trn, model=fve) #used model created above
  p <- predict(gs, tst)
  rmse[k] <- RMSE(tst$kgOrgC.m2, p$var1.pred)
}
rmse
mean(rmse) #0.6219091
soil_0_30cm_shp$orgC_est_krig <- extract(orgC_krigged$var1.pred, soil_0_30cm_shp)
plot(soil_0_30cm_shp$orgC_est_krig, soil_0_30cm_shp$kgOrgC.m2)
summary(lm(kgOrgC.m2 ~ orgC_est_krig, data=soil_0_30cm_df))
sd(soil_0_30cm_shp$kgOrgC.m2) #sd is 0.710 kgC m^-2

#muliple regression prediction
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan'), full.names = TRUE))
names(Mar2017_terrain_3m)
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
Mar2017_terrain_3m$annual_kwh.m2 <- solrad_raster


#krig residuals
#regression krigging model
library(relaimpo)
library(car)
library(gstat)
library(spdep)
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + coords.x1 + coords.x2, data=soil_0_30cm_shp))
orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + annual_kwh.m2, data=soil_0_30cm_shp)
summary(orgC_lm)
vif(orgC_lm)
orgC_terrain_pred <- predict(Mar2017_terrain_3m, orgC_lm, fun=predict)
plot(orgC_terrain_pred)
orgC_terrain_pred <- resample(orgC_terrain_pred, r, method='bilinear')
soil_0_30cm_shp$org_lm_residuals <- orgC_lm$residuals
#check autocorrelation in residuals
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
  results
}
result <- autocorr_test_soil(soil_0_30cm_shp, 'org_lm_residuals', nsim=999)
result #some evidence for autocorrelation in residuals

soil_0_30cm_shp$org_lm_predictions <- orgC_lm$fitted.values
summary(lm(kgOrgC.m2 ~ org_lm_predictions, data = soil_0_30cm_shp))
plot(soil_0_30cm_shp$org_lm_predictions, soil_0_30cm_shp$kgOrgC.m2)
orgC_reg_krig <- gstat(formula=org_lm_residuals~1, locations =  soil_0_30cm_shp)
v <- variogram(orgC_reg_krig, width=21)
plot(v)

fve <- fit.variogram(v, vgm(psill=0.1, model="Sph", range=150, nugget = 0.25)) #finally, convergence!
fve
#model      psill    range
#1   Nug 0.23304615   0.0000
#2   Sph 0.07379567 168.1028
plot(variogramLine(fve, 100), type='l', ylim=c(0,0.7))
points(v[,2:3], pch=20, col='red')
regkrig_model <- gstat(formula=org_lm_residuals ~ 1, locations = soil_0_30cm_shp, model=fve)
# predicted values
orgC_res_regkrig <- predict(regkrig_model, as(r, 'SpatialGrid'))
## [using ordinary kriging]
spplot(orgC_res_regkrig)
orgC_res_regkrig <- brick(orgC_res_regkrig)
plot(orgC_res_regkrig$var1.pred)
plot(soil_0_30cm_shp$org_lm_residuals, soil_0_30cm_shp$kgOrgC.m2)
orgC_regkrig_est <- orgC_res_regkrig$var1.pred + orgC_terrain_pred
plot(orgC_regkrig_est)
soil_0_30cm_shp$orgC_est_regkrig <- extract(orgC_regkrig_est, soil_0_30cm_shp)
plot(soil_0_30cm_shp$orgC_est_regkrig, soil_0_30cm_shp$kgOrgC.m2)
summary(lm(kgOrgC.m2 ~ orgC_est_regkrig, data=soil_0_30cm_shp)) #r^2=0.57
all_forage_sp$orgC_est_regkrig <- extract(orgC_regkrig_est, all_forage_sp)
plot(all_forage_sp$orgC_est_regkrig, all_forage_sp$peak_2017)
abline(lm(peak_2017 ~ orgC_est_regkrig, data = all_forage_sp), lty=2)
summary(lm(peak_2017 ~ orgC_est_regkrig, data = all_forage_sp)) #r^2=0.24

#cross-validate regression krigging
set.seed(20161203)
library(dismo)
kf <- kfold(soil_0_30cm_shp, k=20)
tapply(kf, kf, length)
rmse <- rep(NA, 20)
k <- 1
for (k in 1:20) {
  tst <- soil_0_30cm_shp[kf == k, ]
  trn <- soil_0_30cm_shp[kf != k, ]
  orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope, data=trn)
  trn$residuals <- orgC_lm$residuals
  orgC_terrain_pred <- predict(Mar2017_terrain_3m, orgC_lm, fun=predict)
  orgC_tst_pred <- extract(orgC_terrain_pred, tst)
  orgC_reg_krig <- gstat(formula=residuals ~ 1, locations = trn)
  v <- variogram(orgC_reg_krig, width=15)
  fve <- fit.variogram(v, vgm(psill=0.1, model="Sph", range=150, nugget = 0.2))
  regkrig_model <- gstat(formula=orgC_lm$residuals ~ 1, locations = trn, model=fve)
  p_res_correction <- predict(regkrig_model, tst)
  p <- p_res_correction$var1.pred + orgC_tst_pred
  rmse[k] <- RMSE(tst$kgOrgC.m2, p)
  print(fve)
  plot(variogramLine(fve, 200), type='l', ylim=c(0,0.7), main=paste0(k, '-fold plot'))
  points(v[,2:3], pch=20, col='red')
}
rmse
mean(rmse) #0.5748748


#random forest test
tuneRF(x=as.data.frame(soil_0_10cm_shp)[,c('SAND', 'curvature_mean', 'slope', 'annual_kwh.m2')], soil_0_10cm_shp$orgC.percent, stepFactor = 1)
randomForest(CLAY ~ curvature_mean + slope + annual_kwh.m2 + elevation, data = soilC_0_10cm_df, mtry=1)
tuneRF(x=as.data.frame(soil_0_10cm_shp)[,c('elevation', 'curvature_mean', 'slope', 'annual_kwh.m2')], soil_0_10cm_shp$CLAY, stepFactor = 1.5)
randomForest(SAND ~ curvature_mean + slope + annual_kwh.m2 + elevation, data = soilC_0_10cm_df, mtry=1)
randomForest(x=as.data.frame(soil_0_10cm_shp)[,c('elevation', 'curvature_mean', 'slope', 'annual_kwh.m2')], y=soil_0_10cm_shp$CLAY, mtry=1)
randomForest(x=as.data.frame(soil_0_10cm_shp)[,c('elevation', 'curvature_mean', 'slope', 'annual_kwh.m2')], y=soil_0_10cm_shp$SAND, mtry=1)
head(soilC_0_10cm)
#RF_0_10C <- 
soilC_0_10cm_df <- as.data.frame(soil_0_10cm_shp)
soilC_10_30cm_df <- as.data.frame(soil_10_30cm_shp)
colnames(soilC_0_10cm_df)
RF_0_10C <- randomForest(orgC.percent ~ SAND + curvature_mean + annual_kwh.m2 + slope, data=soilC_0_10cm_df, ntree=500)
RF_10_30C <- randomForest(orgC.percent ~ SAND + curvature_mean + annual_kwh.m2 + slope, data=soilC_10_30cm_df, ntree=500, importance=TRUE)
importance(RF_10_30C)
tuneRF(x=soilC_10_30cm_df[,c('SAND', 'curvature_mean', 'slope', 'annual_kwh.m2')], soilC_10_30cm_df$orgC.percent, mtry=1, improve=0.0005, stepFactor = 1.5)
randomForest(orgC.percent ~ SAND + curvature_mean + annual_kwh.m2 + slope, data=soilC_10_30cm_df, mtry=1, importance=TRUE, ntree=500)
randomForest(x=as.data.frame(soil_0_10cm_shp)[ ,c('SAND', 'curvature_mean', 'slope', 'annual_kwh.m2')], y=soil_0_10cm_shp$orgC.percent, ntree=500)
RF_0_10C
plot(RF_0_10C)
importance(RF_0_10C)
varImpPlot(RF_0_10C)
