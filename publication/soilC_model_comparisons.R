#version 1 of this is embedded in the behemoth, soilC_spatial.R
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
FiguresDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/Figures'
ResultsDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/analysis/stratified random tests'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
modelResults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/analysis/model_results'
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
library(raster)
library(gstat)
library(randomForest)
set.seed(20161203)
library(dismo)
kf <- kfold(1:105, k=20)
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
#read-in intermediate results for 0-30 cm, 0-10 cm, and 10-30 cm
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_0_10cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_10_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

#first calculate NULL cross-validated model
#simplest approach to calculate null RMSE
null <- RMSE(soil_0_30cm_shp$kgOrgC.m2, mean(soil_0_30cm_shp$kgOrgC.m2))
null #0.698

#k-fold approach to calculate null RMSE
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
orgC_0_10_rmse_null <- crossval_null(soil_0_10cm_shp, 'kgOrgC.m2')
orgC_10_30_rmse_null <- crossval_null(soil_10_30cm_shp, 'kgOrgC.m2')


#MLR model selection and CV tests
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

model_selection_MLR <- function(df, varname, depth, varDir) {
  if(!dir.exists(file.path(modelResults, 'MLR_model_selection'))) {
    dir.create(file.path(modelResults, 'MLR_model_selection'))
  }
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
  write.csv(mean_RMSEs_all, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_CV_RMSEs.csv')), row.names = FALSE)
  write.csv(oob_predictions_all, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_oob_pred.csv')), row.names = FALSE)
  write.csv(models_to_test, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_model_test_grid.csv')), row.names = FALSE)
  write.csv(final_summary, file.path(modelResults, 'MLR_model_selection', varDir, paste(varname, '_', depth, '_MLR_8var_selection_BEST_models.csv')), row.names = FALSE)
  list(RMSEs=mean_RMSEs_all, OOBs=oob_predictions_all, test_grid=models_to_test, best_models=final_summary)
}
#SOC content
orgC_0_30_MLR8var <- model_selection_MLR(soil_0_30cm_shp, 'kgOrgC.m2', '0_30cm', 'SOC')
orgC_0_30_MLR8var$best_models
orgC_0_10_MLR8var <- model_selection_MLR(soil_0_10cm_shp, 'kgOrgC.m2', '0_10cm', 'SOC')
orgC_0_10_MLR8var$best_models
orgC_10_30_MLR8var <- model_selection_MLR(soil_10_30cm_shp, 'kgOrgC.m2', '10_30cm', 'SOC')
orgC_10_30_MLR8var$best_models

#inverse distance-weighted approach 
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

#cross-validate ordinary kriging
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

#Random Forest cross-validation exercise
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

model_selection_RF <- function(df, varname, depth, varDir) {
  if(!dir.exists(file.path(modelResults, 'RF_model_selection'))) {
    dir.create(file.path(modelResults, 'RF_model_selection'))
  }
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

#export model testing results to csvs
#soil organic carbon results
if(!dir.exists(file.path(modelResults, 'model_comparison'))) {
  dir.create(file.path(modelResults, 'model_comparison'))
}
orgC_0_30_model_comparison_RMSEs <- data.frame(null=orgC_0_30_rmse_null$rmse.kfold, idw=orgC_0_30_rmse_idw$rmse.kfold, nn=orgC_0_30_rmse_nn$rmse.kfold, ordkrig=orgC_0_30_rmse_ordkrig$rmse.kfold, regkrig=orgC_0_30_rmse_regkrig$rmse.kfold, MLR_terrain=orgC_0_30_rmse_lm$rmse.kfold,  RF_terrain=orgC_0_30_rmse_RF$rmse.kfold)
write.csv(orgC_0_30_model_comparison_RMSEs, file.path(modelResults, 'model_comparison', 'orgC_0_30cm_model_comps_RMSEs.csv'), row.names=TRUE)#row.names will be kfold

orgC_0_10_model_comparison_RMSEs <- data.frame(null=orgC_0_10_rmse_null$rmse.kfold, idw=orgC_0_10_rmse_idw$rmse.kfold, nn=orgC_0_10_rmse_nn$rmse.kfold, ordkrig=orgC_0_10_rmse_ordkrig$rmse.kfold, regkrig=orgC_0_10_rmse_regkrig$rmse.kfold, MLR_terrain=orgC_0_10_rmse_lm$rmse.kfold, RF_terrain=orgC_0_10_rmse_RF$rmse.kfold)
write.csv(orgC_0_10_model_comparison_RMSEs, file.path(modelResults, 'model_comparison', 'orgC_0_10cm_model_comps_RMSEs.csv'), row.names=TRUE)

orgC_10_30_model_comparison_RMSEs <- data.frame(null=orgC_10_30_rmse_null$rmse.kfold, idw=orgC_10_30_rmse_idw$rmse.kfold, nn=orgC_10_30_rmse_nn$rmse.kfold, ordkrig=orgC_10_30_rmse_ordkrig$rmse.kfold, regkrig=orgC_10_30_rmse_regkrig$rmse.kfold, MLR_terrain=orgC_10_30_rmse_lm$rmse.kfold, RF_terrain=orgC_10_30_rmse_RF$rmse.kfold)
write.csv(orgC_10_30_model_comparison_RMSEs, file.path(modelResults, 'model_comparison', 'orgC_10_30cm_model_comps_RMSEs.csv'), row.names=TRUE)

orgC_0_30_model_comparison_OOBs <- data.frame(null=orgC_0_30_rmse_null$oob.predictions, idw=orgC_0_30_rmse_idw$oob.predictions, nn=orgC_0_30_rmse_nn$oob.predictions, ordkrig=orgC_0_30_rmse_ordkrig$oob.predictions, regkrig=orgC_0_30_rmse_regkrig$oob.predictions, MLR_terrain=orgC_0_30_rmse_lm$oob.predictions, RF_terrain=orgC_0_30_rmse_RF$oob.predictions)
write.csv(orgC_0_30_model_comparison_OOBs, file.path(modelResults, 'model_comparison', 'orgC_0_30cm_model_comps_OOBs.csv'), row.names=TRUE)#row.names will be kfold

orgC_0_10_model_comparison_OOBs <- data.frame(null=orgC_0_10_rmse_null$oob.predictions, idw=orgC_0_10_rmse_idw$oob.predictions, nn=orgC_0_10_rmse_nn$oob.predictions, ordkrig=orgC_0_10_rmse_ordkrig$oob.predictions, regkrig=orgC_0_10_rmse_regkrig$oob.predictions, lm_4terrain=orgC_0_10_rmse_lm$oob.predictions, MLR_terrain_clay=orgC_0_10_rmse_lm_clay$oob.predictions, RF_4terrain=orgC_0_10_rmse_RF$oob.predictions)
write.csv(orgC_0_10_model_comparison_OOBs, file.path(modelResults, 'model_comparison', 'orgC_0_10cm_model_comps_OOBs.csv'), row.names=TRUE)

orgC_10_30_model_comparison_OOBs <- data.frame(null=orgC_10_30_rmse_null$oob.predictions, idw=orgC_10_30_rmse_idw$oob.predictions, nn=orgC_10_30_rmse_nn$oob.predictions, ordkrig=orgC_10_30_rmse_ordkrig$oob.predictions, regkrig=orgC_10_30_rmse_regkrig$oob.predictions, lm_4terrain=orgC_10_30_rmse_lm$oob.predictions, lm_3terrain_clay=orgC_10_30_rmse_lm_clay$oob.predictions, RF_4terrain=orgC_10_30_rmse_RF$oob.predictions)
write.csv(orgC_10_30_model_comparison_OOBs, file.path(modelResults, 'model_comparison', 'orgC_10_30cm_model_comps_OOBs.csv'), row.names=TRUE)