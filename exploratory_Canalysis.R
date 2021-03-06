
#multiple linear regression CV test and final map
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + aspect_class + elevation, data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + clay_wtd, data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_0_10cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, data =  soil_10_30cm_shp))

#test some models
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data =  soil_0_30cm_shp)) #r2=0.48
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + CLAY, data =  soil_0_10cm_shp)) #r2=0.23
summary(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + CLAY, data =  soil_10_30cm_shp)) #r2=0.52
summary(lm(kgOrgC.m2 ~ poly(curvature_mean, 2) + poly(slope, 2) + poly(annual_kwh.m2, 2) + poly(clay_wtd, 2), data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ kgTN.m2, data = soil_0_30cm_shp)) #highly correlated: r^2=0.82
summary(lm(kgTN.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data =  soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ gP.m2, data = soil_0_30cm_shp)) #not clearly correlated
plot(lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data =  soil_0_30cm_shp))
lm_terrain3_clay_orgC <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + clay_wtd, data =  soil_0_30cm_shp)

#interaction models
summary(lm(kgOrgC.m2 ~ curvature_mean * clay_wtd, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + clay_wtd, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean * elevation + annual_kwh.m2, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean * slope + annual_kwh.m2 + elevation, data = soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean * slope + annual_kwh.m2 + elevation + NDVI_2017max_1m, data = soil_0_30cm_shp))


#lmterrain3+clay model
tiff(file = file.path(FiguresDir, 'lm_terrain3_clay_vs_orgC_0_30cm.tif', sep = ''), family = 'Times New Roman', width = 3.5, height = 3.5, pointsize = 11, units = 'in', res=150)
par(mar=c(4.5, 4.5, 1, 1))
plot(lm_terrain3_clay_orgC$fitted.values, soil_0_30cm_shp$kgOrgC.m2, xlab='', ylab=expression(paste('soil organic carbon (kg', ' ', m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(0, 1, lty=2)
mtext(text=expression(paste('multiple linear regression fitted values (kg org C ', m^-2, ')')), side=1, line=2.5, at=3.3)
text(x=2.2, y=5.62, labels=expression(paste(r^2, '= 0.48')), adj=c(0, 0))
text(x=2.2, y=5.35, labels=paste('mult. linear reg. model:'), adj=c(0, 0))
text(x=2.2, y=5.05, labels=paste('SOC = mean curv. + aspect'), adj=c(0, 0))
text(x=2.8, y=4.75, labels=paste('+ slope + clay'), adj=c(0, 0))
text(x=2.2, y=4.45, labels='p-val < 0.001', adj=c(0, 0))
dev.off()
summary(lm_terrain3_clay_orgC)

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

#cross-validate orgC regression krigging for 0-30 cm dataset
#turned this into flexible function for multiple depths and varnames
# set.seed(20161203)
# library(dismo)
# kf <- kfold(soil_0_30cm_shp, k=20)
# table(kf)
# rmse <- rep(NA, 20)
# predictions <- rep(NA, 105)
# #kappa list example from autofitVariogram function: c(0.05, seq(0.2, 2, 0.1), 5, 10)
# #k <- 1
# for (k in 1:20) {
#   tst <- soil_0_30cm_shp[kf == k, ]
#   trn <- soil_0_30cm_shp[kf != k, ]
#   orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope + annual_kwh.m2, data=trn)
#   trn$residuals <- orgC_lm$residuals
#   #orgC_terrain_pred <- predict(Mar2017_terrain_3m, orgC_lm, fun=predict)
#   #orgC_tst_pred <- extract(orgC_terrain_pred, tst)
#   #print(summary(varname_lm))
#   orgC_tst_pred <- predict.lm(orgC_lm, tst)
#   soil_0_30cm_shp$kgOrgC.m2.oob.predictions[kf == k] <- orgC_tst_pred
#   orgC_reg_krig <- gstat(formula=residuals ~ 1, locations = trn)
#   v <- variogram(orgC_reg_krig)
#   test <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE)
#   v <- variogram(orgC_reg_krig, boundaries=test$exp_var$dist + 10)
#   #fve <- fit.variogram(v, vgm(c("Exp", "Mat", "Sph", 'Gau', 'Lin', 'Bes', 'Log', 'Ste')), fit.kappa = seq(0.3,30,0.05))
#   fve <- autofitVariogram(formula=residuals~1, input_data = trn, verbose = TRUE, kappa=c(0.05, seq(0.1,5,0.1), seq(5, 250, 5)))
#   fve <- fit.variogram(v, vgm(psill=fve$var_model$psill[2], model=as.character(fve$var_model$model)[2], range=fve$var_model$range[2], nugget = fve$var_model$psill[1], kappa = fve$var_model$kappa[2]))
#   #fve <- fit.variogram(v, vgm(psill=0.3, model="Exp", range=50, nugget = 0.2, kappa = 10))
#   #fve <- fit.variogram(v, vgm(psill=0.3, model="Sph", range=50, nugget = 0.2))
#   regkrig_model <- gstat(formula=residuals ~ 1, locations = trn, model=fve)
#   p_res_correction <- predict(regkrig_model, tst)
#   p <- p_res_correction$var1.pred + orgC_tst_pred
#   rmse[k] <- RMSE(tst$kgOrgC.m2, p)
#   predictions[kf==k] <- p
#   print(fve)
#   plot(variogramLine(fve, 200), type='l', ylim=c(0,0.7), main=paste0(k, '-fold plot'))
#   points(v[,2:3], pch=20, col='red')
# }
# plot(predictions, soil_0_30cm_shp$kgOrgC.m2)
# summary(lm(soil_0_30cm_shp$kgOrgC.m2 ~ predictions)) #OOB prediction is r^2=0.36
# rmse #0.4939221 0.5600666 1.0633800 0.6527512 0.4721887 0.3774705 0.7145091 0.4467885 0.7634694 0.4930032 0.4282195 0.4026534 0.3921021 0.6156166 0.4385753 0.8660794 0.7585025 0.5583349 0.6681401 0.2705451
# mean(rmse) #0.5098548 using v <- variogram(orgC_reg_krig, width = 22) and fve <- fit.variogram(v, vgm(psill=0.3, model="Exp", range=50, nugget = 0.2)), best so far except for multiple linear regression
#0.525458 using autofit function
#0.5266747 using fit.variogram function that considers multiple models

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

#mixed effects test
#linear mixed effects ET vs. P
#doesn't work for some reason
library(nlme)
soil_0_30cm_shp$aspect_class <- as.factor(soil_0_30cm_shp$aspect_class)
table(soil_0_30cm_shp$aspect_class)
lme_orgC_0_30cm <- lme(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation, random = ~1 | aspect_class, data =  as.data.frame(soil_0_30cm_shp[!soil_0_30cm_shp$aspect_class=='east',]))
summary(lme_1)
plot(lme_1$residuals ~ lme_1$fitted)