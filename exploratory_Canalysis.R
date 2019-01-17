
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