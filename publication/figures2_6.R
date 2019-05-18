terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
FiguresDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/Geoderma publication/Figures'
NDVIDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/NDVI'
spatialforageDir <- 'C:/Users/smdevine/Desktop/rangeland project/clip_plots/results'
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win') #once per session
library(raster)
res_plot <- 800

#for kmeans classification Fig 6
kmeans_cluster <- function(classes, vars, writetofile=FALSE, fname) {
  km.out.norm <- kmeans(na.omit(terrain_features_3m_df_norm[,vars]), classes) #verified this omits all rows where any var is NA
  catch_clusters <- rep(NA, nrow(terrain_features_3m_df_norm))
  catch_clusters[!is.na(terrain_features_3m_df_norm$NDVI_2017mean_1m_norm)] <- km.out.norm$cluster
  #Mar2017_terrain_3m_cropped$climate_cluster <- catch_clusters
  raster_object <- raster(extent(Mar2017_terrain_3m_cropped), resolution=res(Mar2017_terrain_3m_cropped), crs=crs(Mar2017_terrain_3m_cropped))
  catch_clusters <- setValues(raster_object, catch_clusters)
  if(writetofile) {
    writeRaster(catch_clusters, filename = file.path(FiguresDir, fname))
  }
  #plot(catch_clusters)
  cluster_ID <- extract(catch_clusters, soil_0_30cm_shp)
  cluster_ID
}

#read-in intermediate results for 0-30 cm, 0-10 cm, and 10-30 cm
soil_0_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_30cm_df.csv'), stringsAsFactors = FALSE)
soil_0_30cm_shp <- SpatialPointsDataFrame(soil_0_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_0_10cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_0_10cm_df.csv'), stringsAsFactors = FALSE)
soil_0_10cm_shp <- SpatialPointsDataFrame(soil_0_10cm_df[,c('coords.x1', 'coords.x2')], data=soil_0_10cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

soil_10_30cm_df <- read.csv(file.path(soilCresults, 'intermediate_results', 'soil_10_30cm_df.csv'), stringsAsFactors = FALSE)
soil_10_30cm_shp <- SpatialPointsDataFrame(soil_10_30cm_df[,c('coords.x1', 'coords.x2')], data=soil_10_30cm_df, proj4string = CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

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

#set up for kmeans unsupervised classification Fig6
Mar2017_terrain_3m_cropped <- crop(Mar2017_terrain_3m, r)
#Mar2017_terrain_3m_cropped$maxNDVI_2017 <- resample(maxNDVI_2017, Mar2017_terrain_3m_cropped)
Mar2017_terrain_3m_cropped$NDVI_2017mean_1m <- resample(meanNDVI_2017, Mar2017_terrain_3m_cropped)

#make plots of direct associations (Fig 2a-f)
#Fig2a
tiff(file = file.path(FiguresDir, 'Fig2a.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$elevation, soil_0_30cm_shp$kgOrgC.m2, xlab=paste('Elevation (m)'), ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #ylim = c(300, 1600), xlim=c(1400, 4700), col=soil_0_30cm_shp$energy_colors
abline(lm(kgOrgC.m2 ~ elevation, data = soil_0_30cm_shp), lty=2)
text(x=475, y=2.5, labels=expression(paste(R^2, '= 0.11')), adj=c(0,0))
text(x=475, y=2.1,labels=paste('p value < 0.001'), adj=c(0,0))
text(x=476, y=5.4, labels='a')
dev.off()
summary(lm(kgOrgC.m2 ~ elevation, data = soil_0_30cm_shp)) #r2=0.11
summary(lm(kgOrgC.m2 ~ elevation, data = soil_0_10cm_shp)) #NS p.val=0.13
summary(lm(kgOrgC.m2 ~ elevation, data = soil_10_30cm_shp)) #r2=0.15

#Fig2b
tiff(file = file.path(FiguresDir, 'Fig2b.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$curvature_mean, soil_0_30cm_shp$kgOrgC.m2, xlab=paste('Mean curvature'), ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ curvature_mean, data = soil_0_30cm_shp), lty=2)
text(x=-1.5, y=2.5, labels=expression(paste(R^2, '= 0.24')), adj=c(0,0))
text(x=-1.5, y=2.1,labels=paste('p value < 0.001'),  adj=c(0,0))
text(x=-1.5, y=5.4, labels='b')
dev.off()
summary(lm(kgOrgC.m2 ~ curvature_mean, data = soil_0_30cm_shp)) #r2=0.24
summary(lm(kgOrgC.m2 ~ curvature_mean, data = soil_0_10cm_shp)) #r2=0.06
summary(lm(kgOrgC.m2 ~ curvature_mean, data = soil_10_30cm_shp)) #r2=0.32

#Fig2c
tiff(file = file.path(FiguresDir, 'Fig2c.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$NDVI_2017mean_1m, soil_0_30cm_shp$kgOrgC.m2, xlab='Mean 2017 NDVI', ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ NDVI_2017mean_1m, data = soil_0_30cm_shp), lty=2)
text(x=0.43, y=5.4, labels='c')
text(x=0.42, y=4.9, labels=expression(paste(R^2, '= 0.31')), adj=c(0,0))
text(x=0.42, y=4.6,labels=paste('p value < 0.001'), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ NDVI_2017mean_1m, data = soil_0_30cm_shp)) #R2=0.31

#Fig2d
tiff(file = file.path(FiguresDir, 'Fig2d.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$NDVI_2018mean_1m, soil_0_30cm_shp$kgOrgC.m2, xlab='Mean 2018 NDVI', ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ NDVI_2018mean_1m, data = soil_0_30cm_shp), lty=2)
text(x=0.52, y=5.2, labels=expression(paste(R^2, '< 0.01')), adj=c(0,0))
text(x=0.48, y=4.8,labels=paste('p value = 0.65'), adj=c(0,0))
text(x=0.36, y=5.4, labels='d')
dev.off()
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_0_30cm_shp)) #r^2<0.01 p-val=0.65
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_0_10cm_shp)) #r2 < 0.01 p-val=0.34
summary(lm(kgOrgC.m2 ~ NDVI_2018max_1m, data = soil_10_30cm_shp)) #r2 < 0.01 p-val=0.97

#Fig2e
tiff(file = file.path(FiguresDir, 'Fig2e.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$slope, soil_0_30cm_shp$kgOrgC.m2, xlab='Slope (%)', ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ slope, data = soil_0_30cm_shp), lty=2)
text(x=0, y=2.7, labels=expression(paste(R^2, '< 0.01')), adj=c(0,0))
text(x=0, y=2.0,labels=paste('p value = 0.61'), adj=c(0,0))
text(x=2, y=5.4, labels='e')
dev.off()
summary(lm(kgOrgC.m2 ~ slope, data = soil_0_30cm_shp)) #NS: p-val=0.61
summary(lm(kgOrgC.m2 ~ slope, data = soil_0_10cm_shp)) #r2=0.05
summary(lm(kgOrgC.m2 ~ slope, data = soil_10_30cm_shp)) #NS: p-val=0.16

#Fig2f
tiff(file = file.path(FiguresDir, 'Fig2f.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(soil_0_30cm_shp$annual_kwh.m2, soil_0_30cm_shp$kgOrgC.m2, xlab=expression(paste('Annual clear sky insolation (kWh', ' ', yr^-1, ')')), ylab=expression(paste('0-30 cm SOC (kg ', ~m^-2, ')')), pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
abline(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_0_30cm_shp), lty=2)
text(x=1100, y=2.5, labels=expression(paste(R^2, '= 0.03')), adj=c(0,0))
text(x=1100, y=2.1,labels=paste('p value = 0.10'), adj=c(0,0))
text(x=1095, y=5.4, labels='f')
dev.off()
summary(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_0_30cm_shp)) #NS; pval=0.10
summary(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_0_10cm_shp)) #NS
summary(lm(kgOrgC.m2 ~ annual_kwh.m2, data = soil_10_30cm_shp)) #r2=0.09

#Fig 3
all_forage_sp <- shapefile(file.path(spatialforageDir, 'all_pts_2017_2018.shp'))
lm_terrain5_0_30cm <- lm(kgOrgC.m2 ~ curvature_mean + slope + annual_kwh.m2 + elevation + NDVI_2017mean_1m, data =  soil_0_30cm_shp)
summary(lm_terrain5_0_30cm)
kgOrgC.m2_terrain5_0_30cm <- predict(Mar2017_terrain_3m_cropped, lm_terrain5_0_30cm)
#check relationship between SOC at forage sampling points and peak forage in both years
all_forage_sp$kgOrgC.m2_lm.terrain5 <- extract(kgOrgC.m2_terrain5_0_30cm, all_forage_sp)

#peak 2017 forage vs. modeled SOC
#Fig 5a
tiff(file = file.path(FiguresDir, 'Fig5a.tif', sep = ''), family = 'Times New Roman', width = 3.25, height = 3.3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(all_forage_sp$kgOrgC.m2_lm.terrain5, all_forage_sp$peak_2017, xlab=expression('Estimated SOC (kg'~m^-2*')'), ylab=expression('Peak forage 2017 (kg'~ha^-1*')'), xlim=c(2.3,4.8), xaxt='n')
axis(side=1, at=c(2.5, 3, 3.5, 4, 4.5), labels = as.character(c(2.5,3,3.5,4,4.5)))
abline(lm(peak_2017 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp), lty=2)
text(3.95, 1700, label=expression(paste(R^2, '= 0.23')), adj=c(0,0))
text(3.6, 1400, label='p value = 0.004', adj=c(0,0))
text(2.3, 1025, label=expression('slope = 1.1 kg forage'~kg^-1~'SOC'~m^-2), adj=c(0,0), cex=0.9)
text(x=2.5, y=4500, label='a', adj=c(0,0))
dev.off()
summary(lm(peak_2017 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp)) #R2=0.23; p-val=0.004

#peak 2018 vs SOC
#Fig5b
tiff(file = file.path(FiguresDir, 'Fig5b.tif', sep = ''), family = 'Times New Roman', width = 3.25, height = 3.3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(all_forage_sp$kgOrgC.m2_lm.terrain5, all_forage_sp$peak_2018, xlab=expression('Estimated SOC (kg'~m^-2*')'), ylab=expression('Peak forage 2018 (kg'~ha^-1*')'), xlim=c(2.3,4.8), xaxt='n', ylim=c(475,1500), cex.axis=0.95)
axis(side=1, at=c(2.5, 3, 3.5, 4, 4.5), labels = as.character(c(2.5,3,3.5,4,4.5)))
abline(lm(peak_2018 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp), lty=2)
text(2.3, 1100, label=expression(paste(R^2, '= 0.22')), adj=c(0,0))
text(2.3, 1010, label='p value = 0.06', adj=c(0.0))
text(2.3, 875, label=expression('slope = 0.37 kg forage'~kg^-1~'SOC'~m^-2), adj=c(0,0), cex=0.9)
text(2.5, 1400, label='b', adj=c(0.0))
dev.off()
summary(lm(peak_2018 ~ kgOrgC.m2_lm.terrain5, data = all_forage_sp)) #r2=0.22; p-val=0.06

#5var predicted vs. observed plot
#Fig 5a
tiff(file = file.path(FiguresDir, 'Fig5a.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values, soil_0_30cm_shp$kgOrgC.m2, xlab='', ylab='', pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
#points(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values[c(2, 5, 82)], soil_0_30cm_shp$kgOrgC.m2[c(2, 5, 82)], col='red', pch=19)
abline(0, 1, lty=2)
mtext(text=expression(paste('Predicted 0-30 cm SOC (kg', ~m^-2, ')')), side=1, line=2.75)
mtext(text=expression(paste('Observed 0-30 cm SOC (kg ', ~m^-2, ')')), side=2, line=2.75)
text(x=2, y=4.8, labels=expression(paste(R^2, '= 0.50, all data')), adj=c(0,0))
text(x=2, y=5.2, labels='b', adj=c(0,0))
#text(x=2, y=4.8,labels=expression(paste(r^2, '= 0.62, red pts removed')), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ curvature_mean + annual_kwh.m2 + slope + elevation + NDVI_2017mean_1m, data = soil_0_30cm_df)) #r2=0.50
summary(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)) #r2=0.41

#2var predicted vs. observed plot
#Fig5b
tiff(file = file.path(FiguresDir, 'Fig5b.tif', sep = ''), family = 'Times New Roman', width = 3, height = 3, pointsize = 11, units = 'in', res=res_plot)
par(mar=c(4.5, 4.5, 1, 1))
plot(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values, soil_0_30cm_shp$kgOrgC.m2, xlab='', ylab='', pch=21, cex.axis=1, cex.lab=1) #col=soil_0_30cm_shp$energy_colors, ylim = c(300, 1600), xlim=c(1400, 4700)
#points(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)$fitted.values[c(2, 5, 82)], soil_0_30cm_shp$kgOrgC.m2[c(2, 5, 82)], col='red', pch=19)
abline(0, 1, lty=2)
mtext(text=expression(paste('Predicted 0-30 cm SOC (kg', ~m^-2, ')')), side=1, line=2.75)
mtext(text=expression(paste('Observed 0-30 cm SOC (kg ', ~m^-2, ')')), side=2, line=2.75)
#text(x=2, y=5.2, labels=expression(paste(r^2, '= 0.50, all data')), adj=c(0,0))
text(x=2.25, y=5.2, labels='a', adj=c(0,0))
text(x=2.25, y=4.8, labels=expression(paste(R^2, '= 0.41, all data')), adj=c(0,0))
dev.off()
summary(lm(kgOrgC.m2 ~ curvature_mean + NDVI_2017mean_1m, data = soil_0_30cm_df)) #r2=0.41

#data for Fig6a
kmeans_cluster(2, c('NDVI_2017mean_1m_norm', 'curvature_mean_norm'), writetofile =  TRUE, fname = 'Fig6a.tif')

#data for Fig6b
kmeans_cluster(3, c('NDVI_2017mean_1m_norm', 'curvature_mean_norm'), writetofile = TRUE, 'Fig6b.tif')
