#TO-DO
#(1) re-check quality of organic C data
#(2)
mainDir <- 'C:/Users/smdevine/Desktop/rangeland project'
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
soilCDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data/soil C'
soilDataDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data'
results <- 'C:/Users/smdevine/Desktop/rangeland project/results'
list.files(soilCDir)
soilmasses <- read.csv(file.path(soilCDir, 'Camatta_soilCN_masses.csv'), stringsAsFactors = FALSE)
dim(soilmasses) #498
head(soilmasses)
soilCN <- read.csv(file.path(soilCDir, 'Camatta_soilCN_data.csv'), stringsAsFactors = FALSE)
head(soilCN)
tail(soilCN)
soilmasses$tray.ID <- paste0(soilmasses$tray, '-', soilmasses$well.plate)
dim(soilCN) #498 samples analyzed (15 extra)
soilCN <- merge(soilCN, soilmasses, by='tray.ID')
dim(soilCN) #498
length(unique(soilCN$sample.ID)) #213 unique
water_contents <- read.csv(file.path(soilCDir, 'airdry_soilmoisture.csv'), stringsAsFactors = FALSE)
soilCN <- merge(soilCN, water_contents, by='sample.ID')
dim(soilCN) #495 (dropped analyses with W in sample.ID, which were orgC that had not been pre-dried prior to weighing)
length(unique(soilCN$sample.ID)) #now 210 unique
soilCtotN <- soilCN[soilCN$analysis.type=='TC',]
soilCorgN <- soilCN[soilCN$analysis.type=='OC',]
dim(soilCtotN) #231, 21 dups
tapply(soilCtotN$C.mg, soilCtotN$sample.ID, length)#57_1 has two data points
soilCtotN$totC.percent <- 100 * soilCtotN$C.mg / soilCtotN$mass.mg
#soilCtotN$totC.percent <- round(100 * (soilCtotN$totC.percent / (100 - soilCtotN$H2O.perc.Adbasis)), 3)
hist(soilCtotN$totC.percent)
summary(soilCtotN$totC.percent)
soilCtotN$totN.percent <- 100 * soilCtotN$N.mg / soilCtotN$mass.mg
#soilCtotN$totN.percent <- round(100 * (soilCtotN$totN.percent / (100 - soilCtotN$H2O.perc.Adbasis)), 4)
soilCtotN[soilCtotN$sample.ID=='57_1',] #both ok
hist(soilCtotN$totN.percent)
summary(soilCtotN$totN.percent)
which(soilCtotN$totC.percent > 4)
soilCtotN[which(soilCtotN$totC.percent > 3), c('sample.ID', 'totC.percent', 'totN.percent')]
sum(duplicated(soilCtotN$sample.ID)) #21 duplicated
duplicated_samples <- which(soilCtotN$sample.ID %in% soilCtotN$sample.ID[duplicated(soilCtotN$sample.ID)])
soilCtotN[duplicated_samples, c('sample.ID', 'totC.percent', 'totN.percent', 'tray.ID')] #all were ok except 14_1
soilCtotN[soilCtotN$sample.ID=='14_1',] #1-B10 was erroneous, all others OK
soilCtotN <- soilCtotN[-which(soilCtotN$tray.ID=='1-B10'), ]
dim(soilCtotN) #230 now
dim(soilCorgN) #264 from draft dataset
length(unique(soilCorgN$sample.ID)) #210 unique as expected
soilCorgN$orgC.percent <- round(100 * soilCorgN$C.mg / soilCorgN$mass.mg, 3)
soilCorgN$totN.percent.v2 <- round(100 * soilCorgN$N.mg / soilCorgN$mass.mg, 4)
hist(soilCorgN$orgC.percent)
summary(soilCorgN$orgC.percent)
duplicated_samples_orgC <- which(soilCorgN$sample.ID %in% soilCorgN$sample.ID[duplicated(soilCorgN$sample.ID)])
soilCorgN[duplicated_samples_orgC, c('sample.ID', 'orgC.percent', 'totN.percent.v2', 'tray.ID')]
QC_orgC <- data.frame(sample.ID =as.character(unique(soilCorgN$sample.ID)), orgC.error = as.numeric(tapply(soilCorgN$orgC.percent, soilCorgN$sample.ID, function(x) round(100*(max(x) - min(x))/mean(x), 1))))
QC_orgC[QC_orgC$orgC.error > 15,]
samples_to_review <- QC_orgC[QC_orgC$orgC.error > 15,]
samples_to_review[order(samples_to_review$orgC.error),]
check <- soilCorgN[soilCorgN$sample.ID %in% samples_to_review$sample.ID, c('sample.ID', 'orgC.percent', 'totN.percent.v2', 'tray.ID')]
check[order(check$sample.ID, check$tray.ID),]
#soilCorgN$sample.ID
colnames(soilCtotN)
colnames(soilCtotN)[2] <- 'tray.ID.TC'
colnames(soilCorgN)
colnames(soilCorgN)[2] <- 'tray.ID.OC'
soilC_all <- merge(soilCtotN[,c('sample.ID', 'tray.ID.TC', 'totC.percent', 'totN.percent')], soilCorgN[ ,c('sample.ID', 'tray.ID.OC', 'orgC.percent', 'totN.percent.v2')], by='sample.ID') #W in org C samples denoted those that had not been dried at 60C first before immediately weighing !grepl('W', soilCorgN$sample.ID)
dim(soilC_all) #299 rows
length(unique(soilC_all$sample.ID)) #210 unique, so 11 duplicates (2 OC, 9 TC as of 11/8/18)
soilC_all$inorgC.percent <- soilC_all$totC.percent - soilC_all$orgC.percent
soilC_all$orgC.to.totC <- round(soilC_all$orgC.percent / soilC_all$totC.percent, 3)
summary(soilC_all$inorgC.percent)
hist(soilC_all$orgC.to.totC)
hist(soilC_all$orgC.to.totC[grepl('_1', soilC_all$sample.ID)])
hist(soilC_all$orgC.to.totC[grepl('_2', soilC_all$sample.ID)])
sum(soilC_all$inorgC.percent < 0) #11 less than 0
soilC_all[soilC_all$inorgC.percent < 0,]
plot(soilC_all$totN.percent, soilC_all$totN.percent.v2)
plot(soilC_all$totN.percent, soilC_all$totN.percent.v2, col=ifelse(soilC_all$sample.ID=='65_1', 'red', 'black'))
soilC_all[soilC_all$sample.ID=='65_1',]
samples_to_review[order(samples_to_review$orgC.error),]
plot(soilC_all$totN.percent, soilC_all$totN.percent.v2, col=ifelse(soilC_all$sample.ID=='65_1', 'red', 'black'))
summary(lm(totN.percent.v2 ~ totN.percent, data = soilC_all))
abline(0, 1, lty=2)
#abline(lm(totN.percent.v2 ~ totN.percent, data = soilC_all), lty=2)
text(soilC_all$totN.percent, soilC_all$totN.percent.v2, labels = soilC_all$sample.ID, pos=1, offset = 0.3)
soilC_all$N.error.percent <- round(100 * (soilC_all$totN.percent.v2 - soilC_all$totN.percent) / soilC_all$totN.percent, 2)
summary(soilC_all$N.error.percent)
hist(soilC_all$N.error.percent)
sum(abs(soilC_all$N.error.percent) > 20) #9 greater than 20% error
sum(abs(soilC_all$N.error.percent) > 15) #27 greater than 20% error
sum(abs(soilC_all$N.error.percent) > 10) #59 greater than 10% error

soilC_all$N.diff.abs <- abs(soilC_all$totN.percent - soilC_all$totN.percent.v2)
summary(soilC_all$N.diff.abs)
hist(soilC_all$N.diff.abs)
soilC_all$N.diff <- soilC_all$totN.percent - soilC_all$totN.percent.v2
summary(soilC_all$N.diff)
hist(soilC_all$N.diff.abs[soilC_all$N.diff.abs < 0.1])
soilC_all[order(abs(soilC_all$N.error.percent)), c('sample.ID', 'N.error.percent', 'N.diff')]
mean(soilC_all$totN.percent) #mean is 0.115 %N
sd(soilC_all$totN.percent) #sd is 0.05
sum(soilC_all$N.diff.abs > 0.01) #72 greater than 0.01% difference
sum(soilC_all$N.diff.abs > 0.015) #35 greater
sum(soilC_all$N.diff.abs > 0.02) #22 greater than this
#N_QC_summary <- summary(lm(totN.percent.v2 ~ totN.percent, data = soilC_all))
#N_QC_summary$residuals[order(N_QC_summary$residuals)]
soilC_all$CaCO3.percent <- round((100/12) * soilC_all$inorgC.percent, 3)
summary(soilC_all$CaCO3.percent)
soilC_all[order(soilC_all$CaCO3.percent), ] #two are 45_2
soilC_all[soilC_all$sample.ID=='45_2', ]
soilCorgN[soilCorgN$sample.ID=='45_2', ]
soilC_reruns <- soilC_all[abs(soilC_all$N.error.percent) > 15,]
dim(soilC_reruns)
soilC_reruns$N.diff.abs <- abs(soilC_reruns$totN.percent - soilC_reruns$totN.percent.v2)
soilC_reruns[order(abs(soilC_reruns$N.error.percent)),]
sum(soilC_reruns$N.diff.abs > 0.02) #11 are greater
sum(soilC_reruns$N.diff.abs > 0.015) #18 are greater
plot(soilC_all$orgC.percent, soilC_all$N.diff.abs)
soilC_reruns.v2 <- soilC_all[soilC_all$N.diff.abs > 0.01 & abs(soilC_all$N.error.percent) > 10,] #if more than 15% relative difference with absolute difference at least 0.010 %N
dim(soilC_reruns.v2) #48 re-runs
soilC_reruns.v2[abs(soilC_reruns.v2$N.error.percent) > 20, ]
sum(abs(soilC_reruns.v2$N.error.percent) > 20) #8 have relative difference of 20% or more
sum(abs(soilC_reruns.v2$N.error.percent) > 15) #25 have relative difference of 15% or more
soilC_reruns.v2
soilC_reruns.v2$depth <- ifelse(grepl('_2', soilC_reruns.v2$sample.ID), 2, 1)
soilC_reruns.v2$location <- gsub('_2', '', soilC_reruns.v2$sample.ID)
soilC_reruns.v2$location <- gsub('_1', '', soilC_reruns.v2$location)
soilC_reruns.v2$location <- as.integer(soilC_reruns.v2$location)
soilC_reruns.v2 <- soilC_reruns.v2[order(soilC_reruns.v2$location),]
#write.csv(soilC_reruns.v2, file.path(soilCDir, 'soilC_reruns.v2.csv'), row.names = FALSE)
summary(soilC_reruns.v2$CaCO3.percent)
summary(soilC_all$CaCO3.percent)
colnames(soilC_all)
soilC_0_10cm <- soilC_all[grepl('_1', soilC_all$sample.ID), c("sample.ID", "totC.percent", "totN.percent", "orgC.percent", "totN.percent.v2", "inorgC.percent", "orgC.to.totC", "CaCO3.percent", 'N.error.percent')]
soilC_0_10cm[which(soilC_0_10cm$sample.ID %in% soilC_0_10cm$sample.ID[duplicated(soilC_0_10cm$sample.ID)]),]
soilC_0_10cm <- data.frame(sample.ID=unique(soilC_0_10cm$sample.ID), totC.percent=as.numeric(tapply(soilC_0_10cm$totC.percent, soilC_0_10cm$sample.ID, mean)), totN.percent=as.numeric(tapply(soilC_0_10cm$totN.percent, soilC_0_10cm$sample.ID, mean)), orgC.percent=as.numeric(tapply(soilC_0_10cm$orgC.percent, soilC_0_10cm$sample.ID, mean)), inorgC.percent=as.numeric(tapply(soilC_0_10cm$inorgC.percent, soilC_0_10cm$sample.ID, mean)), CaCO3.percent=as.numeric(tapply(soilC_0_10cm$CaCO3.percent, soilC_0_10cm$sample.ID, mean)), orgC.to.totC=as.numeric(tapply(soilC_0_10cm$orgC.to.totC, soilC_0_10cm$sample.ID, mean)), N.error.percent=as.numeric(tapply(soilC_0_10cm$N.error.percent, soilC_0_10cm$sample.ID, mean)))
soilC_0_10cm$inorgC.percent[soilC_0_10cm$inorgC.percent < 0] <- 0 #4 were negative, meaning no detectable carbonates
soilC_0_10cm$CaCO3.percent[soilC_0_10cm$CaCO3.percent < 0] <- 0
soilC_10_30cm <- soilC_all[grepl('_2', soilC_all$sample.ID), ]
soilC_10_30cm <- data.frame(sample.ID=unique(soilC_10_30cm$sample.ID), totC.percent=as.numeric(tapply(soilC_10_30cm$totC.percent, soilC_10_30cm$sample.ID, mean)), totN.percent=as.numeric(tapply(soilC_10_30cm$totN.percent, soilC_10_30cm$sample.ID, mean)), orgC.percent=as.numeric(tapply(soilC_10_30cm$orgC.percent, soilC_10_30cm$sample.ID, mean)), inorgC.percent=as.numeric(tapply(soilC_10_30cm$inorgC.percent, soilC_10_30cm$sample.ID, mean)), CaCO3.percent=as.numeric(tapply(soilC_10_30cm$CaCO3.percent, soilC_10_30cm$sample.ID, mean)), orgC.to.totC=as.numeric(tapply(soilC_10_30cm$orgC.to.totC, soilC_10_30cm$sample.ID, mean)), N.error.percent=as.numeric(tapply(soilC_10_30cm$N.error.percent, soilC_10_30cm$sample.ID, mean)))
soilC_10_30cm$inorgC.percent[soilC_10_30cm$inorgC.percent < 0] <- 0 #2 samples had undetectable levels of carbonates
soilC_10_30cm$CaCO3.percent[soilC_10_30cm$CaCO3.percent < 0] <- 0
mean(soilC_0_10cm$totN.percent)
sd(soilC_0_10cm$totN.percent)
mean(soilC_0_10cm$N.diff.abs)
sum(abs(soilC_0_10cm$N.error.percent) > 10) #22 greater than 10% rel. diff
sum(abs(soilC_0_10cm$N.error.percent) > 5) #56 greater than 5% rel. diff
mean(soilC_10_30cm$totN.percent)
sd(soilC_10_30cm$totN.percent)
mean(soilC_10_30cm$N.diff.abs)
sum(abs(soilC_10_30cm$N.error.percent) > 10) #32 greater than 10% rel. diff
hist(soilC_0_10cm$orgC.percent)
dim(soilC_0_10cm)
dim(soilC_10_30cm)
length(unique(soilC_0_10cm$sample.ID))
tapply(soilC_0_10cm$totC.percent, soilC_0_10cm$sample.ID, length) #57_1 and 58_1 are reps
# rownames(soilC_0_10cm) <- 1:nrow(soilC_0_10cm)
# soilC_0_10cm[soilC_0_10cm$sample.ID=='57_1',]
# soilC_0_10cm[soilC_0_10cm$sample.ID=='58_1',]
#temp fix to get rid of duplicate IDs
# soilC_0_10cm <- soilC_0_10cm[-c(36, 39), ]
hist(soilC_10_30cm$orgC.percent)
hist(soilC_10_30cm$orgC.percent)
#write.csv(soilC_0_10cm, file.path(soilCDir, 'soilC_0_10cm.draft.csv'), row.names = FALSE)
#write.csv(soilC_10_30cm, file.path(soilCDir, 'soilC_10_30cm.draft.csv'), row.names = FALSE)

#read in bulk density data
#read-in BD_data
BD_data <- read.csv(file.path(results, 'soil_data', 'soilBD_H2O_frags_2018-09-17.csv'), stringsAsFactors = FALSE)
#split dataset by depth
BD_data_0_10cm <- BD_data[BD_data$depth_code==1 | BD_data$depth_code==3, ]
BD_data_10_30cm <- BD_data[BD_data$depth_code==2,]

#read-in terrain characteristics from S Hogan's drone data
library(randomForest)
library(raster)
#soils characterization data
soil_chars <- read.csv(file.path(soilDataDir, 'A&L.Western.Labs.Results.csv'), stringsAsFactors = FALSE)
soil_chars_0_10 <- soil_chars[grepl('_1', soil_chars$sample.ID),]
soil_chars_0_10$point_no <-  as.integer(gsub('_1', '', soil_chars_0_10$sample.ID))
soil_chars_10_30 <- soil_chars[grepl('_2', soil_chars$sample.ID),]
soil_chars_10_30$point_no <-  as.integer(gsub('_2', '', soil_chars_10_30$sample.ID))
#merge with spatial dataset
sampling_pts <- shapefile(file.path(mainDir, 'sampling points 2018', 'soil_sampling_points.shp'))
#biomass_Apr2017 <- raster(file.path(mainDir, 'sampling strategy April 2018', 'Biomass_2017-04-10_APAR.tif'))
sampling_pts$point_no <- as.integer(gsub('point', '', sampling_pts$Comment))
list.files(file.path(terrainDir, 'filtered_Hogan'))
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan'), full.names = TRUE))
#Mar2017_terrain_3m
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
#read in 3 m solrad raster produced in ArcGIS
list.files(file.path(terrainDir, 'solrad_analysis')) #all files were sky size 500 x 500 and 64 calc directions
 #all v2 files were sky size 500 x 500 and 64 calc directions
solrad <- shapefile(file.path(terrainDir, 'solrad_analysis', 'solrad_105pts.shp'))
names(solrad)
#solrad$point_no <- 1:105
#solrad$annual_kwh.m2 <- apply(as.data.frame(solrad[ ,1:52]), 1, sum) / 1000
solrad_df <- as.data.frame(solrad)
colnames(solrad_df)
head(solrad_df)
solrad_df$point_no <- sampling_pts$point_no
solrad_df$annual_kwh.m2 <- apply(solrad_df[ ,1:52], 1, sum) / 1000
summary(solrad_df$annual_kwh.m2)
length(seq.Date(as.Date('Oct_01_2017', '%b_%d_%Y'), as.Date('Apr_15_2018', '%b_%d_%Y'), by='day')) #197 days
format(as.Date('Oct_01_2017', '%b_%d_%Y'), '%V') #so, T38-T51
format(as.Date('Apr_15_2018', '%b_%d_%Y'), '%V') #and T0-T15
solrad_df$seasonal_kwh.m2 <- apply(solrad_df[ ,c(1:16, 39:52)], 1, sum) / 1000
sampling_pts <- merge(sampling_pts, solrad_df[ ,c('annual_kwh.m2', 'point_no')], by='point_no')
head(sampling_pts)
sampling_pts$slope <- extract(Mar2017_terrain_3m$slope, coordinates(sampling_pts)[,1:2])
sampling_pts$curvature_mean <- extract(Mar2017_terrain_3m$curvature_mean, coordinates(sampling_pts)[,1:2])
sampling_pts$elevation <- extract(Mar2017_terrain_3m$elevation, coordinates(sampling_pts)[,1:2])
sampling_pts$TCI <- extract(Mar2017_terrain_3m$TCI, coordinates(sampling_pts)[,1:2])
sampling_pts$curvature_plan <- extract(Mar2017_terrain_3m$curvature_plan, coordinates(sampling_pts)[,1:2])
sampling_pts$curvature_profile <- extract(Mar2017_terrain_3m$curvature_profile, coordinates(sampling_pts)[,1:2])
sampling_pts$aspect <- extract(Mar2017_terrain_3m$aspect, coordinates(sampling_pts)[,1:2])
soilC_0_10cm$point_no <- as.integer(gsub('_1', '', soilC_0_10cm$sample.ID))
soilC_10_30cm$point_no <- as.integer(gsub('_2', '', soilC_10_30cm$sample.ID))
soil_0_10cm_shp <- merge(sampling_pts, soilC_0_10cm, by='point_no')
soil_0_10cm_shp <- merge(soil_0_10cm_shp, soil_chars_0_10, by='point_no')
soil_0_10cm_shp <- merge(soil_0_10cm_shp, BD_data_0_10cm, by='point_no')
soil_0_10cm_shp$soil.kg.orgC.m2 <- soil_0_10cm_shp$orgC.percent * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100
soil_0_10cm_shp$soil.kg.IC.m2 <- soil_0_10cm_shp$s
soil_0_10cm_shp$soil.kgClay.m2 <- soil_0_10cm_shp$CLAY * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100
hist(soil_0_10cm_shp$soil.kg.orgC.m2)
plot(1:105, soil_0_10cm_shp$soil.kg.orgC.m2)
text(1:105, soil_0_10cm_shp$soil.kg.orgC.m2, labels=soil_0_10cm_shp$point_no, offset=0.2, pos=1)

soil_10_30cm_shp <- merge(sampling_pts, soilC_10_30cm, by='point_no')
soil_10_30cm_shp <- merge(soil_10_30cm_shp, soil_chars_10_30, by='point_no')
soil_10_30cm_shp <- merge(soil_10_30cm_shp, BD_data_10_30cm, by='point_no')
#soil_0_10cm_shp$energy_colors <- 
soil_10_30cm_shp$soil.kg.orgC.m2 <- soil_10_30cm_shp$orgC.percent * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100
soil_10_30cm_shp$soil.kgClay.m2 <- soil_10_30cm_shp$CLAY * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100
hist(soil_10_30cm_shp$soil.kg.orgC.m2)
hist(soil_0_10cm_shp$soil.kg.orgC.m2 + soil_10_30cm_shp$soil.kg.orgC.m2)
soil_0_10cm_shp$point_no - soil_10_30cm_shp$point_no
sampling_pts$sand_wtd <- (10*soil_0_10cm_shp$SAND + 20*soil_10_30cm_shp$SAND) / 30
sampling_pts$clay_wtd <- (10*soil_0_10cm_shp$CLAY + 20*soil_10_30cm_shp$CLAY) / 30
sampling_pts$clay_content <- soil_0_10cm_shp$soil.kgClay.m2 + soil_10_30cm_shp$soil.kgClay.m2
#sampling_pts$sand_content <-
#sampling_pts$Apr2017biomass <- extract(biomass_Apr2017, coordinates(sampling_pts)[,1:2], fun=mean, buffer=1)
hist(sampling_pts$sand_wtd)
hist(sampling_pts$clay_wtd)
sampling_pts$orgC_content <- soil_0_10cm_shp$soil.kg.orgC.m2 + soil_10_30cm_shp$soil.kg.orgC.m2
#write.csv(as.data.frame(sampling_pts), 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/terrain.features/terrain_105pts.csv', row.names=FALSE)

#plot soil C as interpolated map
#see labs 14 and 15 from Quant Geo for tips
#function interpolate is good start
library(gstat)
#r <- raster(extent(Mar2017_terrain_3m), resolution=res(Mar2017_terrain_3m), crs=crs(Mar2017_terrain_3m))
r <- raster(extent(soil_0_10cm_shp), resolution=3, crs=crs(soil_0_10cm_shp))
#values(r) <- 1:ncell(r)
#plot(r)
#plot(soil_0_10cm_shp, add=TRUE)
class(soil_0_10cm_shp$CLAY)
summary(soil_0_10cm_shp$CLAY)
test <- soil_0_10cm_shp[,'CLAY']
summary(soil_0_10cm_shp$CLAY)
names(sampling_pts)
idm <- gstat(id='orgC_content', formula = orgC_content~1, data = sampling_pts, nmax=10)
idp <- interpolate(r, idm)
summary(idp)
plot(idp)
summary(sampling_pts$orgC_content)
#optimize a couple of parameters in interpolate
#from lab 14 line 112
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula=OZDLYAV~1, locations=train, nmax=nmx, set=list(idp=idp))
  p <- predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(test$OZDLYAV, p)
}
set.seed(20150518)
i <- sample(nrow(soil_0_10cm_shp), 0.2 * nrow(soil_0_10cm_shp))
tst <- soil_0_10cm_shp[i,]
trn <- aq[-i,]
opt <- optim(c(8, .5), f1, test=tst, train=trn)
opt
#now use the optimal IDW model
m <- gstat(formula=OZDLYAV~1, locations=aq, nmax=opt$par[1], set=list(idp=opt$par[2]))
idw <- interpolate(r, m)
idw <- mask(idw, ca)
plot(idw)

gs <- gstat(formula=orgC_content ~ 1, locations=sampling_pts)
v <- variogram(gs, width=10, cutoff=150)
v
plot(v)
fve <- fit.variogram(v, vgm(model = "Exp"))
fve
plot(variogramLine(fve, 400), type='l')
points(v[,2:3], pch=20, col='red')
#Try a different type (spherical in stead of exponential)
fvs <- fit.variogram(v, vgm(model="Sph"))
fvs
plot(variogramLine(fvs, 200), type='l', col='blue', lwd=2)
points(v[,2:3], pch=20, col='red')
#Another way to plot the variogram and the model
plot(v, fve)
plot(v, fvs)
#Use variogram fve in a kriging interpolation
k <- gstat(formula=orgC_content ~ 1, locations=sampling_pts, model=fve)
#gstat(formula=orgC_content ~ 1, locations=sampling_pts)
class(k)
# predicted values

g <- as(r, 'SpatialGrid')
kp <- predict(k, g)
class(kp)
spplot(kp)
kp_raster <- raster(kp)
plot(kp_raster)


plot(biomass_Apr2017)
plot(sampling_pts, cex=sampling_pts$orgC_content/2, pch=1, add=TRUE)
summary(lm(orgC_content ~ elevation + slope + curvature_mean + annual_kwh.m2, data= sampling_pts))
summary(lm(orgC_content ~ elevation + slope + curvature_mean + annual_kwh.m2 + sand_wtd, data= sampling_pts))
summary(lm(orgC_content ~ slope + curvature_mean + annual_kwh.m2 + clay_wtd, data= sampling_pts))
plot(lm(orgC_content ~ slope + curvature_mean + annual_kwh.m2 + clay_wtd, data= sampling_pts))
sampling_pts$C_residuals <- lm(orgC_content ~ slope + curvature_mean + annual_kwh.m2 + clay_wtd, data= sampling_pts)$residuals
plot(sampling_pts, cex=abs(sampling_pts$C_residuals))
summary(lm(orgC_content ~ slope + curvature_mean + annual_kwh.m2 + clay_content, data= sampling_pts))

#summary(lm(Apr2017biomass ~ slope + curvature_mean + annual_kwh.m2 + clay_wtd, data= sampling_pts))
#summary(lm(Apr2017biomass ~ orgC_content + clay_wtd, data= sampling_pts))
#summary(lm(Apr2017biomass ~ orgC_content + sand_wtd, data= sampling_pts))
#summary(lm(Apr2017biomass ~ orgC_content + sand_wtd + annual_kwh.m2, data= sampling_pts))
#summary(lm(Apr2017biomass ~ orgC_content + sand_wtd + annual_kwh.m2 + curvature_mean, data= sampling_pts)) #r^2=0.48
#summary(lm(Apr2017biomass ~ orgC_content + sand_wtd + annual_kwh.m2 + curvature_mean + slope, data= sampling_pts)) 
#summary(lm(Apr2017biomass ~ orgC_content + sand_wtd + annual_kwh.m2 + curvature_mean + slope + elevation, data= sampling_pts))#r^2 =0.67
#summary(lm(Apr2017biomass ~ annual_kwh.m2 + curvature_mean + slope + elevation, data= sampling_pts)) #r^2=0.67
#summary(lm(Apr2017biomass ~ curvature_mean + slope + elevation, data= sampling_pts)) #r2=0.66
#summary(lm(Apr2017biomass ~ curvature_mean, data= sampling_pts)) #r^2=0.22
#summary(lm(Apr2017biomass ~ slope, data= sampling_pts)) #r^2=0.08
#summary(lm(Apr2017biomass ~ elevation, data= sampling_pts)) #r^2= 0.36 strongest relationship
#summary(lm(Apr2017biomass ~ curvature_mean + elevation, data= sampling_pts)) #r2=0.41
#vif(lm(Apr2017biomass ~ orgC_content + sand_wtd + annual_kwh.m2 + curvature_mean + slope + elevation, data= sampling_pts))
#summary(lm(Apr2017biomass ~ slope + curvature_mean + annual_kwh.m2 + CLAY, data= soilC_10_30cm_df))
#summary(lm(Apr2017biomass ~ slope + curvature_mean + annual_kwh.m2, data= soilC_10_30cm_df))
#summary(lm(Apr2017biomass ~ slope + curvature_mean + annual_kwh.m2, data= soilC_0_10cm_df)) #double-check nothing funky
#summary(lm(Apr2017biomass ~ slope + curvature_mean + annual_kwh.m2 + SAND, data= soilC_0_10cm_df))
#summary(lm(Apr2017biomass ~ slope + curvature_mean + annual_kwh.m2 + SAND + orgC.percent, data= soilC_0_10cm_df))
#summary(lm(Apr2017biomass ~ orgC.percent + SAND, data= soilC_0_10cm_df))
plot(orgC_contents, soil_0_10cm_shp$Apr2017biomass, color)
plot(biomass_Apr2017)
plot(soil_0_10cm_shp, cex=soil_0_10cm_shp$orgC.percent, pch=1, add=T)
plot(soil_0_10cm_shp, cex=soil_0_10cm_shp$CLAY/20, pch=1)
soil_0_10cm_shp$Apr2017biomass <- extract(biomass_Apr2017, coordinates(soil_0_10cm_shp)[,1:2], buffer=1, fun=mean)
soil_10_30cm_shp$Apr2017biomass <- extract(biomass_Apr2017, coordinates(soil_10_30cm_shp)[,1:2], buffer=1, fun=mean)
plot(soil_0_10cm_shp$orgC.percent, soil_0_10cm_shp$Apr2017biomass, col=ifelse(soil_0_10cm_shp$N.error.percent > 10, 'red', 'black'))
plot(soil_0_10cm_shp$slope, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$annual_kwh.m2, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$curvature_mean, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$CLAY, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$SAND, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$SILT, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$SAND, soil_0_10cm_shp$Apr2017biomass)
plot(soil_0_10cm_shp$CLAY, soil_0_10cm_shp$Apr2017biomass)
plot(soil_0_10cm_shp$P1, soil_0_10cm_shp$orgC.percent)
plot(soil_0_10cm_shp$HCO3_P, soil_0_10cm_shp$orgC.percent)
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
summary(lm(orgC.percent ~ SAND + curvature_mean + slope + annual_kwh.m2, data = as.data.frame(soil_0_10cm_shp)))
lm_0_10cm <- lm(orgC.percent ~ SAND + curvature_mean + slope + annual_kwh.m2, data = as.data.frame(soil_0_10cm_shp))
plot(lm_0_10cm$fitted.values, soilC_10_30cm$orgC.percent)
plot(lm(orgC.percent ~ SAND + curvature_mean + slope + annual_kwh.m2, data = as.data.frame(soil_0_10cm_shp))) #59, 79, and 83 are problematic
summary(lm(orgC.percent ~ slope + annual_kwh.m2 + curvature_mean, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(orgC.percent ~ curvature_mean + slope, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(orgC.percent ~ curvature_profile + slope, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(orgC.percent ~ curvature_mean + slope + elevation, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(orgC.percent ~ curvature_mean + slope + elevation + Apr2017biomass, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(orgC.percent ~ SAND + curvature_mean + slope + annual_kwh.m2, data=as.data.frame(soil_10_30cm_shp))) #all params significant and 
lm_10_30cm <- lm(orgC.percent ~ SAND + curvature_mean + slope + annual_kwh.m2, data=as.data.frame(soil_10_30cm_shp))
plot(lm_10_30cm$fitted.values, soilC_10_30cm_df$orgC.percent)
plot(lm(orgC.percent ~ SAND + curvature_mean + slope + annual_kwh.m2, data=as.data.frame(soil_10_30cm_shp))) #5, 64, and 101 are problematic

plot(soil_0_10cm_shp$orgC.percent, soil_0_10cm_shp$Apr2017biomass)
text(soil_0_10cm_shp$orgC.percent, soil_0_10cm_shp$Apr2017biomass, labels = soil_0_10cm_shp$point_no, pos=1, offset =0.3)
abline(lm(Apr2017biomass ~ orgC.percent, data = soil_0_10cm_shp), lty=2)
summary(lm(Apr2017biomass ~ orgC.percent, data = soil_0_10cm_shp))
lm.result <- lm(Apr2017biomass ~ orgC.percent, data = soil_0_10cm_shp)
lm.result$residuals[order(lm.result$residuals)]
plot(lm.result)
summary(lm(Apr2017biomass ~ inorgC.percent, data = soil_0_10cm_shp))
summary(lm(Apr2017biomass ~ orgC.percent, data = soil_10_30cm_shp))
summary(lm(Apr2017biomass ~ inorgC.percent, data = soil_10_30cm_shp))
summary(lm(Apr2017biomass ~ orgC.percent + inorgC.percent, data = soil_10_30cm_shp))
plot(soil_10_30cm_shp$orgC.percent, soil_10_30cm_shp$Apr2017biomass, col=ifelse(soil_10_30cm_shp$N.error.percent > 5, 'red', 'black'))
abline(lm(Apr2017biomass ~ orgC.percent, data = soil_10_30cm_shp))
text(soil_10_30cm_shp$orgC.percent, soil_10_30cm_shp$Apr2017biomass, labels=soil_10_30cm_shp$point_no, pos=1, offset=0.3, cex=0.8)
plot(biomass_Apr2017)
plot(soil_10_30cm_shp, cex=soil_10_30cm_shp$orgC.percent*2, pch=1, add=T)
plot(lm(Apr2017biomass ~ orgC.percent, data = soil_10_30cm_shp))
plot(soil_0_10cm_shp$orgC.percent, soil_10_30cm_shp$orgC.percent)
abline(glm(soil_10_30cm_shp$orgC.percent ~ soil_0_10cm_shp$orgC.percent))
summary(lm(soil_10_30cm_shp$orgC.percent ~ soil_0_10cm_shp$orgC.percent))
abline(0.5179, 0.1365, col='red')
plot(lm(soil_10_30cm_shp$orgC.percent ~ soil_0_10cm_shp$orgC.percent))
