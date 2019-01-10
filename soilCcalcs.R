#TO-DO
#(1) re-check quality of organic C data [DONE]
#(2) re-aggregate data because of TN, SIC, and P content calc mistakes.  Reran from line 184 on 1/7/19 [DONE]
#fix points 2 and 3 0-10 cm data so that they are corrected for the fact they were just 0-5 cm samples
mainDir <- 'C:/Users/smdevine/Desktop/rangeland project'
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
solradDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/solrad_analysis'
soilCDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data/soil C'
soilDataDir <- 'C:/Users/smdevine/Desktop/rangeland project/soils_data'
soilCresults <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/soilCresults'
results <- 'C:/Users/smdevine/Desktop/rangeland project/results'
forageDir <- 'C:/Users/smdevine/Desktop/rangeland project/clip_plots/results'
list.files(forageDir)
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
hist(soilCtotN$totN.percent)
summary(soilCtotN$totN.percent)
which(soilCtotN$totC.percent > 4)
soilCtotN[which(soilCtotN$totC.percent > 3), c('sample.ID', 'totC.percent', 'totN.percent')]
sum(duplicated(soilCtotN$sample.ID)) #21 duplicated
duplicated_samples <- which(soilCtotN$sample.ID %in% soilCtotN$sample.ID[duplicated(soilCtotN$sample.ID)])
soilCtotN[duplicated_samples, c('sample.ID', 'totC.percent', 'totN.percent', 'tray.ID')] #all were ok except 14_1
soilCtotN[soilCtotN$sample.ID=='72_1',] #1-B10, 1-C6 were erroneous, all others OK; verified 12/1/18
soilCtotN <- soilCtotN[-which(soilCtotN$tray.ID=='1-B10'), ]
soilCtotN <- soilCtotN[-which(soilCtotN$tray.ID=='1-C6'), ]
dim(soilCtotN) #230 now
dim(soilCorgN) #264 from draft dataset
length(unique(soilCorgN$sample.ID)) #210 unique as expected
soilCorgN$orgC.percent <- round(100 * soilCorgN$C.mg / soilCorgN$mass.mg, 3)
soilCorgN$totN.percent.v2 <- round(100 * soilCorgN$N.mg / soilCorgN$mass.mg, 4)
hist(soilCorgN$orgC.percent)
summary(soilCorgN$orgC.percent)
soilCorgN[soilCorgN$sample.ID=='72_1',]
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

soilCorgN <- soilCorgN[-which(soilCorgN$tray.ID.OC=='6-D7'),] #total N way off compared to total C analysis and had written a note that was not positive that re-analysis label was 5-2
soilCorgN <- soilCorgN[-which(soilCorgN$tray.ID.OC=='6-D5'),]
soilCorgN <- soilCorgN[-which(soilCorgN$tray.ID.OC=='4-D1'),]
soilCorgN$orgC.percent[soilCorgN$sample.ID=='79_1'] <- NA #concluded that this analysis was erroneous; will base organic C estimate as avg% of total C estimate from adjancent points which all had 46-77% of total C
soilCorgN$totN.percent.v2[soilCorgN$sample.ID=='79_1'] <- NA
soilC_all <- merge(soilCtotN[,c('sample.ID', 'tray.ID.TC', 'totC.percent', 'totN.percent')], soilCorgN[ ,c('sample.ID', 'tray.ID.OC', 'orgC.percent', 'totN.percent.v2')], by='sample.ID') #W in org C samples denoted those that had not been dried at 60C first before immediately weighing !grepl('W', soilCorgN$sample.ID)
dim(soilC_all) #294 rows
length(unique(soilC_all$sample.ID)) #210 unique as expected
soilC_all$orgC.percent[soilC_all$sample.ID=='79_1']
soilC_all$inorgC.percent <- soilC_all$totC.percent - soilC_all$orgC.percent
soilC_all$orgC.to.totC <- round(soilC_all$orgC.percent / soilC_all$totC.percent, 3)
summary(soilC_all$inorgC.percent)
hist(soilC_all$orgC.to.totC)
hist(soilC_all$orgC.to.totC[grepl('_1', soilC_all$sample.ID)])
hist(soilC_all$orgC.to.totC[grepl('_2', soilC_all$sample.ID)])
sum(soilC_all$inorgC.percent < 0) #10 less than 0
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
sum(soilC_all$N.diff.abs > 0.01) #113 greater than 0.01% difference
sum(soilC_all$N.diff.abs > 0.015) #59 greater
sum(soilC_all$N.diff.abs > 0.02) #37 greater than this
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
sum(soilC_0_10cm$inorgC.percent < 0, na.rm=TRUE) #4 were negative, meaning no detectable carbonates
soilC_0_10cm$inorgC.percent[soilC_0_10cm$inorgC.percent < 0] <- 0 
soilC_0_10cm$CaCO3.percent[soilC_0_10cm$CaCO3.percent < 0] <- 0
soilC_0_10cm$orgC.to.totC[soilC_0_10cm$orgC.to.totC > 1] <- 1
#manual corrections for 79-1
soilC_0_10cm[soilC_0_10cm$sample.ID=='79_1',]
soilC_0_10cm$orgC.to.totC[soilC_0_10cm$sample.ID=='79_1'] <- round(mean(c(0.543, 0.592, 0.566, 0.462, 0.590, 0.766, 0.696)), 3)
soilC_0_10cm$orgC.percent[soilC_0_10cm$sample.ID=='79_1'] <- soilC_0_10cm$totC.percent[soilC_0_10cm$sample.ID=='79_1'] * mean(c(0.543, 0.592, 0.566, 0.462, 0.590, 0.766, 0.696))
soilC_0_10cm$inorgC.percent[soilC_0_10cm$sample.ID=='79_1'] <- soilC_0_10cm$totC.percent[soilC_0_10cm$sample.ID=='79_1'] * (1 - mean(c(0.543, 0.592, 0.566, 0.462, 0.590, 0.766, 0.696)))
soilC_0_10cm$CaCO3.percent[soilC_0_10cm$sample.ID=='79_1'] <- round((100/12) * soilC_0_10cm$inorgC.percent[soilC_0_10cm$sample.ID=='79_1'], 3)
soilC_0_10cm$orgC.percent[soilC_0_10cm$orgC.percent > soilC_0_10cm$totC.percent] <- soilC_0_10cm$totC.percent[soilC_0_10cm$orgC.percent > soilC_0_10cm$totC.percent]
#now soilC_10_30cm
soilC_10_30cm <- soilC_all[grepl('_2', soilC_all$sample.ID), ]
soilC_10_30cm <- data.frame(sample.ID=unique(soilC_10_30cm$sample.ID), totC.percent=as.numeric(tapply(soilC_10_30cm$totC.percent, soilC_10_30cm$sample.ID, mean)), totN.percent=as.numeric(tapply(soilC_10_30cm$totN.percent, soilC_10_30cm$sample.ID, mean)), orgC.percent=as.numeric(tapply(soilC_10_30cm$orgC.percent, soilC_10_30cm$sample.ID, mean)), inorgC.percent=as.numeric(tapply(soilC_10_30cm$inorgC.percent, soilC_10_30cm$sample.ID, mean)), CaCO3.percent=as.numeric(tapply(soilC_10_30cm$CaCO3.percent, soilC_10_30cm$sample.ID, mean)), orgC.to.totC=as.numeric(tapply(soilC_10_30cm$orgC.to.totC, soilC_10_30cm$sample.ID, mean)), N.error.percent=as.numeric(tapply(soilC_10_30cm$N.error.percent, soilC_10_30cm$sample.ID, mean)))
sum(soilC_10_30cm$inorgC.percent < 0)
soilC_10_30cm$inorgC.percent[soilC_10_30cm$inorgC.percent < 0] <- 0 #2 samples had undetectable levels of carbonates
soilC_10_30cm$CaCO3.percent[soilC_10_30cm$CaCO3.percent < 0] <- 0
soilC_10_30cm$orgC.to.totC[soilC_10_30cm$orgC.to.totC > 1] <- 1
soilC_10_30cm$orgC.percent[soilC_10_30cm$orgC.percent > soilC_10_30cm$totC.percent] <- soilC_10_30cm$totC.percent[soilC_10_30cm$orgC.percent > soilC_10_30cm$totC.percent]
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
hist(soilC_10_30cm$orgC.percent)
write.csv(soilC_0_10cm, file.path(soilCresults, 'soilC_0_10cm_Camatta.csv'), row.names = FALSE)
write.csv(soilC_10_30cm, file.path(soilCresults, 'soilC_10_30cm_Camatta.csv'), row.names = FALSE)
soilC_0_10cm <- read.csv(file.path(soilCresults, 'soilC_0_10cm_Camatta.csv'), stringsAsFactors = FALSE)
soilC_10_30cm <- read.csv(file.path(soilCresults, 'soilC_10_30cm_Camatta.csv'), stringsAsFactors = FALSE)
dim(soilC_0_10cm)
dim(soilC_10_30cm)
sum(soilC_0_10cm$orgC.percent > soilC_0_10cm$totC.percent) #0
sum(soilC_10_30cm$orgC.percent > soilC_10_30cm$totC.percent) #0

#read in bulk density data
#read-in BD_data
BD_data <- read.csv(file.path(results, 'soil_data', 'soilBD_H2O_frags_2018-09-17.csv'), stringsAsFactors = FALSE)
#split dataset by depth
BD_data_0_10cm <- BD_data[BD_data$depth_code==1 | BD_data$depth_code==3, ]
BD_data_10_30cm <- BD_data[BD_data$depth_code==2,]

#read-in terrain characteristics from S Hogan's drone data
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
names(Mar2017_terrain_3m)
names(Mar2017_terrain_3m) <- c('aspect', 'curvature_mean', 'curvature_plan', 'curvature_profile', 'elevation', 'slope', 'TCI')
#read in 3 m solrad results produced in ArcGIS
list.files(file.path(solradDir)) #all files were sky size 500 x 500 and 64 calc directions
 #all v2 files were sky size 500 x 500 and 64 calc directions
solrad_raster <- raster(file.path(solradDir, 'solrad_3m_filtered.tif'))
solrad_raster <- solrad_raster / 1000
solrad <- shapefile(file.path(solradDir, 'solrad_105pts.shp'))
plot(solrad)
names(solrad)
solrad$point_no <- sampling_pts$point_no #this is correct labelling
text(solrad, labels=solrad$point_no, offset=0.1, pos=1)
#solrad$annual_kwh.m2 <- apply(as.data.frame(solrad[ ,1:52]), 1, sum) / 1000
solrad_df <- as.data.frame(solrad)
colnames(solrad_df)
head(solrad_df)
solrad_df$annual_kwh.m2 <- apply(solrad_df[ ,1:52], 1, sum) / 1000
summary(solrad_df$annual_kwh.m2)
#check values against raster
solrad_df$check_sums <- extract(solrad_raster, coordinates(solrad)[,1:2])
plot(solrad_df$annual_kwh.m2, solrad_df$check_sums) #computational perfection!
Mar2017_terrain_3m$solrad <- solrad_raster
#length(seq.Date(as.Date('Oct_01_2017', '%b_%d_%Y'), as.Date('Apr_15_2018', '%b_%d_%Y'), by='day')) #197 days
#format(as.Date('Oct_01_2017', '%b_%d_%Y'), '%V') #so, T38-T51
#format(as.Date('Apr_15_2018', '%b_%d_%Y'), '%V') #and T0-T15
#solrad_df$seasonal_kwh.m2 <- apply(solrad_df[ ,c(1:16, 39:52)], 1, sum) / 1000
sampling_pts <- merge(sampling_pts, solrad_df[ ,c('annual_kwh.m2', 'point_no')], by='point_no')
#head(sampling_pts)
sampling_pts$slope <- extract(Mar2017_terrain_3m$slope, coordinates(sampling_pts)[,1:2])
sampling_pts$curvature_mean <- extract(Mar2017_terrain_3m$curvature_mean, coordinates(sampling_pts)[,1:2])
sampling_pts$elevation <- extract(Mar2017_terrain_3m$elevation, coordinates(sampling_pts)[,1:2])
sampling_pts$TCI <- extract(Mar2017_terrain_3m$TCI, coordinates(sampling_pts)[,1:2])
sampling_pts$curvature_plan <- extract(Mar2017_terrain_3m$curvature_plan, coordinates(sampling_pts)[,1:2])
sampling_pts$curvature_profile <- extract(Mar2017_terrain_3m$curvature_profile, coordinates(sampling_pts)[,1:2])
sampling_pts$aspect <- extract(Mar2017_terrain_3m$aspect, coordinates(sampling_pts)[,1:2])
soilC_0_10cm$point_no <- as.integer(gsub('_1', '', soilC_0_10cm$sample.ID))
soilC_10_30cm$point_no <- as.integer(gsub('_2', '', soilC_10_30cm$sample.ID))

#now 10-30 cm soils
soil_10_30cm_shp <- merge(sampling_pts, soilC_10_30cm, by='point_no')
soil_10_30cm_shp <- merge(soil_10_30cm_shp, soil_chars_10_30, by='point_no')
soil_10_30cm_shp <- merge(soil_10_30cm_shp, BD_data_10_30cm, by='point_no')
soil_10_30cm_shp$kgOrgC.m2 <- soil_10_30cm_shp$orgC.percent * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100
soil_10_30cm_shp$kgClay.m2 <- soil_10_30cm_shp$CLAY * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100
soil_10_30cm_shp$kgIC.m2 <- soil_10_30cm_shp$inorgC.percent * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100
soil_10_30cm_shp$gP.m2 <- soil_10_30cm_shp$HCO3_P * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100 * 0.1 #see notes
soil_10_30cm_shp$kgTN.m2 <- soil_10_30cm_shp$totN.percent * soil_10_30cm_shp$bulk_density_g_cm3 * 2 * (100 - soil_10_30cm_shp$frags_vol_perc) / 100
soil_10_30cm_shp$energy_colors <- ifelse(soil_10_30cm_shp$annual_kwh.m2 <= 1200, 'blue', ifelse(soil_10_30cm_shp$annual_kwh.m2 > 1200 & soil_10_30cm_shp$annual_kwh.m2 < 1410, 'orange2', 'red3'))
soil_10_30cm_shp$WMPD_mm <- (soil_10_30cm_shp$CLAY * 0.001 ) / 100 + (soil_10_30cm_shp$SILT * 0.026) / 100 + (soil_10_30cm_shp$SAND * 1.025) / 100
plot(soil_10_30cm_shp$kgOrgC.m2, soil_10_30cm_shp$kgTN.m2)
summary(lm(kgOrgC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.50, p.val=.002, all params sig.
summary(lm(orgC.percent ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.39
summary(lm(kgIC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.48
summary(lm(kgIC.m2 ~ annual_kwh.m2, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.11
summary(lm(kgIC.m2 ~ elevation, data=as.data.frame(soil_10_30cm_shp)))
summary(lm(kgClay.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.31, p.val=0.004; elevation most sig
summary(lm(WMPD_mm ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.28, p.val < 0.001; elevation most sig
summary(lm(gP.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.28, p.val=0.01; aspect most sig
summary(lm(gP.m2 ~ kgIC.m2, data = as.data.frame(soil_10_30cm_shp))) #r^2=0.35
summary(lm(kgOrgC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation + kgClay.m2, data=as.data.frame(soil_10_30cm_shp))) #r2=0.55 p.val < .001
summary(lm(kgOrgC.m2 ~ slope + curvature_mean + kgClay.m2, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.43 p.val < 0.001; all params sig
summary(lm(kgOrgC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation + WMPD_mm, data=as.data.frame(soil_10_30cm_shp))) #r^2=0.53
shapefile(soil_10_30cm_shp, file.path(soilCresults, 'shapefiles', 'soil_10_30cm.shp'), overwrite=TRUE)
write.csv(as.data.frame(soil_10_30cm_shp), file.path(soilCresults, 'shapefiles', 'soil_10_30cm_df.csv'), row.names = FALSE)

#now 0_10 since points 2 and 3 will depend on 10_30 data
soil_0_10cm_shp <- merge(sampling_pts, soilC_0_10cm, by='point_no')
soil_0_10cm_shp <- merge(soil_0_10cm_shp, soil_chars_0_10, by='point_no')
soil_0_10cm_shp <- merge(soil_0_10cm_shp, BD_data_0_10cm, by='point_no')
#correct points 2 and 3 data for these:
names(soil_0_10cm_shp)
vars_to_fix <- c("totC.percent", "totN.percent", "orgC.percent", "inorgC.percent", "CaCO3.percent", "OM", "ENR", "P1", "HCO3_P", "PH", "K", "MG", "CA", "NA.", "CEC", "K_PCT", "MG_PCT", "CA_PCT", "H_PCT", "NA_PCT", "S", "SAND", "SILT", "CLAY", "bulk_density_g_cm3")
soil_0_10cm_shp$orgC.percent[2]
soil_10_30cm_shp$orgC.percent[2]
# coeffs <- lm(c(soil_10_30cm_shp$orgC.percent[2], soil_0_10cm_shp$orgC.percent[2]) ~ c(2.5, 20))$coefficients
# coeffs
# coeffs[1] + coeffs[2] * 5 #to get mid-point value, assuming linear interpolation from 2.5 cm to 20 cm
fix_0_5_samples <- function(varname, point) {
  coeffs <- lm(c(soil_10_30cm_shp[[varname]][point], soil_0_10cm_shp[[varname]][point]) ~ c(2.5, 20))$coefficients #2.5 is midpoint for 0-5 cm sample; 20 is midpoint for 10-30 cm sample
  coeffs[1] + coeffs[2] * 5 #5 is midpoint for 0-10 cm sample; so data is corrected to a 0-10 cm slice based on 0-5 and 10-30 cm data
}
#fix point 2, which was a 0-5 cm sample
for (i in seq_along(vars_to_fix)) {
  soil_0_10cm_shp[2, vars_to_fix[i]] <- fix_0_5_samples(vars_to_fix[i], 2)
}
#fix point 3, which was a 0-5 cm sample
for (i in seq_along(vars_to_fix)) {
  soil_0_10cm_shp[3, vars_to_fix[i]] <- fix_0_5_samples(vars_to_fix[i], 3)
}
as.data.frame(soil_0_10cm_shp)[2,]
soil_0_10cm_shp$orgC.percent[2] + soil_0_10cm_shp$inorgC.percent[2]
soil_0_10cm_shp$totC.percent[2]
as.data.frame(soil_0_10cm_shp)[3,]
soil_0_10cm_shp$kgOrgC.m2 <- soil_0_10cm_shp$orgC.percent * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100
soil_0_10cm_shp$kgIC.m2 <- soil_0_10cm_shp$inorgC.percent * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100
soil_0_10cm_shp$kgClay.m2 <- soil_0_10cm_shp$CLAY * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100
soil_0_10cm_shp$gP.m2 <- soil_0_10cm_shp$HCO3_P * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100 * 0.1 #see notes; P1 is Bray and HCO3-P is Olsen, better for high pH soils
soil_0_10cm_shp$energy_colors <- ifelse(soil_0_10cm_shp$annual_kwh.m2 <= 1200, 'blue', ifelse(soil_0_10cm_shp$annual_kwh.m2 > 1200 & soil_0_10cm_shp$annual_kwh.m2 < 1410, 'orange2', 'red3'))
soil_0_10cm_shp$WMPD_mm <- (soil_0_10cm_shp$CLAY * 0.001 ) / 100 + (soil_0_10cm_shp$SILT * 0.026) / 100 + (soil_0_10cm_shp$SAND * 1.025) / 100
soil_0_10cm_shp$kgTN.m2 <- soil_0_10cm_shp$totN.percent * soil_0_10cm_shp$bulk_density_g_cm3 * (100 - soil_0_10cm_shp$frags_vol_perc) / 100
plot(soil_0_10cm_shp$kgOrgC.m2, soil_0_10cm_shp$kgTN.m2)
text(soil_0_10cm_shp$kgOrgC.m2, soil_0_10cm_shp$kgTN.m2, labels=soil_0_10cm_shp$point_no, offset=0.1, pos=1)
summary(lm(kgOrgC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_0_10cm_shp))) #r^2=0.21, p.val=.002, slope, solrad, and curvature all sig
summary(lm(orgC.percent ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_0_10cm_shp))) #r2=0.16
summary(lm(kgIC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_0_10cm_shp))) #NS p.val=0.28
summary(lm(kgIC.m2 ~ annual_kwh.m2, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(kgIC.m2 ~ elevation, data=as.data.frame(soil_0_10cm_shp)))
summary(lm(kgClay.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_0_10cm_shp))) #r^2=0.14, p.val=0.004; elevation most sig
summary(lm(WMPD_mm ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_0_10cm_shp))) #r^2=0.23, p.val < 0.001; elevation most sig
summary(lm(gP.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation, data=as.data.frame(soil_0_10cm_shp))) #r^2=0.12, p.val=0.01; aspect most sig
summary(lm(kgOrgC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation + kgClay.m2, data=as.data.frame(soil_0_10cm_shp))) #r2=0.26 p.val < .001
summary(lm(kgOrgC.m2 ~ slope + curvature_mean + kgClay.m2, data=as.data.frame(soil_0_10cm_shp))) #r^2=0.21 p.val < 0.001; all params sig
summary(lm(kgOrgC.m2 ~ slope + annual_kwh.m2 + curvature_mean + elevation + WMPD_mm, data=as.data.frame(soil_0_10cm_shp))) #r^2=0.24
shapefile(soil_0_10cm_shp, file.path(soilCresults, 'shapefiles', 'soil_0_10cm.shp'), overwrite=TRUE)
write.csv(as.data.frame(soil_0_10cm_shp), file.path(soilCresults, 'shapefiles', 'soil_0_10cm_df.csv'), row.names = FALSE)
hist(soil_0_10cm_shp$kgOrgC.m2)
plot(soil_0_10cm_shp, cex=soil_0_10cm_shp$kgOrgC.m2/2, pch=2, col=soil_0_10cm_shp$energy_colors)
text(soil_0_10cm_shp, labels=soil_0_10cm_shp$point_no, offset=0.2, pos=1)
plot(soil_0_10cm_shp[soil_0_10cm_shp$orgC.to.totC>0.9,], add=TRUE, pch=20, cex=0.8)


#read-in csv files there to combine into one master file
names(soil_0_10cm_shp)
soil_0_10cm_shp$point_no - soil_10_30cm_shp$point_no
soil_0_30cm_shp <- soil_0_10cm_shp[,1:15]
names(soil_0_30cm_shp)
soil_0_30cm_shp$kgOrgC.m2 <- soil_0_10cm_shp$kgOrgC.m2 + soil_10_30cm_shp$kgOrgC.m2
soil_0_30cm_shp$kgTN.m2 <- soil_0_10cm_shp$kgTN.m2 + soil_10_30cm_shp$kgTN.m2
soil_0_30cm_shp$kgClay.m2 <- soil_0_10cm_shp$kgClay.m2 + soil_10_30cm_shp$kgClay.m2
soil_0_30cm_shp$kgIC.m2 <- soil_0_10cm_shp$kgIC.m2 + soil_10_30cm_shp$kgIC.m2
soil_0_30cm_shp$gP.m2 <- soil_0_10cm_shp$gP.m2 + soil_10_30cm_shp$gP.m2
soil_0_30cm_shp$WMPD_mm <- (soil_0_10cm_shp$WMPD_mm * 10 + soil_10_30cm_shp$WMPD_mm * 20) / 30 
soil_0_30cm_shp$sand_wtd <- (10*soil_0_10cm_shp$SAND + 20*soil_10_30cm_shp$SAND) / 30
soil_0_30cm_shp$silt_wtd <- (10*soil_0_10cm_shp$SILT + 20*soil_10_30cm_shp$SILT) / 30
soil_0_30cm_shp$clay_wtd <- (10*soil_0_10cm_shp$CLAY + 20*soil_10_30cm_shp$CLAY) / 30
plot(soil_0_30cm_shp, cex=soil_0_30cm_shp$kgOrgC.m2/3, pch=20)
plot(soil_0_30cm_shp, cex=soil_0_30cm_shp$kgIC.m2/2, pch=20)
hist(soil_0_30cm_shp$kgOrgC.m2)
summary(soil_0_30cm_shp$kgOrgC.m2)
hist(soil_0_30cm_shp$kgIC.m2)
hist(soil_0_30cm_shp$kgTN.m2)
hist(soil_0_30cm_shp$gP.m2)
hist(soil_0_30cm_shp$WMPD_mm)
hist(soil_0_30cm_shp$sand_wtd)
hist(soil_0_30cm_shp$silt_wtd)
hist(soil_0_30cm_shp$clay_wtd)
plot(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$kgTN.m2)
plot(soil_0_30cm_shp$annual_kwh.m2, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$elevation, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$curvature_mean, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$slope, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$TCI, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$curvature_plan, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$curvature_profile, soil_0_30cm_shp$kgOrgC.m2)
plot(soil_0_30cm_shp$kgOrgC.m2, soil_0_30cm_shp$kgIC.m2)
plot(soil_0_30cm_shp$TCI, soil_0_30cm_shp$curvature_mean)
plot(soil_0_30cm_shp$slope, soil_0_30cm_shp$annual_kwh.m2)
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + annual_kwh.m2 + slope, data=soil_0_30cm_shp)) #r^2 =0.42, all sig
summary(lm(kgOrgC.m2 ~ kgIC.m2, data=soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope+ annual_kwh.m2 + WMPD_mm, data=soil_0_30cm_shp)) #r^2=0.48
plot(lm(kgOrgC.m2 ~ curvature_mean + elevation + annual_kwh.m2 + slope, data=soil_0_30cm_shp))
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope+ annual_kwh.m2 + WMPD_mm, data=soil_0_30cm_shp[-2,])) #point 2 is an outlier
plot(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope+ annual_kwh.m2 + WMPD_mm, data=soil_0_30cm_shp[-2,])) #point 2 is an outlier
orgC_lm <- lm(kgOrgC.m2 ~ curvature_mean + elevation + slope+ annual_kwh.m2 + WMPD_mm, data=soil_0_30cm_shp)
summary(lm(kgOrgC.m2 ~ curvature_mean + elevation + slope+ annual_kwh.m2 + kgClay.m2, data=soil_0_30cm_shp))
summary(lm(kgIC.m2 ~ curvature_mean + elevation + slope+ annual_kwh.m2 + kgClay.m2, data=soil_0_30cm_shp))
write.csv(as.data.frame(soil_0_30cm_shp), file.path(soilCresults, 'shapefiles', 'soil_0_30cm_df.csv'), row.names = FALSE)
shapefile(soil_0_30cm_shp, file.path(soilCresults, 'shapefiles', 'soil_0_30cm.shp'), overwrite=TRUE)
