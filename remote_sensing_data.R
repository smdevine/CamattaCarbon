library(raster)
mainDir <- 'C:/Users/smdevine/Desktop/rangeland project/Grace forage estimates'
forageDir <- 'C:/Users/smdevine/Desktop/rangeland project/clip_plots/summaries'
sensorDir <- 'C:/Users/smdevine/Desktop/rangeland project/soilmoisture/sensor_coordinates'
list.files(forageDir)
list.files(file.path(mainDir, 'APAR only', '2017'))
biomass_2017 <- stack(lapply(list.files(file.path(mainDir, 'APAR only', '2017'), full.names = TRUE), raster))
biomass_2018 <- stack(lapply(list.files(file.path(mainDir, 'APAR only', '2018'), full.names = TRUE), raster))

biomass_2017_full_model <- stack(lapply(list.files(file.path(mainDir, 'Model II', '2017'), full.names = TRUE), raster))
biomass_2018_full_model <- stack(lapply(list.files(file.path(mainDir, 'Model II', '2018'), full.names = TRUE), raster))
cellStats(biomass_2017, 'mean')
#Biomass_2017.02.15_APAR Biomass_2017.03.14_APAR Biomass_2017.04.10_APAR 
#1039.508                1816.917                2708.223 
#Biomass_2017.04.30_APAR 
#3083.160 
cellStats(biomass_2017_full_model, 'mean')
#Biomass_2017.02.15_APAR Biomass_2017.03.14_APAR Biomass_2017.04.10_APAR 
#1168.512                1824.496                3272.777 
#Biomass_2017.04.29_APAR 
#3586.956
cellStats(biomass_2018_full_model, 'mean')
#Biomass_2018.01.15_APAR Biomass_2018.02.15_APAR Biomass_2018.03.22_APAR 
#280.3398                482.4730                746.0248 
#Biomass_2018.04.14_APAR 
#1054.3800
forage_data <- read.csv(file.path(forageDir, 'forage2017_2018_summary.csv'), stringsAsFactors = FALSE)
sensor_locations <- shapefile(file.path(sensorDir, '5TM_sensor_locations_Camatta.shp')) 
biomass_2017_sensors <- extract(biomass_2017, coordinates(sensor_locations)[,1:2], fun=mean, buffer=3, df=TRUE)
biomass_2018_sensors <- extract(biomass_2018, coordinates(sensor_locations)[,1:2], fun=mean, buffer=2, df=TRUE)
biomass_2017_full_model_sensors <- extract(biomass_2017_full_model, coordinates(sensor_locations)[,1:2], fun=mean, buffer=3, df=TRUE)
biomass_2018_full_model_sensors <- extract(biomass_2018_full_model, coordinates(sensor_locations)[,1:2], fun=mean, buffer=2, df=TRUE)
plot(forage_data$clp031417, biomass_2017_sensors$Biomass_2017.03.14_APAR)
abline(0,1,lty=2)
abline(lm(biomass_2017_sensors$Biomass_2017.03.14_APAR[-13] ~ forage_data$clp031417[-13]))
summary(lm(biomass_2017_sensors$Biomass_2017.03.14_APAR[-13] ~ forage_data$clp031417[-13])) #r2=0.52 with location 13
plot(forage_data$clp041017, biomass_2017_sensors$Biomass_2017.04.10_APAR)
plot(forage_data$clp050117, biomass_2017_sensors$Biomass_2017.04.30_APAR)
AprGrowth <- biomass_2017$Biomass_2017.04.10_APAR - biomass_2017$Biomass_2017.03.14_APAR

#full model examination
plot(forage_data$clp031417, biomass_2017_full_model_sensors$Biomass_2017.03.14_APAR)
plot(forage_data$clp041017, biomass_2017_full_model_sensors$Biomass_2017.04.10_APAR)
plot(forage_data$clp050117, biomass_2017_full_model_sensors$Biomass_2017.04.29_APAR) #something wrong with 4/30 raster; no data

#2018 plotting, APAR only model
plot(forage_data$clp021518, biomass_2018_sensors$Biomass_2018.02.15_APAR)
plot(forage_data$clp032218, biomass_2018_sensors$Biomass_2018.03.22_APAR)
plot(forage_data$clp041518, biomass_2018_sensors$Biomass_2018.04.14_APAR)

#2018 plotting, full model
plot(forage_data$clp021518, biomass_2018_full_model_sensors$Biomass_2018.02.15_APAR)
plot(forage_data$clp032218, biomass_2018_full_model_sensors$Biomass_2018.03.22_APAR)
plot(forage_data$clp041518, biomass_2018_full_model_sensors$Biomass_2018.04.14_APAR)

#stats for 2017
summary(lm(biomass_2017_sensors$Biomass_2017.02.15_APAR ~ forage_data$clp021517))
summary(lm(biomass_2017_sensors$Biomass_2017.03.14_APAR ~ forage_data$clp031417))
summary(lm(biomass_2017_sensors$Biomass_2017.03.14_APAR[-13] ~ forage_data$clp031417[-13]))
summary(lm(biomass_2017_sensors$Biomass_2017.04.06_APAR ~ forage_data$clp041017))
summary(lm(biomass_2017_sensors$Biomass_2017.04.10_APAR ~ forage_data$clp041017))
summary(lm(biomass_2017_sensors$Biomass_2017.04.30_APAR ~ forage_data$clp050117))
summary(lm(biomass_2017_full_model_sensors$Biomass_2017.02.15_APAR ~ forage_data$clp021517))
summary(lm(biomass_2017_full_model_sensors$Biomass_2017.03.14_APAR ~ forage_data$clp031417))
summary(lm(biomass_2017_full_model_sensors$Biomass_2017.03.14_APAR[-13] ~ forage_data$clp031417[-13]))
summary(lm(biomass_2017_full_model_sensors$Biomass_2017.04.06_APAR ~ forage_data$clp041017))
summary(lm(biomass_2017_full_model_sensors$Biomass_2017.04.10_APAR ~ forage_data$clp041017))
summary(lm(biomass_2017_full_model_sensors$Biomass_2017.04.29_APAR ~ forage_data$clp050117))

#stats for 2018
summary(lm(biomass_2018_sensors$Biomass_2018.02.15_APAR ~ forage_data$clp021518))
summary(lm(biomass_2018_sensors$Biomass_2018.03.19_APAR ~ forage_data$clp032218))
summary(lm(biomass_2018_sensors$Biomass_2018.03.22_APAR ~ forage_data$clp032218))
summary(lm(biomass_2018_sensors$Biomass_2018.04.14_APAR ~ forage_data$clp041518))

summary(lm(biomass_2018_full_model_sensors$Biomass_2018.02.15_APAR ~ forage_data$clp021518))
summary(lm(biomass_2018_full_model_sensors$Biomass_2018.03.19_APAR ~ forage_data$clp032218))
summary(lm(biomass_2018_full_model_sensors$Biomass_2018.03.22_APAR ~ forage_data$clp032218))
summary(lm(biomass_2018_full_model_sensors$Biomass_2018.04.14_APAR ~ forage_data$clp041518))

#quick test with sampling points produced in soilCcalcs
sampling_pts$Mar2017_forage <- extract(biomass_2017$Biomass_2017.03.14_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
sampling_pts$Mar2017_forage_II <- extract(biomass_2017_full_model$Biomass_2017.03.14_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
sampling_pts$Apr2017_forage <- extract(biomass_2017$Biomass_2017.04.10_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
sampling_pts$Apr2017_forage_II <- extract(biomass_2017_full_model$Biomass_2017.04.10_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
plot(sampling_pts$orgC_content, sampling_pts$Mar2017_forage, col=ifelse(sampling_pts$annual_kwh.m2 > 1337, 'red', 'blue'))
plot(sampling_pts$orgC_content, sampling_pts$Mar2017_forage_II, col=ifelse(sampling_pts$annual_kwh.m2 > 1337, 'red', 'blue'))
plot(sampling_pts$orgC_content, sampling_pts$Apr2017_forage, col=ifelse(sampling_pts$annual_kwh.m2 > 1337, 'red', 'blue'))
plot(sampling_pts$orgC_content, sampling_pts$Apr2017_forage_II, col=ifelse(sampling_pts$annual_kwh.m2 > 1337, 'red', 'blue'))

sampling_pts$Mar2018_forage <- extract(biomass_2018$Biomass_2018.03.22_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
plot(sampling_pts$orgC_content, sampling_pts$Mar2018_forage, col=ifelse(sampling_pts$annual_kwh.m2 > 1337, 'red', 'blue'))
sampling_pts$Mar2018_forage_II <- extract(biomass_2018_full_model$Biomass_2018.03.22_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
sampling_pts$Apr2018_forage <- extract(biomass_2018$Biomass_2018.04.14_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)
plot(sampling_pts$orgC_content, sampling_pts$Apr2018_forage, col=ifelse(sampling_pts$annual_kwh.m2 > 1337, 'red', 'blue'))
sampling_pts$Apr2018_forage_II <- extract(biomass_2018_full_model$Biomass_2018.04.14_APAR, coordinates(sampling_pts)[,1:2], fun=mean, buffer=2)

summary(lm(Mar2017_forage ~ orgC_content, data = sampling_pts))
summary(lm(Mar2017_forage_II ~ orgC_content, data = sampling_pts)) #r2=0.05, p=0.02
summary(lm(Apr2017_forage ~ orgC_content, data = sampling_pts))
summary(lm(Apr2017_forage_II ~ orgC_content, data = sampling_pts)) #r2=0.105, p<0.001
summary(lm(Mar2018_forage ~ orgC_content, data = sampling_pts))
summary(lm(Mar2018_forage_II ~ orgC_content, data = sampling_pts))
summary(lm(Apr2018_forage ~ orgC_content, data = sampling_pts))
summary(lm(Apr2018_forage_II ~ orgC_content, data = sampling_pts))