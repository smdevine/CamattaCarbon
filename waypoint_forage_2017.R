forageDir <- 'C:/Users/smdevine/Desktop/rangeland project/clip_plots'
sensorDir <- 'C:/Users/smdevine/Desktop/rangeland project/soilmoisture/sensor_coordinates'
terrainDir <- 'C:/Users/smdevine/Desktop/rangeland project/terrain_analysis_r_v3'
library(raster)
#read in shapefiles for waypoint and sensor forage data
sensor_forage_sp <- shapefile(file.path(forageDir, 'results', "sensor_forage2017.shp"))
waypoint_forage_sp <- shapefile(file.path(forageDir, "waypoint_forage2017.shp"))
waypoint_forage_data <- read.csv(file.path(forageDir, 'CamattaBiomassWaypointPlotsOnly2017.csv'))
#waypoint forage has erroneous locations for D01 and D02
waypoint_forage_data

all_forage_sp <- rbind(sensor_forage_sp, waypoint_forage_sp)
as.data.frame(all_forage_sp)
shapefile(all_forage_sp, file.path(forageDir, 'summaries', 'forage2017_allpts.shp'))
forage_data <- read.csv(file.path(forageDir, 'summaries', 'forage2017_2018_summary.csv'), stringsAsFactors = FALSE)
forage_data$clp021517 - all_forage_sp$clp021517[1:16]
forage_data$clp031417 - all_forage_sp$clp031417[1:16]
forage_data$clp041017 - all_forage_sp$clp041017[1:16]
rm(forage_data)
Mar2017_terrain_3m <- stack(list.files(file.path(terrainDir, 'filtered_Hogan'), full.names = TRUE))
