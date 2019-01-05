dataDir <- 'C:/Users/smdevine/Desktop/rangeland project/SoilCarbonProject/OtherSites'
list.files(dataDir)
SFREC <- read.csv(file.path(dataDir, 'sfrec_lab_combined.csv'), stringsAsFactors = FALSE)
head(SFREC)
SFREC_0_30 <- SFREC[which(SFREC$top < 30),]
dim(SFREC)
dim(SFREC_0_30)
length(unique(SFREC_0_30$pedon_id)) #119 pedon IDs
SFREC_0_30$weighting_factor <- (ifelse(SFREC_0_30$bottom < 30, SFREC_0_30$bottom, 30) - SFREC_0_30$top) / 30
summary(SFREC_0_30$weighting_factor)
QC_check <- tapply(SFREC_0_30$weighting_factor, SFREC_0_30$pedon_id, sum)
pedons_to_remove <- names(QC_check[QC_check < 1])
which(QC_check < 1)
SFREC_0_30[which(SFREC_0_30$pedon_id=='d10x6'),]
SFREC_0_30[which(SFREC_0_30$pedon_id=='d10x6'),]
SFREC_0_30_clean <- SFREC_0_30[!(SFREC_0_30$pedon_id %in% pedons_to_remove),]
soilorgC_0_30 <- tapply(SFREC_0_30_clean$tc * SFREC_0_30_clean$weighting_factor, SFREC_0_30_clean$pedon_id, sum)
sum(!is.na(soilorgC_0_30)) #71 data points
hist(soilorgC_0_30)
summary(soilorgC_0_30) #mean is 1.6%; Camatta mean is in lower quartile compared to SFREC

SJER <- read.csv(file.path(dataDir, 'sjer-lab-data.csv'), stringsAsFactors = FALSE)
head(SJER)
SJER_0_30 <- SJER[which(SJER$top < 30),]
dim(SJER)
dim(SJER_0_30)
length(unique(SJER_0_30$pedon_id)) #16 pedon IDs as expected
SJER_0_30$weighting_factor <- (ifelse(SJER_0_30$bottom < 30, SJER_0_30$bottom, 30) - SJER_0_30$top) / 30
summary(SJER_0_30$weighting_factor)
QC_check <- tapply(SJER_0_30$weighting_factor, SJER_0_30$pedon_id, sum)
which(QC_check < 1) #all passed
SJER_0_30_clean <- SJER_0_30
#SJER_0_30_clean <- SJER_0_30[!(SJER_0_30$pedon_id %in% pedons_to_remove),]
soilorgC_0_30_SJER <- tapply(SJER_0_30_clean$tc * SJER_0_30_clean$weighting_factor, SJER_0_30_clean$pedon_id, sum)
soilorgC_0_30_SJER
sum(!is.na(soilorgC_0_30_SJER)) #12 data points
hist(soilorgC_0_30_SJER)
summary(soilorgC_0_30_SJER) #mean is 1.07 closer to Camatta
summary(SJER_0_30$tc) #four NAs
SJER_0_30[is.na(SJER_0_30$tc),]
