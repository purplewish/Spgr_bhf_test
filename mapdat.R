####
library(raster)
library(dplyr)
library(plyr)
library(ggplot2)
library(rgdal)
library(maptools)

setwd("/Users/Xin/Research/NRI_urban/")
### county boundry ###
county_boundry <- readOGR("data/cb_2016_us_county_5m/",layer = "cb_2016_us_county_5m")
ia_boundry <- county_boundry[county_boundry$STATEFP==19,]
ia_dat <- ia_boundry@data
ia_dat <- cbind(ID = 1:nrow(ia_dat), ia_dat)
ia_fort <- fortify(ia_boundry, region="NAME") 
ia_boundry <- as(ia_boundry,"SpatialPolygons")



ggplot(data=ia_fort, aes(long, lat, group = group)) + 
  geom_polygon(color = "black", fill = "white") + theme_bw()


### map all data to county level for all years
#datall <- read.csv("Research/survey/data/IA_TS_800m/IA_UrbanArea_TS.csv")
year <- seq(1985,2015,by = 1)
nyear <- length(year)

crsboundry <- crs(county_boundry)

### get index based on 1985
dat85r <- raster("data/IA_TS_800m/IA_UrbanArea_1985.tif")
dat85 <- rasterToPoints(dat85r,spatial = TRUE)
dat85 <- spTransform(dat85,CRSobj = crs(county_boundry))

crsts <- crs(dat85r)


index85 <- over(dat85, ia_boundry) ### find county 
## the numbers are the orders in the data not the fips
indexnn <- !is.na(index85) ## not NA index
index85n <- index85[indexnn] # county index 
df85 <- as.data.frame(dat85)[indexnn,] ## data without NA

ia_urban_area <- as.data.frame(matrix(0, length(index85n), 3 + nyear))
colnames(ia_urban_area) <- c("index","x","y",paste("A",year,sep=""))
ia_urban_area$index <- index85n
ia_urban_area[,c("x","y","A1985")] <- df85[,c(2,3,1)]


for(j in 2:nyear)
{
  namej <- paste("data/IA_TS_800m/IA_UrbanArea_",year[j],".tif",sep="")
  datjr <- raster(namej)
  datj <- rasterToPoints(datjr,spatial = TRUE)
  datj <- spTransform(datj,CRSobj = crsboundry)
  dfj <- as.data.frame(datj)
  
  ia_urban_area[,j+3] <- dfj[indexnn,1]
}

### replace county numbers with fips 

ia_urban_area <- merge(ia_urban_area, ia_dat[,c("ID","COUNTYFP","GEOID","NAME")],by.x = "index",by.y ="ID")


ia_urban_area$COUNTYFP <- as.numeric(as.character(ia_urban_area$COUNTYFP))
ia_urban_area$GEOID <- as.numeric(as.character(ia_urban_area$GEOID))

ia_urban_area <- arrange(ia_urban_area, COUNTYFP)
write.csv(ia_urban_area, "data/ia_urban_area_800m.csv",row.names = FALSE)


#write.csv(dat, "Research/survey/data/datia_800m.csv",row.names = FALSE)
#write.csv(dat1, "Research/survey/data/datia1_800m.csv",row.names = FALSE)

#datcounty <- ddply(dat, .(county), summarize, area01 = sum(A2001), 
#                 area11 = sum(A2011), numseg = length(county)) ## population 



temp <- raster("Research/survey/IA_TS_800m/IA_UrbanArea_2001.tif")
temp <- rasterToPoints(temp)
temp1 <- as.data.frame(temp, xy = TRUE)



####  map crs based on ts #####
dat85r <- raster("data/IA_TS_800m/IA_UrbanArea_1985.tif")
dat85 <- rasterToPoints(dat85r,spatial = TRUE)

crsts <- crs(dat85r)

ia_boundry1 <- spTransform(ia_boundry, crsts)

index85 <- over(dat85, ia_boundry1) ### find county 
## the numbers are the orders in the data not the fips
indexnn <- !is.na(index85) ## not NA index
index85n <- index85[indexnn] # county index 
df85 <- as.data.frame(dat85)[indexnn,] ## data without NA

ia_urban_area <- as.data.frame(matrix(0, length(index85n), 3 + nyear))
colnames(ia_urban_area) <- c("index","x","y",paste("A",year,sep=""))
ia_urban_area$index <- index85n
ia_urban_area[,c("x","y","A1985")] <- df85[,c(2,3,1)]


for(j in 2:nyear)
{
  namej <- paste("data/IA_TS_800m/IA_UrbanArea_",year[j],".tif",sep="")
  datjr <- raster(namej)
  datj <- rasterToPoints(datjr,spatial = TRUE)
  dfj <- as.data.frame(datj)
  
  ia_urban_area[,j+3] <- dfj[indexnn,1]
}

### replace county numbers with fips 

ia_urban_area <- merge(ia_urban_area, ia_dat[,c("ID","COUNTYFP","GEOID","NAME")],by.x = "index",by.y ="ID")


ia_urban_area$COUNTYFP <- as.numeric(as.character(ia_urban_area$COUNTYFP))
ia_urban_area$GEOID <- as.numeric(as.character(ia_urban_area$GEOID))

ia_urban_area <- arrange(ia_urban_area, COUNTYFP)
write.csv(ia_urban_area, "data/ia_urban_area_800m_ts.csv",row.names = FALSE)

temp2 <- read.csv("data/ia_urban_area_800m_ts.csv",stringsAsFactors = FALSE)
temp3 <- read.csv("data/ia_urban_area_800m.csv",stringsAsFactors = FALSE)

