############## segment information from NRI
library(plyr)
library(dplyr)
datia <- read.csv("/Volumes/cssm_groups$/Grads/XinWang/NRI/data/2012/PgenAllIA.csv")
nseg <- ddply(datia,.(COUNTY),summarize, nseg = length(unique(PSU)))

year = 2000:2012

nsegmat <- matrix(0, nrow(nseg), length(year))
for(j in 1:length(year))
{
  namej <- paste("SAMPLE",year[j],sep="")
  nsegj <- ddply(datia[datia[,namej]==1,],.(COUNTY),summarize, nseg = length(unique(PSU)))
  nsegmat[,j] <- nsegj[,2]
}

total1 <- ddply(datia, .(COUNTY), summarize, area = sum(WEIGHT))

IAset <- read.csv("/Volumes/cssm_groups$/Grads/XinWang/NRI/data/IASet5.csv",header = TRUE)
IAset$HUCCO <- substr(IAset$HUCCO,3,5)
IAset$HUCCO <- as.numeric(IAset$HUCCO)
IAarea <- ddply(IAset,.(HUCCO),summarize,area = sum(HUCCOACRES))

IAarea$area = IAarea$area * sum(total1[,2])/sum(IAarea$area)

write.csv(IAarea, file = "Research/survey/IAarea.csv",row.names = FALSE)


IAarea <- read.csv("data/IAarea.csv",stringsAsFactors = FALSE)

