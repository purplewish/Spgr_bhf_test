########  

library(plyr)
dat <- read.csv("data/datia_800m.csv",stringsAsFactors = FALSE)
dat1 <- read.csv("data/datia1_800m.csv",stringsAsFactors = FALSE)

datcounty <- ddply(dat1, .(county), summarize, numseg = length(county),
                   area2001 = sum(A2001), area2011 = sum(A2011), 
                   mean2001 = mean(A2001), mean2011 = mean(A2011),
                   mean1985 = mean(A1985), mean2015 = mean(A2015)) ## population 

ia_urban_area <- read.csv("Research/survey/data/ia_urban_area_800m.csv",stringsAsFactors = FALSE)

### srs 

source("Research/survey/code/simsrs.R")

Nv <- datcounty$numseg
sr <- 0.03
popn <- datcounty

M <- 100


res1 <- simsrs(100, dat,popn, 0.03)
res2 <- simsrs(100, dat, popn, 0.02)
res3 <- simsrs(100, dat, popn, 0.01)

save(res1,res2, res3, res21, file = "Research/survey/result/urbansrs.RData")

boxplot(res1$rmse)
boxplot(res2$rmse)
boxplot(res3$rmse)

boxplot(res1$rmse[,c(2,3,5,6)])
boxplot(res2$rmse[,c(2,3,5,6)])
boxplot(res3$rmse[,c(2,3,5,6)])



##### 1985 vs 2015 
ia_urban_area <- read.csv("Research/survey/data/ia_urban_area_800m.csv",stringsAsFactors = FALSE)
datcounty1 <- ddply(ia_urban_area, .(county), summarize, numseg = length(county),
                   mean1985 = mean(A1985), mean2015 = mean(A2015), 
                   mean2001 = mean(A2001), mean2011 = mean(A2011))

popn <- datcounty1[,c("county","numseg")]
meanx <- datcounty1[,c("county","mean1985")]
meany <- datcounty1[,c("county","mean2015")]
group <- c(rep(1,50),rep(2,49))
srh <- 0.03
srl <- 0.03/5

meanx01 <- datcounty1[,c("county","mean2001")]
meany11 <- datcounty1[,c("county","mean2011")]



res21 <- simsrs2(M = 100,dat = ia_urban_area,popn = popn,xname = "A2001",yname = "A2011",
                 meanx = meanx01,meany = meany11,srh = 0.03,srl = 0.03/5)
boxplot(res21$rmse)
boxplot(res21$rmse[,c(2,3,5,6)])


res22 <- simsrs2(M = 100,dat = ia_urban_area,popn = popn,xname = "A1985",yname = "A2015",
                 meanx = meanx,meany = meany,srh = 0.02,srl = 0.02/5)
boxplot(res22$rmse)

res23 <- simsrs2(M = 100,dat = ia_urban_area,popn = popn,xname = "A1985",yname = "A2015",
                 meanx = meanx,meany = meany,srh = 0.01,srl = 0.01/5)
boxplot(res23$rmse)

save(res21, res22, res23, file = "Research/survey/result/urbansrs2.RData")
