##### iowa urban county level estimates #####
library(reshape2)
library(plyr)
setwd("/Users/Xin/Research/survey/")
dat_nriia <- read.csv("/Volumes/cssm_groups$/Grads/XinWang/NRI/data/Pgen/PgenAllWithLonLatIDIA.csv")

#### calculate urban area in each year and each county 
names <- paste("BU",c(1997,2000:2012),sep="")
urbania <- dat_nriia[,c("COUNTY","WEIGHT",names)]
urbania[,names] <- (urbania[,names]==7)*urbania$WEIGHT

urbania_nri <- ddply(urbania, .(COUNTY), numcolwise(sum)) ### urban total in 100 acreage

#### calculate jk variance 
years <- c(1997, 2000:2012)
for(j in 1:length(years))
{
  yearj <- years[j]
  if(yearj <= 97)
  {
    weight_names <- paste("WEIGHT",1:29,"_","82","_","97",sep="")
  }
  if(yearj >=2000 & yearj <= 2003)
  {
    weight_names <- paste("WEIGHT",1:29,"_","00","_","03",sep="")
  }
  if(yearj >= 2004 & yearj <= 2012)
  {
    weight_names <- paste("WEIGHT",1:29,"_","04","_","12",sep="")
  }
  buj <- paste("BU",years[j],sep="")
  datj <- dat_nriia[,c("COUNTY",buj,weight_names)]
  datj[,weight_names] <- datj[,weight_names]*(datj[,buj]==7)/100
  datcountyj <- ddply(datj[,-2],.(COUNTY), numcolwise(sum))
  
  meanj <- rowMeans(datcountyj[,-1])
  varj <- 28*rowSums((datcountyj[,-1] - meanj)^2)/29
  urbania_nri[,paste("BUvar",yearj,sep="")] <- varj
}

ia_urban_area <- read.csv("data/ia_urban_area_800m_ts.csv",stringsAsFactors = FALSE)



urbania_ts <- ddply(ia_urban_area[,c(paste("A",1985:2015,sep=""),"COUNTYFP","NAME")],.(COUNTYFP,NAME),numcolwise(sum))

urbania_ts <- arrange(urbania_ts, COUNTYFP)
urbania_nri <- arrange(urbania_nri, COUNTY)

constkm2ac <- 247.105/100

plot(urbania_ts$A1997*constkm2ac, urbania_nri$BU1997)
abline(0,1)

ia_district <- read.csv("data/ia_district.csv",stringsAsFactors = FALSE)
sum(urbania_nri$COUNTY %in% ia_district$COUNTY)
sum(ia_district$COUNTY %in% urbania_nri$COUNTY)


colnames(urbania_ts) <- c("COUNTYFP","NAME",1985:2015)


urbania_nri_est <- urbania_nri[,c("COUNTY",paste("BU",c(1997,2000:2012),sep=""))]
colnames(urbania_nri_est) <- c("COUNTY",c(1997,2000:2012))
urbania_nri_est_mt <- melt(urbania_nri_est,id.vars = "COUNTY",variable.name = "year",value.name = "area_nri")

urbania_nri_var <- urbania_nri[,c("COUNTY",paste("BUvar",c(1997,2000:2012),sep=""))]
colnames(urbania_nri_var) <- c("COUNTY",c(1997,2000:2012))
urbania_nri_var_mt <- melt(urbania_nri_var,id.vars = "COUNTY",variable.name = "year",value.name = "var_nri")


urbania_ts_mt <- melt(urbania_ts,id.vars = c("COUNTYFP","NAME"),variable.name = "year",value.name = "area_ts")
urbania_ts_mt$area_ts <- urbania_ts_mt$area_ts*constkm2ac


urbania <- merge(urbania_ts_mt, urbania_nri_est_mt, by.y =c("COUNTY","year"),by.x = c("COUNTYFP","year"),all.y =TRUE)
urbania <- merge(urbania, urbania_nri_var_mt, by.x = c("COUNTYFP","year"),by.y = c("COUNTY","year"))

#urbania <- merge(urbania, ia_district, by = "COUNTY")
urbania <- arrange(urbania, year, COUNTYFP)

write.csv(urbania, "data/urbania.csv",row.names = FALSE)

### regression in each regression 

urbania <- read.csv("data/urbania.csv",stringsAsFactors = FALSE)

ggplot(urbania[urbania$year==2012,], aes(x= area_ts, y = area_nri)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x)
  geom_abline(slope = 1,intercept = 0) + 
  theme_bw()

pdf("doc/figures/NRI2012.pdf",width = 6,height = 5)
ggplot(urbania[urbania$year==2012,], aes(x= area_ts, y = area_nri)) + 
  geom_point()+  ggtitle("NRI 2012")+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept = 0, slope = 1)
dev.off()

pdf("doc/figures/NRI2011.pdf",width = 6,height = 5)
ggplot(urbania[urbania$year==2011,], aes(x= area_ts, y = area_nri)) + 
  geom_point()+  ggtitle("NRI 2011")+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept = 0, slope = 1)
dev.off()




numdistrict <- length(unique(urbania$DISTRICT))
district <- sort(unique(urbania$DISTRICT))

coefm <- matrix(0, numdistrict, 2)
for(j in 1:numdistrict)
{
  coefm[j,] <- coef(lm(area_nri ~ area_ts, urbania[urbania$year==1997 & urbania$DISTRICT==district[j],]))
}


library(Spgr)

cvec0 <- rep(1,99*98/2)

z0 <- matrix(1, ncol = 1, nrow = 99)

dat2012 <- dplyr::filter(urbania, year == 2012)
betam01 <- cal_initialx(y = scale(dat2012$area_nri), x = cbind(1,scale(dat2012$area_ts)),lam0 = 0.0001)

betam02 <- cal_initial(y = scale(dat2012$area_nri), z = z0,x = as.matrix(scale(dat2012$area_ts)),lam0 = 0.0001)
res2 <- Spgr(y = scale(dat2012$area_nri),z = z0,x = as.matrix(scale(dat2012$area_ts)),weights = cvec0,betam0 = betam02,maxiter = 2000,lam = 0.1)


lamv <- seq(0.05,0.2,by = 0.01)
bicv1 <- rep(0, length(lamv))  
eta.mat1 <- rep(0, length(lamv))
group1 <- matrix(0, 99, length(lamv))
beta.array1 <- array(0, dim = c(99,1,length(lamv)))
flagv1 <- rep(0,length(lamv))
niterv1 <- rep(0,length(lamv))

betam01 <- cal_initialx(y = scale(dat2012$area_nri), x = cbind(1,scale(dat2012$area_ts)),lam0 = 0.0001)

betam02 <- cal_initial(y = scale(dat2012$area_nri), z = z0,x = as.matrix(scale(dat2012$area_ts)),lam0 = 0.0001)

for(l in 1:length(lamv))
{
  # resl <- Spgrx(y = scale(dat2012$area_nri),x = cbind(1, scale(dat2012$area_ts)),weights = cvec0,betam0 = betam01,maxiter = 2000,lam = lamv[l])
  
  resl <- Spgr(y = scale(dat2012$area_nri),z = z0,x = as.matrix(scale(dat2012$area_ts)),weights = cvec0,betam0 = betam02,maxiter = 2000,lam = lamv[l])
  
  #bicv1[l] <- BICcx(obj = resl,y = scale(dat2012$area_nri),x = cbind(1, scale(dat2012$area_ts)),c0 = 1)
  
  #bicv1[l] <- BIClogx(obj = resl,y = scale(dat2012$area_nri),x = cbind(1, scale(dat2012$area_ts)))
  
  bicv1[l] <- BIClog(obj = resl,y = scale(dat2012$area_nri),z = z0,x =  as.matrix(scale(dat2012$area_ts)))
  
  beta.array1[,,l] <- resl$beta
  betam02 <- resl$beta
  eta.mat1[l] <- resl$eta
  group1[,l] <- resl$group
  flagv1[l] <- resl$flag
  niterv1[l] <- resl$niteration
  
  print(l)
}

plot(bicv1)

plot(beta.array1[,,which.min(bicv1)])


### without intercept

lamv <- seq(0.05,0.3,by = 0.01)
bicv3 <- rep(0, length(lamv))  
group3 <- matrix(0, 99, length(lamv))
beta.array3 <- array(0, dim = c(99,1,length(lamv)))
flagv3 <- rep(0,length(lamv))
niterv3 <- rep(0,length(lamv))

betam03 <- cal_initialx(y = scale(dat2012$area_nri), x = as.matrix(scale(dat2012$area_ts)),lam0 = 0.0001)


for(l in 1:length(lamv))
{
  resl <- Spgrx(y = scale(dat2012$area_nri),x = as.matrix(scale(dat2012$area_ts)),weights = cvec0,betam0 = betam03,maxiter = 2000,lam = lamv[l])
  
  bicv3[l] <- BIClogx(obj = resl,y = scale(dat2012$area_nri),x =  as.matrix(scale(dat2012$area_ts)))
  
  beta.array3[,,l] <- resl$beta
  betam03 <- resl$beta

  group3[,l] <- resl$group
  flagv3[l] <- resl$flag
  niterv3[l] <- resl$niteration
  
  print(l)
}

plot(bicv3)

plot(beta.array3[,,which.min(bicv3)])


groupest <- group3[,which.min(bicv3)]

plot(dat2012$area_ts, dat2012$area_nri, col=groupest)


