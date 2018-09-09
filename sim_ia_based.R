#### iowa regions ###
setwd("Research/NRI_urban/")
ia_district <- read.csv("data/ia_district.csv",stringsAsFactors = FALSE)
# daturban <- read.csv("data/urbania.csv",stringsAsFactors = FALSE)
# 
# ia_district <- merge(ia_district, daturban[daturban$year==2012,c("COUNTYFP","NAME")],by.x = "COUNTY",by.y="COUNTYFP")
# 
# ia_district$NAME <- tolower(ia_district$NAME)
# ia_district[ia_district$NAME == "o'brien","NAME"]<- "obrien"
# 
# write.csv(ia_district, "data/ia_district.csv",row.names = FALSE)

source("code/mapfun.R")
source("code/mcp_scad.R")
source("code/Spgr_sampling_em.R")
library(ggplot2)
library(gridExtra)
library(survey)
library(sae)
library(plyr)
library(Spgr)




### create new groups ####

ia_district$newgroup <- rep(0, nrow(ia_district))
ia_district$newgroup[ia_district$DISTRICT %in% c(10,20,40)] <- 1
ia_district$newgroup[ia_district$DISTRICT %in% c(30,50,60)] <- 2
ia_district$newgroup[ia_district$DISTRICT %in% c(70,80,90)] <- 3

write.csv(ia_district, "data/ia_district.csv",row.names = FALSE)

g1 <- mapfun(group = ia_district$DISTRICT, countyname = ia_district$NAME) +
  ggtitle("USDA district") +
  theme(plot.title = element_text(hjust = 0.5))
g2 <- mapfun(group = ia_district$newgroup, countyname = ia_district$NAME) +
  ggtitle("Three groups") +
  theme(plot.title = element_text(hjust = 0.5))
pdf("doc/figures/iagroup.pdf",width = 10,height = 5)
grid.arrange(g1,g2,ncol=2)
dev.off()

mapfun(group = ia_district$DISTRICT, countyname = ia_district$NAME)+
  geom_point(data = urbants[urbants$COUNTYFP==27,], aes(x=x,y=y))

summary(daturban$area_ts[daturban$year==2011])

urbants <- read.csv("/Users/Xin/Research/NRI_urban/data/ia_urban_area_800m_ts.csv",stringsAsFactors = FALSE)
nsegtotal <- ddply(urbants, .(COUNTYFP,NAME), summarize, nseg = length(COUNTYFP))
nsegtotal <- merge(nsegtotal, ia_district[,c("COUNTY","DISTRICT","newgroup")], by.x = "COUNTYFP",by.y = "COUNTY")
nsegtotal$NAME <- tolower(nsegtotal$NAME)
nsegtotal[nsegtotal$NAME == "o'brien","NAME"]<- "obrien"
write.csv(nsegtotal,"data/nsegtotal.csv",row.names = FALSE)

urban11 <- ddply(urbants, .(COUNTYFP, NAME), summarize, area11 = sum(A2011)) ### 

daturban <- read.csv("data/urbania.csv",stringsAsFactors = FALSE) ### county level, unit 100 acres

urban11 <- merge(nsegtotal, daturban[daturban$year==2011,c("COUNTYFP","area_ts","area_nri")],by = "COUNTYFP")

mapfun(group = urban11$newgroup, countyname = urban11$NAME)

## change units to 160 acres in x 

urban11$area_ts <- urban11$area_ts*100/160

load("data/adjMat.RData")
Cmatia <- adjMat[["19"]]



#### simulation#####

beta <- matrix(c(0.5,0.5,2.5,2.5,4.5,4.5),ncol=2,byrow = TRUE)
#beta <- matrix(c(0.5,0.5,2,2,3.5,3.5),ncol=2,byrow = TRUE)
sdu <- 0.2
rate <- 0.2

sdx <- sd(urban11$area_ts)
meanx <- mean(urban11$area_ts)

nobs <- nrow(urban11)  
nsegcounty <- urban11$nseg
domsize <- urban11[,c("COUNTYFP","nseg")]
nsampled <- round(nsegcounty*rate)

xmat <- cbind(1, (urban11$area_ts - meanx)/sdx) 


# second level, unit is 160 acres

msemat <- matrix(0,50,4)
remsemat <- matrix(0,50, 4)
colnames(msemat) <- colnames(remsemat) <- c("direct","regression","equal","sp")
for(m in 1:50)
{
  set.seed(631109 + m )
  #set.seed(631109)
  repeat{
    ypop <- sdx*(rowSums(xmat * beta[urban11$newgroup,]) + rnorm(nobs)*sdu) + meanx
    if(sum(ypop<=0)==0){break}
  }
  
  nsegurban <- round(ypop)  # number of segments have urban
  
  datsample <- as.data.frame(matrix(0, sum(nsampled), 3))
  colnames(datsample) <- c("y","dom","weights")
  datsample$dom <- rep(urban11$COUNTYFP, nsampled)
  
  for(i in 1:nobs)
  {
    popseg <- rep(0, nsegcounty[i])
    urbanindex <- sample.int( nsegcounty[i],nsegurban[i],replace = FALSE) # segemnts are urban
    popseg[urbanindex] <- 1
    
    sampled <- sample.int(nsegcounty[i], nsampled[i])
    
    datsample$y[datsample$dom == urban11$COUNTYFP[i]] <- popseg[sampled]
    
    datsample$weights[datsample$dom == urban11$COUNTYFP[i]] <- nsegcounty[i]/nsampled[i]
  }
  
  estdir <- direct(y = y,dom = dom,sweight = weights,domsize = domsize,data = datsample)
  estdir$Direct <- estdir$Direct*nsegcounty
  estdir$SD <- estdir$SD*nsegcounty
  
  datxy <- merge(estdir, urban11[,c("COUNTYFP","area_ts")], by.x = "Domain", by.y = "COUNTYFP")
  datxy$var <- datxy$SD^2
  
  estfh <- eblupFH(formula = Direct ~ area_ts, vardir = var,method = "REML",data = datxy)
  
  countyfips <- as.character(19000+ datxy$Domain)
  Cmatia <- Cmatia[countyfips, countyfips]
  
  ordermat <- getorder(Matrix(Cmatia))
  
  
  xus <- scale(datxy$area_ts)
  yus <- (datxy$Direct - meanx)/sdx
  sdus <- datxy$SD/sdx
  
  betamu0 <- cal_initialrx2(indexy = 1:nobs,y = yus,x = cbind(1,xus),K0 = 10,lam0 = 0.0001)
  
  lamv <- seq(0.3,1,by = 0.025)
  
  betamu01 <- betamu0
  j <- 0
  repeat{
    j <- j + 1
    resu1 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = rep(1,nobs*(nobs-1)/2),betam0 = betamu01,lam = lamv[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
    betamu01 <- resu1$betaest
    Khat <- length(unique(resu1$group))
    print(Khat)
    if(round(Khat)==3){break}
  }
  
  
  
  
  # lamv2 <- seq(0.3,1,by = 0.025)
  # alp <- c(0.1,0.5,1)
  # likv <- rep(0, length(alp))
  # 
  # beta_array <- array(0, dim = c(nobs, 2, length(alp)))
  # group_mat <- matrix(0, nobs, length(alp))
  # muhat_mat <- matrix(0, nobs, length(alp))
  # 
  # for(l in 1:length(alp))
  # {
  #   cvec1 <- exp(alp[l]*(1 - ordermat))
  #   betamu02 <- betamu0
  #   j <- 0
  #   repeat{
  #     j <- j + 1
  #     resu2 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = cvec1,betam0 = betamu02,lam = lamv2[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.01)
  #     betamu02 <- resu2$betaest
  #     Khat <- length(unique(resu2$group))
  #     #print(Khat)
  #     if(round(Khat)==3){break}
  #   }
  # 
  #   likv[l] <- resu2$loglikvalue
  # #
  # #
  #   beta_array[,,l] <- resu2$betaest
  #   group_mat[,l] <- resu2$group
  #   muhat_mat[,l] <- resu2$mhat
  # 
  # 
  # }


  # msemat[m,] <- c(mean((estdir$Direct - ypop)^2),
  #              mean((estfh$eblup - ypop)^2),
  #              mean((sdx*resu1$mhat + meanx - ypop)^2),
  #              mean((sdx*muhat_mat[,which.max(likv)] + meanx - ypop)^2))
  # 
  # rmsemat[m,] <- c(mean((estdir$Direct - ypop)^2/ypop^2),
  #                 mean((estfh$eblup - ypop)^2/ypop^2),
  #                 mean((sdx*resu1$mhat + meanx - ypop)^2/ypop^2),
  #                 mean((sdx*muhat_mat[,which.max(likv)] + meanx - ypop)^2/ypop^2))

 
   lamv2 <- seq(0.3,1,by = 0.025)
    cvec1 <- exp(0.5*(1 - ordermat))
    betamu02 <- betamu0
    j <- 0
    repeat{
      j <- j + 1
      resu2 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = cvec1,betam0 = betamu02,lam = lamv2[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
      betamu02 <- resu2$betaest
      Khat <- length(unique(resu2$group))
      #print(Khat)
      if(round(Khat)==3){break}
    }
  

  
  msemat[m,] <- c(mean((estdir$Direct - ypop)^2),
                  mean((estfh$eblup - ypop)^2),
                  mean((sdx*resu1$mhat + meanx - ypop)^2),
                  mean((sdx*resu2$mhat + meanx - ypop)^2))
  
  
 remsemat[m,] <- c(mean((estdir$Direct - ypop)^2/ypop^2),
               mean((estfh$eblup - ypop)^2/ypop^2),
               mean((sdx*resu1$mhat + meanx - ypop)^2/ypop^2),
               mean((sdx*resu2$mhat + meanx - ypop)^2/ypop^2))
 

        
 
  print(m)
}



boxplot(msemat[rowSums(msemat)>0,])
boxplot(remsemat[rowSums(msemat)>0,])

save(msemat,remsemat, file = "result/msemat50_diff2_r05_sd2.RData")

load("result/msemat50_diff15_r05_sd1.RData")
pdf("doc/figures/mse_diff15_r05_sd1.pdf",width = 6,height = 5)
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.1)))
dev.off()

load("result/msemat50_diff15_r1_sd1_v2.RData")
pdf("doc/figures/mse_diff15_r1_sd1.pdf",width = 6,height = 5)
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.1)))
dev.off()
load("result/msemat50_diff15_r2_sd1_v2.RData")
pdf("doc/figures/mse_diff15_r2_sd1.pdf",width = 6,height = 5)
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.1)))
dev.off()


load("result/msemat50_diff2_r1_sd1.RData")
pdf("doc/figures/mse_diff2_r1_sd1.pdf",width = 6,height = 5)
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.1)))
dev.off()

load("result/msemat50_diff2_r2_sd1.RData")
pdf("doc/figures/mse_diff2_r2_sd1.pdf",width = 6,height = 5)
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.1)))
dev.off()

load("result/msemat50_diff2_r05_sd1.RData")
pdf("doc/figures/mse_diff2_r05_sd1.pdf",width = 6,height = 5)
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.1)))
dev.off()


##### pictures togeter

### difference 1.5

pdf("doc/figures/mse_diff15_sd05.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat50_diff15_r05_sd05.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.05)))
load("result/msemat50_diff15_r1_sd05.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.05)))
load("result/msemat50_diff15_r2_sd05.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.05)))
dev.off()



pdf("doc/figures/mse_diff15_sd1.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat50_diff15_r05_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.1)))
load("result/msemat50_diff15_r1_sd1_v2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.1)))
load("result/msemat50_diff15_r2_sd1_v2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.1)))
dev.off()


pdf("doc/figures/mse_diff15_sd2.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat50_diff15_r05_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.2)))
load("result/msemat50_diff15_r1_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.2)))
load("result/msemat50_diff15_r2_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.2)))
dev.off()


#### diff2

pdf("doc/figures/mse_diff2_sd05.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat50_diff2_r05_sd05.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.05)))
load("result/msemat50_diff2_r1_sd05.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.05)))
load("result/msemat50_diff2_r2_sd05.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.05)))
dev.off()



pdf("doc/figures/mse_diff2_sd1.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat50_diff2_r05_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.1)))
load("result/msemat50_diff2_r1_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.1)))
load("result/msemat50_diff2_r2_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.1)))
dev.off()


pdf("doc/figures/mse_diff2_sd2.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat50_diff2_r05_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.05," and ",sigma,":",0.2)))
load("result/msemat50_diff2_r1_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.1," and ",sigma,":",0.2)))
load("result/msemat50_diff2_r2_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("rate:",0.2," and ",sigma,":",0.2)))
dev.off()



#### use lik to choose alp #####
#beta <- matrix(c(0.5,0.5,2.5,2.5,4.5,4.5),ncol=2,byrow = TRUE)
beta <- matrix(c(0.5,0.5,2,2,3.5,3.5),ncol=2,byrow = TRUE)
sdu <- 0.2
rate <- 0.2

msemat <- matrix(0,100,4)
remsemat <- matrix(0,100, 4)
colnames(msemat) <- colnames(remsemat) <- c("direct","regression","equal","sp")
for(m in 1:100)
{
  set.seed(631109 + m )
  #set.seed(631109)
  repeat{
    ypop <- sdx*(rowSums(xmat * beta[urban11$newgroup,]) + rnorm(nobs)*sdu) + meanx
    if(sum(ypop<=0)==0){break}
  }
  
  nsegurban <- round(ypop)  # number of segments have urban
  
  datsample <- as.data.frame(matrix(0, sum(nsampled), 3))
  colnames(datsample) <- c("y","dom","weights")
  datsample$dom <- rep(urban11$COUNTYFP, nsampled)
  
  for(i in 1:nobs)
  {
    popseg <- rep(0, nsegcounty[i])
    urbanindex <- sample.int( nsegcounty[i],nsegurban[i],replace = FALSE) # segemnts are urban
    popseg[urbanindex] <- 1
    
    sampled <- sample.int(nsegcounty[i], nsampled[i])
    
    datsample$y[datsample$dom == urban11$COUNTYFP[i]] <- popseg[sampled]
    
    datsample$weights[datsample$dom == urban11$COUNTYFP[i]] <- nsegcounty[i]/nsampled[i]
  }
  
  estdir <- direct(y = y,dom = dom,sweight = weights,domsize = domsize,data = datsample)
  estdir$Direct <- estdir$Direct*nsegcounty
  estdir$SD <- estdir$SD*nsegcounty
  
  datxy <- merge(estdir, urban11[,c("COUNTYFP","area_ts")], by.x = "Domain", by.y = "COUNTYFP")
  datxy$var <- datxy$SD^2
  
  estfh <- eblupFH(formula = Direct ~ area_ts, vardir = var,method = "REML",data = datxy)
  
  countyfips <- as.character(19000+ datxy$Domain)
  Cmatia <- Cmatia[countyfips, countyfips]
  
  ordermat <- getorder(Matrix(Cmatia))
  
  
  xus <- scale(datxy$area_ts)
  yus <- (datxy$Direct - meanx)/sdx
  sdus <- datxy$SD/sdx
  
  betamu0 <- cal_initialrx2(indexy = 1:nobs,y = yus,x = cbind(1,xus),K0 = 10,lam0 = 0.0001)
  
  lamv <- seq(0.3,1,by = 0.025)
  
  betamu01 <- betamu0
  j <- 0
  repeat{
    j <- j + 1
    resu1 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = rep(1,nobs*(nobs-1)/2),betam0 = betamu01,lam = lamv[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
    betamu01 <- resu1$betaest
    Khat <- length(unique(resu1$group))
    print(Khat)
    if(round(Khat)==3){break}
    if(round(Khat) < 3){break}
  }
  
  

  lamv2 <- seq(0.3,1.5,by = 0.025)
  alp <- c(0.1,0.5,1)
  likv <- rep(-10000, length(alp))

  beta_array <- array(0, dim = c(nobs, 2, length(alp)))
  group_mat <- matrix(0, nobs, length(alp))
  muhat_mat <- matrix(0, nobs, length(alp))
  
  Khatv <- rep(0, length(alp))

  for(l in 1:length(alp))
  {
    cvec1 <- exp(alp[l]*(1 - ordermat))
    betamu02 <- betamu0
    j <- 0
    repeat{
      j <- j + 1
      resu2 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = cvec1,betam0 = betamu02,lam = lamv2[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.01)
      betamu02 <- resu2$betaest
      Khatv[l] <- length(unique(resu2$group))
      #print(Khat)
      if(round(Khatv[l])==3){break}
      if(round(Khatv[l]) < 3){break}
    }
    #
    #
    if(Khatv[l]==3){likv[l] <- resu2$loglikvalue}

    beta_array[,,l] <- resu2$betaest
    group_mat[,l] <- resu2$group
    muhat_mat[,l] <- resu2$mhat


  }
  
  
  if(Khat ==3 & sum(Khatv ==3)>0)
  {
    msemat[m,] <- c(mean((estdir$Direct - ypop)^2),
                    mean((estfh$eblup - ypop)^2),
                    mean((sdx*resu1$mhat + meanx - ypop)^2),
                    mean((sdx*muhat_mat[,which.max(likv)] + meanx - ypop)^2))
    
    remsemat[m,] <- c(mean((estdir$Direct - ypop)^2/ypop^2),
                     mean((estfh$eblup - ypop)^2/ypop^2),
                     mean((sdx*resu1$mhat + meanx - ypop)^2/ypop^2),
                     mean((sdx*muhat_mat[,which.max(likv)] + meanx - ypop)^2/ypop^2))
  }

  
  
  print(m)
}

boxplot(msemat[rowSums(msemat)>0,])
save(msemat,remsemat, file = "result/msemat100_diff15_r2_sd2.RData")


pdf("doc/figures/mse100_4.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat100_diff15_r2_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","1.5 ","rate:",0.2," and ",sigma,":",0.1)))
load("result/msemat100_diff15_r2_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","1.5 ","rate:",0.2," and ",sigma,":",0.2)))

load("result/msemat100_diff2_r2_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","2 ","rate:",0.2," and ",sigma,":",0.1)))
load("result/msemat100_diff2_r2_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","2 ","rate:",0.2," and ",sigma,":",0.2)))
dev.off()


###### redefine of MSE of previous ####

beta <- matrix(c(0.5,0.5,2.5,2.5,4.5,4.5),ncol=2,byrow = TRUE)
#beta <- matrix(c(0.5,0.5,2,2,3.5,3.5),ncol=2,byrow = TRUE)
sdu <- 0.2
rate <- 0.2

sdx <- sd(urban11$area_ts)
meanx <- mean(urban11$area_ts)

nobs <- nrow(urban11)  
nsegcounty <- urban11$nseg
domsize <- urban11[,c("COUNTYFP","nseg")]
nsampled <- round(nsegcounty*rate)

xmat <- cbind(1, (urban11$area_ts - meanx)/sdx) 



simfun1 <- function(beta, sdu, rate, M=100)
{
  msemat <- matrix(0,nobs,4)
  colnames(msemat)  <- c("direct","regression","equal","sp")
  nrnzero <- 0 
  for(m in 1:M)
  {
    set.seed(631109 + m )
    #set.seed(631109)
    repeat{
      ypop <- sdx*(rowSums(xmat * beta[urban11$newgroup,]) + rnorm(nobs)*sdu) + meanx
      if(sum(ypop<=0)==0){break}
    }
    
    nsegurban <- round(ypop)  # number of segments have urban
    
    datsample <- as.data.frame(matrix(0, sum(nsampled), 3))
    colnames(datsample) <- c("y","dom","weights")
    datsample$dom <- rep(urban11$COUNTYFP, nsampled)
    
    for(i in 1:nobs)
    {
      popseg <- rep(0, nsegcounty[i])
      urbanindex <- sample.int( nsegcounty[i],nsegurban[i],replace = FALSE) # segemnts are urban
      popseg[urbanindex] <- 1
      
      sampled <- sample.int(nsegcounty[i], nsampled[i])
      
      datsample$y[datsample$dom == urban11$COUNTYFP[i]] <- popseg[sampled]
      
      datsample$weights[datsample$dom == urban11$COUNTYFP[i]] <- nsegcounty[i]/nsampled[i]
    }
    
    estdir <- direct(y = y,dom = dom,sweight = weights,domsize = domsize,data = datsample)
    estdir$Direct <- estdir$Direct*nsegcounty
    estdir$SD <- estdir$SD*nsegcounty
    
    datxy <- merge(estdir, urban11[,c("COUNTYFP","area_ts")], by.x = "Domain", by.y = "COUNTYFP")
    datxy$var <- datxy$SD^2
    
    estfh <- eblupFH(formula = Direct ~ area_ts, vardir = var,method = "REML",data = datxy)
    
    countyfips <- as.character(19000+ datxy$Domain)
    Cmatia <- Cmatia[countyfips, countyfips]
    
    ordermat <- getorder(Matrix(Cmatia))
    
    
    xus <- scale(datxy$area_ts)
    yus <- (datxy$Direct - meanx)/sdx
    sdus <- datxy$SD/sdx
    
    betamu0 <- cal_initialrx2(indexy = 1:nobs,y = yus,x = cbind(1,xus),K0 = 10,lam0 = 0.0001)
    
    lamv <- seq(0.3,1,by = 0.025)
    
    betamu01 <- betamu0
    j <- 0
    repeat{
      j <- j + 1
      resu1 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = rep(1,nobs*(nobs-1)/2),betam0 = betamu01,lam = lamv[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
      betamu01 <- resu1$betaest
      Khat <- length(unique(resu1$group))
      print(Khat)
      if(round(Khat)==3){break}
      if(round(Khat) < 3){break}
    }
    
    
    
    lamv2 <- seq(0.3,1.5,by = 0.025)
    alp <- c(0.1,0.5,1)
    likv <- rep(-10000, length(alp))
    
    beta_array <- array(0, dim = c(nobs, 2, length(alp)))
    group_mat <- matrix(0, nobs, length(alp))
    muhat_mat <- matrix(0, nobs, length(alp))
    
    Khatv <- rep(0, length(alp))
    
    for(l in 1:length(alp))
    {
      cvec1 <- exp(alp[l]*(1 - ordermat))
      betamu02 <- betamu0
      j <- 0
      repeat{
        j <- j + 1
        resu2 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde =  sdus, cvec = cvec1,betam0 = betamu02,lam = lamv2[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.01)
        betamu02 <- resu2$betaest
        Khatv[l] <- length(unique(resu2$group))
        #print(Khat)
        if(round(Khatv[l])==3){break}
        if(round(Khatv[l]) < 3){break}
      }
      #
      #
      if(Khatv[l]==3){likv[l] <- resu2$loglikvalue}
      
      beta_array[,,l] <- resu2$betaest
      group_mat[,l] <- resu2$group
      muhat_mat[,l] <- resu2$mhat
      
      
    }
    
    
    if(Khat ==3 & sum(Khatv ==3)>0)
    {
      nrnzero <- nrnzero + 1
      msemat <-  msemat + cbind((estdir$Direct - ypop)^2,
                                (estfh$eblup - ypop)^2,
                                (sdx*resu1$mhat + meanx - ypop)^2,
                                (sdx*muhat_mat[,which.max(likv)] + meanx - ypop)^2)
    }
    
    print(m)
  }

  msemat <- msemat/nrnzero
  return(list(msemat=msemat, nrnzero = nrnzero))
}


msemat100_diff2_r2_sd2 <- simfun1(beta, sdu = 0.2, rate = 0.2)
msemat100_diff2_r2_sd1 <- simfun1(beta, sdu = 0.1, rate = 0.2)
msemat100_diff2_r1_sd1 <- simfun1(beta, sdu = 0.1, rate = 0.1)
msemat100_diff2_r1_sd2 <- simfun1(beta, sdu = 0.2, rate = 0.1)

save(msemat100_diff2_r1_sd1,msemat100_diff2_r1_sd2, msemat100_diff2_r2_sd1,msemat100_diff2_r2_sd2,file = "result/msematnew.RData")

boxplot(sqrt(msemat100_diff2_r1_sd2$msemat))

boxplot(msemat[rowSums(msemat)>0,])
save(msemat,remsemat, file = "result/msemat100_diff1_r1_sd2.RData")

msemat1 <- msemat

msemat[,1]
msemat1[,1]


pdf("doc/figures/mse100_4.pdf",width = 8,height = 8)
par(mfrow = c(2,2))
load("result/msemat100_diff15_r2_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","1.5 ","rate:",0.2," and ",sigma,":",0.1)))
load("result/msemat100_diff15_r2_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","1.5 ","rate:",0.2," and ",sigma,":",0.2)))

load("result/msemat100_diff2_r2_sd1.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","2 ","rate:",0.2," and ",sigma,":",0.1)))
load("result/msemat100_diff2_r2_sd2.RData")
boxplot(msemat[rowSums(msemat)>0,],main = expression(paste("diff:","2 ","rate:",0.2," and ",sigma,":",0.2)))
dev.off()
