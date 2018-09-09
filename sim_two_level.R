###### simulation to two level model, check code ###
# without spatial group structure  only x included different intercept and different slope 
setwd("Research/NRI_urban/")
source("code/mcp_scad.R")
source("code/mapfun.R")
n <- 200
mux <- 1
sdx <- 0.5
sdu <- 0.5
sde <- rnorm(n)*0.1 + 0.5

beta <- matrix(c(1,1,3,3),ncol = 2) ## two groups, number of columns

ng <- ncol(beta)

set.seed(231315)
group <- sample.int(ng,n,replace = TRUE)
x <- rnorm(n) * sdx + mux
y <- beta[1,group] + x*beta[2,group] + rnorm(n)*sdu
yhat <- y + rnorm(n)*sde

Omega0 <-  1/(sde^2 + sdu^2)
dat  <- data.frame(y = y, x = x, group = group, weights = Omega0)

#x <- cbind(1,x)

#### if sig2 and group information are known

lm(y~x, data = dat[dat$group==1,],weights = weights)
lm(y~x, data = dat[dat$group==2,],weights = weights)




#### algorithm ####
source("code/Spgr_sampling.R")
source("code/Spgr_sampling_em.R")
source("code/Spgr_sampling_emp.R")
library(Spgr)
betam0 <- cal_initialx(y = yhat,x = cbind(1,x),lam0 = 0.0001)

### BIC is defined as -2loglikevalue + Cn logn (Kp)
cvec <- rep(1, n*(n-1)/2)
res1 <- Spgr_sampling(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam0,
                      lam = 0.7,tolabs = 1e-4, tolrel = 0.05)

res2 <- Spgr_sampling_em(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam0,
                      lam = 2,tolabs = 1e-4, tolrel = 0.05)

res21 <- Spgr_sampling_em(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam0,
                         lam = 0.5,tolabs = 1e-4, tolrel = 0.05)

res3 <- Spgr_sampling_emp(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam0, 
                          phi =1,
                         lam = 0.5,tolabs = 1e-4, tolrel = 0.05)

res31 <- Spgr_sampling_emp(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam0, 
                          phi =1,
                          lam = 2,tolabs = 1e-4, tolrel = 0.05)


mean(((yhat - y)^2)/sde^2)
mean(((yhat - res2$mhat)^2 + res2$Vhat)/sde^2)
mean(((yhat - res21$mhat)^2 + res21$Vhat)/sde^2)

plot(yhat,res2$mhat)
plot(yhat,res21$mhat)

plot(y, res2$mhat)
abline(0,1)
plot(y, res21$mhat)
abline(0,1)

mean(((yhat - res2$mhat)^2 + res2$Vhat)/sde^2)

betam01 <- betam0
lam1 <- seq(0.2,2,by = 0.05)
bicv1 <- rep(0, length(lam1))
beta_array1 <- array(0, dim = c(n,2,length(lam1)))
groupmat1 <- matrix(0, n, length(lam1))

for(j in 1:length(lam1))
{
  resj <- Spgr_sampling(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam01,lam = lam1[j],maxiter = 1000, tolabs = 1e-4, tolrel = 0.05)
  betam01 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicv1[j] <- -2*resj$loglikvalue + log(n*2)*log(n)*(Khat*2)
  beta_array1[,,j] <- resj$betaest
  groupmat1[,j] <- resj$group
}




betam02 <- betam0
lam2 <- seq(0.2,2,by = 0.05)
bicv2 <- rep(0, length(lam2))
beta_array2 <- array(0, dim = c(n,2,length(lam2)))
groupmat2 <- matrix(0, n, length(lam2))
phi2 <- khat2<- rep(0, length(lam2))


for(j in 1:length(lam2))
{
  resj <- Spgr_sampling_em(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam02,lam = lam2[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
  betam02 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicv2[j] <- -2*resj$loglikvalue + log(n*2)*log(n)*(Khat*2)
  beta_array2[,,j] <- resj$betaest
  groupmat2[,j] <- resj$group
  phi2[j] <- mean(((yhat - resj$mhat)^2 + resj$Vhat)/sde^2)
  khat2[j] <- Khat
}

plot(khat2, phi2, xlab = expression(hat(K)),ylab = expression(hat(phi)))

plot(beta_array2[,,13])

betam021 <- betam0
lam21 <- seq(0.2,2,by = 0.05)
bicv21 <- rep(0, length(lam2))
beta_array21 <- array(0, dim = c(n,2,length(lam2)))
groupmat21 <- matrix(0, n, length(lam2))
phi21 <- khat21 <- rep(0, length(lam2))

for(j in 1:length(lam2))
{
  resj <- Spgr_sampling_emp(yhat = yhat, x = cbind(1,x), sde = sde, cvec = cvec,betam0 = betam021,lam = lam2[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
  betam021 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicv21[j] <- -2*resj$loglikvalue + log(n*2)*log(n)*(Khat*2)
  beta_array21[,,j] <- resj$betaest
  groupmat21[,j] <- resj$group
  phi21[j] <- resj$phi
  khat21[j]<- Khat
}

plot(khat21, phi21,xlab = expression(hat(K)),ylab = expression(hat(phi)))
plot(beta_array2[,,5])



##### cross validation #####
library(cvTools)
set.seed(1048)
cvgroup <- cvFolds(n = 200,K = 10)


lam22 <- seq(0.2,1.5,by = 0.1)
mse22 <- rep(0, length(lam22))

for(j in 1:length(lam22))
{
  sse <- 0
  for(l in 1:10)
  {
    trainind <- cvgroup$subsets[cvgroup$which!=l]
    testind <- cvgroup$subsets[cvgroup$which==l]
    ntrain <- length(trainind)
    betam0 <- cal_initialx(y = yhat[trainind],x = cbind(1,x[trainind]),lam0 = 0.0001)
    
    resl <- Spgr_sampling_em(yhat = yhat[trainind], x = cbind(1,x[trainind]), sde = sde[trainind], cvec = rep(1,ntrain*(ntrain-1)/2),betam0 = betam0,lam = lam22[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
    
    groupnuml <- as.numeric(table(resl$group))
    ngestl <- length(groupnuml)
    betaestl <- unique(resl$betaest)
    
    xtest <- cbind(1, x[testind])
    
    ypredl <- (xtest %*% t(betaestl))%*%(groupnuml/sum(groupnuml))
    
    sse = sse + sum((ypredl - yhat[testind])^2)
  }
  mse22[j] <- sse/n
}

betam0 <- cal_initialx(y = yhat,x = cbind(1,x),lam0 = 0.0001)
res22 <- Spgr_sampling_em(yhat = yhat, x = cbind(1,x), sde = sde, 
                          cvec = rep(1,n*(n-1)/2),betam0 = betam0,
                         lam = lam22[which.min(mse22)],tolabs = 1e-4, tolrel = 0.05)



############## cross validation for 100 times #####
nrep <- 100
group_rep <- matrix(0, n, nrep)
msebeta <- rep(0, n , rep)
beta_array_nrep <- array(0, dim = c(n,2, nrep))

for(mm in 1:nrep)
{
  set.seed(231315 + mm)
  group <- sample.int(ng,n,replace = TRUE)
  x <- rnorm(n) * sdx + mux
  y <- beta[1,group] + x*beta[2,group] + rnorm(n)*sdu
  yhat <- y + rnorm(n)*sde
  
  Omega0 <-  1/(sde^2 + sdu^2)
  dat  <- data.frame(y = y, x = x, group = group, weights = Omega0)

  cvgroup <- cvFolds(n = 200,K = 10)
  
  lam22 <- seq(0.2,1.5,by = 0.1)
  mse22 <- rep(0, length(lam22))
  
  for(j in 1:length(lam22))
  {
    sse <- 0
    for(l in 1:10)
    {
      trainind <- cvgroup$subsets[cvgroup$which!=l]
      testind <- cvgroup$subsets[cvgroup$which==l]
      ntrain <- length(trainind)
      betam0 <- cal_initialx(y = yhat[trainind],x = cbind(1,x[trainind]),lam0 = 0.0001)
      
      resl <- Spgr_sampling_em(yhat = yhat[trainind], x = cbind(1,x[trainind]), sde = sde[trainind], cvec = rep(1,ntrain*(ntrain-1)/2),betam0 = betam0,lam = lam22[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
      
      groupnuml <- as.numeric(table(resl$group))
      ngestl <- length(groupnuml)
      betaestl <- unique(resl$betaest)
      
      xtest <- cbind(1, x[testind])
      
      ypredl <- (xtest %*% t(betaestl))%*%(groupnuml/sum(groupnuml))
      
      sse = sse + sum((ypredl - yhat[testind])^2)
    }
    mse22[j] <- sse/n
  }
  
  betam0 <- cal_initialx(y = yhat,x = cbind(1,x),lam0 = 0.0001)
  res22 <- Spgr_sampling_em(yhat = yhat, x = cbind(1,x), sde = sde, 
                            cvec = rep(1,n*(n-1)/2),betam0 = betam0,
                            lam = lam22[which.min(mse22)],tolabs = 1e-4, tolrel = 0.05)
  group_rep[,mm] <- res22$group
  
  msebeta[mm] <- mean(rowSums((res22$betaest - t(beta[,group]))^2))
  beta_array_nrep[,,mm] <- res22$betaest
  print(mm)
}


save(group_rep, msebeta, beta_array_nrep, file = "result/simcv.RData")

sum(apply(group_rep[,1:40],2,function(x){length(unique(x))})==2)



#################real data ###################

daturban <- read.csv("data/urbania.csv",stringsAsFactors = FALSE)
daturban$NAME <- tolower(daturban$NAME)
daturban[daturban$NAME == "o'brien","NAME"]<- "obrien"
dat12 <- daturban[daturban$year==2012,]
dat12 <- dat12[dat12$var_nri!=0,] ## remove zeros

n12 <- nrow(dat12)

plot(log(dat12$area_nri) - lm(log(area_nri)~log(area_ts), data = dat12)$fitted.values)

sdy <- sd(dat12$area_nri)
sd12 <- sqrt(dat12$var_nri)
sdx <- sd(dat12$area_ts)
sd12 <- sd12/sdy

betam0u <- cal_initialrx2(indexy = 1:n12,y = dat12$area_nri/sdx,x = cbind(1,scale(dat12$area_ts)),K0 = 15,lam0 = 0.0001)


res_u1 <- Spgr_sampling(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betam0u, lam = 0.1,maxiter = 2000,tolabs = 1e-4,tolrel = 0.05)

res_u2 <- Spgr_sampling_em(yhat = dat12$area_nri/sdx, x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betam0u, lam = 0.2,maxiter = 2000,tolabs = 1e-4,tolrel = 0.05)

res_u22 <- Spgr_sampling_em(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sqrt(dat12$var_nri)/sdy,cvec = rep(1,n12*(n12-1)/2),betam0 =  betam0u, lam = 0.2,maxiter = 2000,tolabs = 1e-4,tolrel = 0.05)

plot(dat12$area_ts,dat12$area_nri, col = res_u2$group)


betamu1 <- betam0u
lamu1 <- seq(1,2,by = 0.01)
bicu1 <- rep(0, length(lamu1))
beta_arrayu1 <- array(0, dim = c(n12,2,length(lamu1)))
groupmatu1 <- matrix(0, n12, length(lamu1))

for(j in 1:length(lam2))
{
  resj <- Spgr_sampling(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betamu1, lam = lamu1[j],maxiter = 2000,tolabs = 1e-4,tolrel = 0.05)
  betamu1 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicu1[j] <- -2*resj$loglikvalue + log(n12*2)*log(n12)*(Khat*2)
  beta_arrayu1[,,j] <- resj$betaest
  groupmatu1[,j] <- resj$group
}


betamu2 <- betam0u
lamu2 <- seq(0.05,0.2,by = 0.005)
bicu2 <- rep(0, length(lamu2))
beta_arrayu2 <- array(0, dim = c(n12,2,length(lamu2)))
groupmatu2 <- matrix(0, n12, length(lamu2))
sig2u <- rep(0, length(lamu2))

for(j in 1:length(lamu2))
{
  resj <- Spgr_sampling_em(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betamu2, lam = lamu2[j],maxiter = 1000,tolabs = 1e-4,tolrel = 0.05)
  betamu2 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicu2[j] <- -2*resj$loglikvalue + 1.1*log(log(2*n12))*log(n12)*(Khat*2)
  beta_arrayu2[,,j] <- resj$betaest
  groupmatu2[,j] <- resj$group
  sig2u[j] <- resj$sig2
}
plot(bicu2)
plot(beta_arrayu2[,,which.min(bicu2)])
mapfun(group = groupmatu2[,which.min(bicu2)],countyname = dat12$NAME)

plot(dat12$area_ts,dat12$area_nri, col = groupmatu2[,18])
table(groupmatu2[,which.min(bicu2)])


load("/Volumes/cssm_groups$/Grads/XinWang/county_adjacency/data/county_adjacency.RData")
Cmatia <- adjMat[["19"]]
countyfips <- as.character(19000+ dat12$COUNTYFP)
Cmatia <- Cmatia[countyfips, countyfips]
save(adjMat, file = "data/adjMat.RData")

load("data/adjMat.RData")

ordermat <- getorder(Matrix(Cmatia))
cvec1 <- exp((1 - ordermat))
res_u3 <- Spgr_sampling_em(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = cvec1,betam0 =  betam0u, lam = 0.12,maxiter = 1000,tolabs = 1e-4,tolrel = 0.05)

mapfun(group = res_u3$group,countyname = dat12$NAME)


betamu3 <- betam0u
lamu3 <- seq(0.05,0.2,by = 0.005)
bicu3 <- rep(0, length(lamu3))
beta_arrayu3 <- array(0, dim = c(n12,2,length(lamu3)))
groupmatu3 <- matrix(0, n12, length(lamu3))
sig2u3 <- rep(0, length(lamu3))

cvec1 <- exp((1 - ordermat))

for(j in 1:length(lamu3))
{
  resj <- Spgr_sampling_em(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = cvec1,betam0 =  betamu3, lam = lamu3[j],maxiter = 1000,tolabs = 1e-4,tolrel = 0.05)
  betamu3 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicu3[j] <- -2*resj$loglikvalue + log(2*n12)*log(n12)*(Khat*2)
  beta_arrayu3[,,j] <- resj$betaest
  groupmatu3[,j] <- resj$group
  sig2u3[j] <- resj$sig2
}
plot(bicu3)
plot(beta_arrayu3[,,which.min(bicu3)])
mapfun(group = groupmatu3[,which.min((bicu3))],countyname = dat12$NAME)


#### cross validation #####
library(cvTools)
set.seed(1048)
cvgroupu <- cvFolds(n = n12,K = 5)

lamu22 <- seq(0.05,0.2,by = 0.005)
mseu22 <- rep(0, length(lamu22))

yus <- scale(dat12$area_nri)
xus <- scale(dat12$area_ts)

for(j in 1:length(lamu22))
{
  sse <- 0
  for(l in 1:5)
  {
    trainind <- cvgroupu$subsets[cvgroupu$which!=l]
    testind <- cvgroupu$subsets[cvgroupu$which==l]
    ntrain <- length(trainind)
    betamu0 <- cal_initialrx2(indexy = 1:ntrain,y = yus[trainind],x = cbind(1,xus[trainind]),K0 = 10,lam0 = 0.0001)
  
    resl <- Spgr_sampling_em(yhat = yus[trainind], x = cbind(1,xus[trainind]), sde = sd12[trainind], cvec = rep(1,ntrain*(ntrain-1)/2),betam0 = betamu0,lam = lamu22[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
    
    groupnuml <- as.numeric(table(resl$group))
    ngestl <- length(groupnuml)
    betaestl <- unique(resl$betaest)
    
    xtest <- cbind(1, xus[testind])
    
    ypredl <- (xtest %*% t(betaestl))%*%(groupnuml/sum(groupnuml))
    
    sse = sse + sum((ypredl - yus[testind])^2)
  }
  mseu22[j] <- sse/n12
}

betamu0 <- cal_initialrx2(indexy = 1:n12,y = yus,x = cbind(1,xus),K0 = 10,lam0 = 0.0001)

resu22 <- Spgr_sampling_em(yhat = yus, x = cbind(1,xus), sde = sd12, cvec = rep(1,n12*(n12-1)/2),betam0 = betamu0,lam = lamu22[9],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)

plot(resu22$betaest)

plot(dat12$area_ts, dat12$area_nri, col = resu22$group,xlab = "ts",ylab ="nri")
qplot(dat12$area_ts,dat12$area_nri,color = as.factor(resu22$group)) 

save(mseu22, file ="data/msecv_urban2012.RData")

mapfun(group = resu22$group,countyname = dat12$NAME)

load("data/msecv_urban2012.RData")
### refit ###

unique(resu22$betaest)

resufit <- Spgr_sampling_em_fixed(yhat = dat12$area_nri, x = cbind(1,dat12$area_ts),sde = sqrt(dat12$var_nri),groupest = resu22$group,betav0 = c(t(unique(resu22$betaest))))

print(xtable(resufit$betaest))

sqrt(diag(resufit$covm))

pdf("doc/figures/NRIvsmodel.pdf",width = 6,height = 5)
qplot(dat12$area_nri,resufit$mhat) + geom_abline(slope = 1,intercept = 0) + theme_bw()+
  xlab("NRI") + ylab("model based")
dev.off()




resufit$sig2

#### take log first ?? 



#### spatial  

lamu22 <- seq(0.05,0.2,by = 0.005)
alp <- c(0.5,1,2)
mseu23 <- matrix(0,length(lamu22), length(alp))


for(l in 1:length(alp))
{
  cvec1 <- exp(alp[l]*(1 - ordermat))
  
  for(j in 1:length(lamu22))
  {
    sse <- 0
    for(l in 1:5)
    {
      trainind <- cvgroupu$subsets[cvgroupu$which!=l]
      testind <- cvgroupu$subsets[cvgroupu$which==l]
      ntrain <- length(trainind)
      betamu0 <- cal_initialrx2(indexy = 1:ntrain,y = yus[trainind],x = cbind(1,xus[trainind]),K0 = 10,lam0 = 0.0001)
      
      resl <- Spgr_sampling_em(yhat = yus[trainind], x = cbind(1,xus[trainind]), sde = sd12[trainind], cvec = cvec1,betam0 = betamu0,lam = lamu22[j],maxiter = 2000, tolabs = 1e-4, tolrel = 0.05)
      
      groupnuml <- as.numeric(table(resl$group))
      ngestl <- length(groupnuml)
      betaestl <- unique(resl$betaest)
      
      xtest <- cbind(1, xus[testind])
      
      ypredl <- (xtest %*% t(betaestl))%*%(groupnuml/sum(groupnuml))
      
      sse = sse + sum((ypredl - yus[testind])^2)
    }
    mseu23[j,l] <- sse/n12
  }
}






