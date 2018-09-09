###### simulation to two level model, check code ###
# without spatial group structure  only x included different intercept and different slope 
setwd("Research/NRI_urban/")
source("code/mcp_scad.R")
n <- 200
mux <- 1
sdx <- 0.5
sdu <- 0.5
sde <- rnorm(n)*0.1 + 0.5
mu0 <- 5

beta <- matrix(c(1,4),ncol = 2) ## two groups, number of columns

ng <- ncol(beta)

set.seed(231315)
group <- sample.int(ng,n,replace = TRUE)
x <- rnorm(n) * sdx + mux
y <- 5 + x*beta[1,group] + rnorm(n)*sdu
yhat <- y + rnorm(n)*sde

Omega0 <-  1/(sde^2 + sdu^2)
dat  <- data.frame(y = y, x = x, group = group, weights = Omega0)

#x <- cbind(1,x)

#### if sig2 and group information are known

lm(y~x, data = dat[dat$group==1,],weights = weights)
lm(y~x, data = dat[dat$group==2,],weights = weights)




#### algorithm ####
source("code/Spgr_sampling2.R")
source("code/Spgr_sampling_em2.R")
library(Spgr)
z0 <- matrix(1, nrow = n, ncol=1)
xmat <- matrix(x, nrow = n, ncol=1)
betam0 <- cal_initialr2(indexy = 1:n,y = yhat,z = z0,x = matrix(x,nrow=n, ncol=1),K0 = 20,lam0 = 0.0001)

### BIC is defined as -2loglikevalue + Cn logn (Kp)
cvec <- rep(1, n*(n-1)/2)
res1 <- Spgr_sampling2(yhat = yhat, z= z0, x = xmat, sde = sde, cvec = cvec,betam0 = betam0,
                      lam = 0.6,tolabs = 1e-4, tolrel = 0.05)

res2 <- Spgr_sampling_em2(yhat = yhat, z= z0,x = xmat, sde = sde, cvec = cvec,betam0 = betam0,
                         lam = 0.6,tolabs = 1e-4, tolrel = 0.05)

betam01 <- betam0
lam1 <- seq(0.2,2,by = 0.05)
bicv1 <- rep(0, length(lam1))
beta_array1 <- array(0, dim = c(n,1,length(lam1)))
groupmat1 <- matrix(0, n, length(lam1))
loglik1 <- rep(0, length(lam1))

for(j in 1:length(lam1))
{
  resj <- Spgr_sampling2(yhat = yhat, z= z0, x = xmat, sde = sde, cvec = cvec,betam0 = betam0,
                         lam = lam1[j],tolabs = 1e-4, tolrel = 0.05)
  betam01 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicv1[j] <- -2*resj$loglikvalue + 2*log(n + 1)*log(n)*(Khat + 1)
  beta_array1[,,j] <- resj$betaest
  groupmat1[,j] <- resj$group
  loglik1[j] <- resj$loglikvalue
}

which.min(bicv1)

plot(beta_array1[,,9])



betam02 <- betam0
lam2 <- seq(0.2,2,by = 0.05)
bicv2 <- rep(0, length(lam2))
beta_array2 <- array(0, dim = c(n,2,length(lam2)))
groupmat2 <- matrix(0, n, length(lam2))

for(j in 1:length(lam2))
{
  resj <- Spgr_sampling_em2(yhat = yhat, z= z0,x = xmat, sde = sde, cvec = cvec,betam0 = betam0,
                            lam = lam2[j],tolabs = 1e-4, tolrel = 0.05)
  betam02 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicv2[j] <- -2*resj$loglikvalue + 2*log(n + 1)*log(n)*(Khat + 1)
  beta_array2[,,j] <- resj$betaest
  groupmat2[,j] <- resj$group
}

plot(bicv2)
plot(beta_array2[,,9])


############real data #######

daturban <- read.csv("data/urbania.csv",stringsAsFactors = FALSE)
daturban$NAME <- tolower(daturban$NAME)
daturban[daturban$NAME == "o'brien","NAME"]<- "obrien"
dat12 <- daturban[daturban$year==2012,]
dat12 <- dat12[dat12$var_nri!=0,] ## remove zeros

n12 <- nrow(dat12)

plot(log(dat12$area_nri) - lm(log(area_nri)~log(area_ts), data = dat12)$fitted.values)

sdy <- sd(dat12$area_nri)
sd12 <- sqrt(dat12$var_nri)
sd12 <- sd12/sdy

betam0u <- cal_initialr2(indexy = 1:n12,y = scale(dat12$area_nri),z = matrix(1,n12),x = matrix(scale(dat12$area_ts)),K0 = 15,lam0 = 0.0001)


res_u1 <- Spgr_sampling2(yhat = scale(dat12$area_nri), z=matrix(1,n12), x = matrix(scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betam0u, lam = 0.2,maxiter = 2000,tolabs = 1e-4,tolrel = 0.01)

res_u2 <- Spgr_sampling_em2(yhat = scale(dat12$area_nri), z=matrix(1,n12), x = matrix(scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betam0u, lam = 0.2,maxiter = 2000,tolabs = 1e-4,tolrel = 0.05)

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
  bicu2[j] <- -2*resj$loglikvalue + log(log(2*n12))*log(n12)*(Khat*2)
  beta_arrayu2[,,j] <- resj$betaest
  groupmatu2[,j] <- resj$group
  sig2u[j] <- resj$sig2
}
plot(bicu2)
plot(beta_arrayu2[,,which.min(bicu2)])
mapfun(group = groupmatu2[,which.min(bicu2)],countyname = dat12$NAME)

plot(dat12$area_ts,dat12$area_nri, col = groupmatu2[,18])


load("/Volumes/cssm_groups$/Grads/XinWang/county_adjacency/data/county_adjacency.RData")
Cmatia <- adjMat[["19"]]
countyfips <- as.character(19000+ dat12$COUNTYFP)
Cmatia <- Cmatia[countyfips, countyfips]


ordermat <- getorder(Matrix(Cmatia))
cvec1 <- exp(0.8*(1 - ordermat))
res_u3 <- Spgr_sampling_em(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = cvec1,betam0 =  betam0u, lam = 0.12,maxiter = 1000,tolabs = 1e-4,tolrel = 0.05)

mapfun(group = res_u3$group,countyname = dat12$NAME)


betamu3 <- betam0u
lamu3 <- seq(0.05,0.2,by = 0.005)
bicu3 <- rep(0, length(lamu3))
beta_arrayu3 <- array(0, dim = c(n12,2,length(lamu3)))
groupmatu3 <- matrix(0, n12, length(lamu3))
sig2u3 <- rep(0, length(lamu3))

for(j in 1:length(lamu3))
{
  resj <- Spgr_sampling_em(yhat = scale(dat12$area_nri), x = cbind(1,scale(dat12$area_ts)),sde = sd12,cvec = rep(1,n12*(n12-1)/2),betam0 =  betamu2, lam = lamu3[j],maxiter = 1000,tolabs = 1e-4,tolrel = 0.05)
  betamu3 <- resj$betaest
  Khat <- length(unique(resj$group))
  bicu3[j] <- -2*resj$loglikvalue + 2*log(log(2*n12))*log(n12)*(Khat*2)
  beta_arrayu3[,,j] <- resj$betaest
  groupmatu3[,j] <- resj$group
  sig2u3[j] <- resj$sig2
}
plot(bicu3)
plot(beta_arrayu3[,,17])
mapfun(group = groupmatu3[,17],countyname = dat12$NAME)

