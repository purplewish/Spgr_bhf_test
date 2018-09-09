source("Research/survey/code/simdat.R")
library(plyr)
IAarea <- read.csv("Research/survey/data/IAarea.csv",stringsAsFactors = FALSE)
library(sae)
library(glmnet)
Nv = seq(50,300,by = 1)
srh <- 0.3
srl <- 0.1
Ax <- IAarea$area
ax <- 2
bx <- 20
beta <- c(100,1)
N <- 99
M <- 100
sde <- 25



##### test lasso weights ####
sim1 <- simdat(N = N,Nv = Nv,srh = srh, srl = srl,
               Ax = Ax,ax = ax,bx = bx,model = 1,beta = beta,
               sde = sde, sdv = sdv, sdb0 = sdb0, sdb1 = sdb1, seed = 2330)

dats <- sim1$sample
datp <- sim1$pop
popn <- sim1$domsize

meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))

xmp <- matrix(0, nrow(meanx), 2*N)
xmp[,1] <- 1
xmp[,2] <- meanx$meanx
xmp[1,seq(3,2*N,by = 2)] <- -1
xmp[1,seq(4,2*N,by = 2)] <- - meanx$meanx[1]

xm <- matrix(0, nrow(dats), 2*(N))
xm[,1] <- 1
xm[,2] <- dats$x
xm[1:popn$sampsize[1],seq(3,ncol(xm),by = 2)] <- -1
xm[1:popn$sampsize[1],seq(4,ncol(xm),by = 2)] <- -dats$x[1:popn$sampsize[1]]

for(i in 2:N)
{
  index1 <- sum(popn$sampsize[1:(i-1)]) + 1
  index2 <- sum(popn$sampsize[1:i])
  xm[index1:index2,(i-1)*2 +1] <- 1
  xm[index1:index2,2*i] <- dats$x[dats$dom==i]
  
  xmp[i, (i-1)*2 +1]  <- 1
  xmp[i,2*i]  <- meanx$meanx[i]
}


#cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
cv2 <- cv.glmnet(x = xm[,-1]*sqrt(dats$weights),y = dats$y*sqrt(dats$weights),alpha = 1)
betaest2 <- coef(cv2,s = cv2$lambda.min)

cv21 <- cv.glmnet(x = xm[,-1],y = dats$y, weights = dats$weights, alpha = 1)
betaest21 <- coef(cv21,s = cv21$lambda.min)
estlasso <- as.matrix(xmp %*% betaest2)[,1]



### model based 
res_m1 <- sim_model(M = 100,N = N,srh = 0.3,srl = 0.1,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 1,beta = beta,sde = sde)

res_m2 <- sim_model(M = 100,N = N,srh = 0.3,srl = 0.1,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 2,beta = beta,sde = sde, sdv = 25)

res_m3 <- sim_model(M = 100,N = N,srh = 0.3,srl = 0.1,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 3,beta = beta,sde = sde,sdv = 25,
                    sdb0 = 10, sdb1 = 0.2)

save(res_m1, res_m2, res_m3, file ="Research/survey/result/res_model.RData")


load("/Users/Xin/Research/survey/result/res_model.RData")

boxplot(res_m1$rmse[,-c(1,4)])
boxplot(res_m2$rmse[,-c(1)])
boxplot(res_m3$rmse[,-c(1,3)])

##### lower sampling rate 
res_m11 <- sim_model(M = 100,N = N,srh = 0.2,srl = 0.05,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 1,beta = beta,sde = sde)

res_m21 <- sim_model(M = 100,N = N,srh = 0.2,srl = 0.05,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 2,beta = beta,sde = sde, sdv = 25)

res_m31 <- sim_model(M = 100,N = N,srh = 0.2,srl = 0.05,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 3,beta = beta,sde = sde,sdv = 25,
                    sdb0 = 10, sdb1 = 0.2)

save(res_m11, res_m21, res_m31, file ="Research/survey/result/res_model1.RData")
load("Research/survey/result/res_model1.RData")


### design based 

res_d1 <- sim_design(M = 100,N = N,srh = 0.3,srl = 0.1,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 1,beta = beta,sde = sde)

res_d2 <- sim_design(M = 100,N = N,srh = 0.3,srl = 0.1,Nv = Nv,Ax = Ax,
                    ax = ax,bx = bx,model = 2,beta = beta,sde = sde, sdv = 25)

res_d3 <- sim_design(M = 100,N = N,srh = 0.3,srl = 0.1,Nv = Nv,Ax = Ax,
                     ax = ax,bx = bx,model = 3,beta = beta,sde = sde, sdv = 25,
                     sdb0 = 10, sdb1 = 0.2)

save(res_d1, res_d2, res_d3, file = "Research/survey/result/res_design.RData")

