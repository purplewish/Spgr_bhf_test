##### model 3
source("Research/survey/code/simdat.R")
IAarea <- read.csv("Research/survey/IAarea.csv",stringsAsFactors = FALSE)
library(sae)
library(glmnet)
Nv = seq(50,300,by = 1)
srh <- 0.3
srl <- 0.1
Ax <- IAarea$area
ax <- 2
bx <- 20
beta <- c(1,1)
N <- 99

M <- 100


### for model 

rmse_mat3 <- matrix(0, M, 6)
colnames(rmse_mat3) <- c("direct","eblup","reg","ridge","lasso","Alasso")
numcoef3 <- rep(0,M)

for(m in 1: M)
{
  sim1 <- simdat(N = N,Nv = Nv,Ax = Ax,ax = ax,bx = bx,model = 3,beta = c(100,1),
                 sde = 25,sdv = 25,sdb0 = 10, sdb1 = 0.2,seed = 2330 + m)
  dats <- sim1$sample
  datp <- sim1$pop
  popn <- sim1$domsize
  meany <- ddply(datp,.(dom),summarize, meany = mean(y))
  
  set.seed(m + 1605)
  estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
  meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
  
  estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx,popnsize = popn, method = "REML",data = dats)
  
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
  
  
  ## regression 
  
  res0 <- lm(y~x, weights = dats$weights, data = dats)
  betaest0 <- coef(res0)
  estreg <- cbind(1,meanx$meanx) %*% betaest0
  
  
  ## ridge regression 
  cv1 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights, alpha = 0, standardize = TRUE)
  betaest1 <- coef(cv1,s = cv1$lambda.min)
  estridge <- as.matrix(xmp %*% betaest1)[,1]
  
  
  ## lasso
  cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
  betaest2 <- coef(cv2,s = cv2$lambda.min)
  estlasso <- as.matrix(xmp %*% betaest2)[,1]
  numcoef3[m] <- sum(betaest2!=0)
  
  ## adaptive lasso 
  ww <- 1/abs(betaest1[-1])
  cv3 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1,
                   penalty.factor = ww)
  betaest3 <- coef(cv3,s = cv3$lambda.min)
  estalasso <- as.matrix(xmp%*% betaest3)[,1]
  
  
  rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                 mean((meany[,2] - estbhf$eblup$eblup)^2),
                 mean((meany[,2] - estreg)^2),
                 mean((meany[,2] - estridge)^2),
                 mean((meany[,2] - estlasso)^2),
                 mean((meany[,2] - estalasso)^2)))
  rmse_mat3[m,] <- rmse
  
}

save(rmse_mat3, numcoef3, file = "Research/survey/result/rmse_mat3.RData")


boxplot(rmse_mat3)




