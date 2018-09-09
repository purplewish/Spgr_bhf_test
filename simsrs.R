# M number of simulation
library(sae)
library(glmnet)
simsrs <- function(M, dat, popn, sr = 0.03, seed = 451920)
{
  
  numcoef1 <- rep(0, M)
  rmse_mat1 <- matrix(0, M, 6)
  colnames(rmse_mat1) <- c("direct","eblup","reg","ridge","lasso","Alasso")
  
  Nv <- popn$numseg
  
  N <- length(Nv) ## total number of counties
  nv <- round(Nv * sr) ## sample size 
  n0 <- sum(nv) ## total sample size
  
  popsize <- popn[,c("county","numseg")]
  meanx <- popn[,c("county","mean01")]
  meany <- popn$mean11
  
  
  dats <- as.data.frame(matrix(0, n0, 4)) # sampled data
  colnames(dats) <- c("dom", "y", "x", "weights")
  dats$dom <- rep(1:N, nv)
  
  
  for(m in 1:M)
  {
    set.seed(seed + m)
    for(i in 1:N)
    {
      indpi <- dat$county == i 
      ypi <- dat$A2011[indpi] ## population y
      xpi <- dat$A2001[indpi] ## population x
      
      indi <- dats$dom == i
      
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
      dats$weights[indi] <- Nv[i]/nv[i]
    }
    
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popsize)
    
    
    estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popsize, 
                       method = "REML",data = dats)
    
    
    res0 <- lm(y~x, weights = dats$weights, data = dats)
    betaest0 <- coef(res0)
    estreg <- cbind(1,meanx$mean01) %*% betaest0
    
    
    ########### expand matrix ##############
    
    xmp <- matrix(0, nrow(meanx), 2*N)
    xmp[,1] <- 1
    xmp[,2] <- meanx$mean01
    xmp[1,seq(3,2*N,by = 2)] <- -1
    xmp[1,seq(4,2*N,by = 2)] <- - meanx$mean01[1]
    
    xm <- matrix(0, nrow(dats), 2*(N))
    xm[,1] <- 1
    xm[,2] <- dats$x
    xm[1:nv[1],seq(3,ncol(xm),by = 2)] <- -1
    xm[1:nv[1],seq(4,ncol(xm),by = 2)] <- -dats$x[1:nv[1]]
    
    for(i in 2:N)
    {
      index1 <- sum(nv[1:(i-1)]) + 1
      index2 <- sum(nv[1:i])
      xm[index1:index2,(i-1)*2 +1] <- 1
      xm[index1:index2,2*i] <- dats$x[dats$dom==i]
      
      xmp[i, (i-1)*2 +1]  <- 1
      xmp[i,2*i]  <- meanx$mean01[i]
    }
    
    cv1 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights, alpha = 0, standardize = TRUE)
    betaest1 <- coef(cv1,s = cv1$lambda.min)
    estridge <- as.matrix(xmp %*% betaest1)[,1]
    
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    
    ww <- 1/abs(betaest1[-1])
    cv3 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1,
                     penalty.factor = ww)
    betaest3 <- coef(cv3,s = cv3$lambda.min)
    estalasso <- as.matrix(xmp%*% betaest3)[,1]
    
    rmse <- sqrt(c(mean((meany - estdir$Direct)^2),
                   mean((meany - estbhf$eblup$eblup)^2),
                   mean((meany - estreg)^2),
                   mean((meany - estridge)^2),
                   mean((meany - estlasso)^2),
                   mean((meany - estalasso)^2)))
    rmse_mat1[m,] <- rmse
  }
  
  out <- list(rmse = rmse_mat1, numcoef = numcoef1)
  return(out)
}



### first 50 counties 
# xname is covariate 
# yname is response 
# group: for counties, group 1 with lower sampling rate, group 2 with higher sampling rate
# popn(2 columns): county index and segment size
# meanx(2 columns): county index and total mean in x
# meany(2 columns): county index and total mean in y 

#
simsrs2 <- function(M, dat, popn, xname, yname, meanx, meany, 
                    srh = 0.03, srl = 0.03/5, seed = 451920)
{
  
  numcoef1 <- rep(0, M)
  rmse_mat1 <- matrix(0, M, 6)
  colnames(rmse_mat1) <- c("direct","eblup","reg","ridge","lasso","Alasso")
  Nv <- popn[,2]
  
  srv <- c(srl,srh)[group]
  
  N <- length(Nv) ## total number of counties
  nv <- round(Nv * srv) ## sample size 
  n0 <- sum(nv) ## total sample size
  
  
  dats <- as.data.frame(matrix(0, n0, 4)) # sampled data
  colnames(dats) <- c("dom", "y", "x", "weights")
  dats$dom <- rep(1:N, nv)
  
  
  for(m in 1:M)
  {
    set.seed(seed + m)
    for(i in 1:N)
    {
      indpi <- dat$county == i 
      ypi <- dat[indpi,yname] ## population y
      xpi <- dat[indpi,xname] ## population x
      
      indi <- dats$dom == i
      
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
      dats$weights[indi] <- Nv[i]/nv[i]
    }
    
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                       method = "REML",data = dats)
    
    res0 <- lm(y~x, weights = dats$weights, data = dats)
    
    betaest0 <- coef(res0)
    
    estreg <- cbind(1,meanx[,2]) %*% betaest0
    
    
    ########### expand matrix ##############
    
    xmp <- matrix(0, nrow(meanx), 2*N)
    xmp[,1] <- 1
    xmp[,2] <- meanx[,2]
    xmp[1,seq(3,2*N,by = 2)] <- -1
    xmp[1,seq(4,2*N,by = 2)] <- - meanx[1,2]
    
    xm <- matrix(0, nrow(dats), 2*(N))
    xm[,1] <- 1
    xm[,2] <- dats$x
    xm[1:nv[1],seq(3,ncol(xm),by = 2)] <- -1
    xm[1:nv[1],seq(4,ncol(xm),by = 2)] <- -dats$x[1:nv[1]]
    
    for(i in 2:N)
    {
      index1 <- sum(nv[1:(i-1)]) + 1
      index2 <- sum(nv[1:i])
      xm[index1:index2,(i-1)*2 +1] <- 1
      xm[index1:index2,2*i] <- dats$x[dats$dom==i]
      
      xmp[i, (i-1)*2 +1]  <- 1
      xmp[i,2*i]  <- meanx[i,2]
    }
    
    cv1 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights, alpha = 0, standardize = TRUE)
    betaest1 <- coef(cv1,s = cv1$lambda.min)
    estridge <- as.matrix(xmp %*% betaest1)[,1]
    
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    
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
    rmse_mat1[m,] <- rmse
  }
  
  out <- list(rmse = rmse_mat1, numcoef = numcoef1, samplesize = nv)
  return(out)
}