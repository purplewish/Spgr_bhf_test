###### regression simulation ##### 
#### subgroup analysis can also be used 
## simulated them from regression 

library(glmnet)
## group for regression coefficients is randomly generated 
betam <- matrix(c(1,0.5,1,1),2)

# group is from sampling group ####

setwd("/Users/Xin/Research/survey/")
nseg <- read.csv("data/nseg_urban.csv",stringsAsFactors = FALSE)
popn <- nseg[,c("COUNTY","nseg")]

sigx <- 1
sige <- 1
sr <- 0.2
M <- 10

comp_epla_reg <- function(betam, mux = 1, sigx, sige, M, sr = 0.2, popn,
                      srh = NULL, srl = NULL, group = NULL, seed = 1630)
{
  
  N <- nrow(popn)  #### number of domains
  doms <- popn[,1]
  rmse_mat1 <- matrix(0, M, 8)
  colnames(rmse_mat1) <- c("direct","eblup","reg","reg1","reg2","lasso","lasso1","lasso2")
  numcoef1 <- rep(0,M)
  
  ng <- ncol(betam)
  
  groupdom <- sample.int(n = ng,size = N,replace = TRUE)
  betamdom <- betam[,groupdom]
  
  Nv <- popn[,2]
  Nt <- sum(Nv)

  
  if(!is.null(group)) ## two groups
  {
    srv <- c(srh, srl)[group]
    nv <- round(Nv*srv)
  }else{
    nv <- round(Nv * sr) ## sample size 
  }
  
  n0 <- sum(nv) ## total sample size
  
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(popn[,1], Nv)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  
  meanx <- cbind(popn[,1], 1)
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x","weights")
  dats$dom <- rep(popn[,1], nv)
  
  fn <- nv/Nv
  
  for(m in 1:M)
  {
    set.seed(m)
    
    x <- rnorm(Nt)*sigx + mux
    datp$x <- x
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      xpi <- x[indpi]
      ypi <- cbind(1,xpi)%*%betamdom[,i] + rnorm(Nv[i])*sige
      datp$y[indpi] <- ypi
    
      indi <- dats$dom == doms[i]
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
      dats$weights[indi] <- Nv[i]/nv[i]
    }
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    meansx <- ddply(dats,.(dom),summarize, meansx = mean(x)) ## sample mean x
    
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                       method = "REML",data = dats)
    
    res0 <- lm(y~x, weights = dats$weights, data = dats)
    betaest0 <- coef(res0)

    estreg <- cbind(1,meanx$meanx)%*% betaest0
    
    estreg1 <- meansy$meansy*fn + cbind(1,(meanx$meanx - fn*meansx$meansx))%*%betaest0
    
    estreg2 <- meansy$meansy + cbind(1,(meanx$meanx - meansx$meansx))%*%betaest0
    
    xmp <- matrix(0, N, 2*N)
    xmp[,1] <- 1
    xmp[,2] <- meanx$meanx
    xmp[1,seq(3,2*N,by = 2)] <- -1
    xmp[1,seq(4,2*N,by = 2)] <- - meanx$meanx[1]
    
    xms <- matrix(0, N, 2*N)
    xms[,1] <- 1
    xms[,2] <- meansx$meansx
    xms[1,seq(3,2*N,by = 2)] <- -1
    xms[1,seq(4,2*N,by = 2)] <- -meansx$meansx[1]
    
    
    xm <- matrix(0, nrow(dats), 2*(N))
    xm[,1] <- 1
    xm[,2] <- dats$x
    xm[1:nv[1],seq(3,ncol(xm),by = 2)] <- -1
    xm[1:nv[1],seq(4,ncol(xm),by = 2)] <- -dats$x[dats$dom == doms[1]]
    
    for(i in 2:N)
    {
      index1 <- sum(nv[1:(i-1)]) + 1
      index2 <- sum(nv[1:i])
      xm[index1:index2,(i-1)*2 +1] <- 1
      xm[index1:index2,2*i] <- dats$x[dats$dom==doms[i]]
      
      xmp[i, (i-1)*2 +1]  <- 1
      xmp[i,2*i]  <- meanx$meanx[i]
      
      xms[i, (i-1)*2 +1]  <- 1
      xms[i,2*i]  <- meansx$meansx[i]
    }
    
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    estlasso1 <- meansy$meansy*fn + estlasso - fn*as.matrix(xms %*% betaest2)[,1]
    
    estlasso2 <- meansy$meansy + estlasso - as.matrix(xms %*% betaest2)[,1]
    
    rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg)^2),
                   mean((meany[,2] - estreg1)^2),
                   mean((meany[,2] - estreg2)^2),
                   mean((meany[,2] - estlasso)^2),
                   mean((meany[,2] - estlasso1)^2),
                   mean((meany[,2] - estlasso2)^2))
                 )
    rmse_mat1[m,]<- rmse
  }
  
  return(list(rmse = rmse_mat1, numcoef = numcoef1))
}

res1 <- comp_epla_reg(betam = betam, mux = 1,sigx = 1,sige = 1,M = 10,sr = 0.2,popn = popn)


#### a sequence of values of sigx 
sigxv <- seq(0.5,3,by = 0.5)
M <- 100

rmse_array <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam, mux = 1,sigx = sigxv[j],sige = 1,M = 100,sr = 0.2,popn = popn)
  rmse_array[,,j] <- resj$rmse
  numcoef_mat[,j] <- resj$numcoef
}


boxplot(rmse_array[,6,1:5]/rmse_array[,2,1:5])
abline(h=1)
boxplot(rmse_array[,7,1:5]/rmse_array[,2,1:5])
abline(h=1)
boxplot(rmse_array[,8,1:5]/rmse_array[,2,1:5])


rmse_array2 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat2 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam, mux = 1,sigx = sigxv[j],sige = 1,M = 100,sr = 0.1,popn = popn)
  rmse_array2[,,j] <- resj$rmse
  numcoef_mat2[,j] <- resj$numcoef
}

boxplot(rmse_array2[,6,1:5]/rmse_array2[,2,1:5])
boxplot(rmse_array2[,7,1:5]/rmse_array2[,2,1:5])
boxplot(rmse_array2[,8,1:5]/rmse_array2[,2,1:5])


set.seed(1040)
groups <- sample.int(2, size = nrow(popn),replace = TRUE)

rmse_array3 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat3 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam, mux = 1,sigx = sigxv[j],sige = 1,M = 100,popn = popn, srh = 0.2,srl = 0.1,group = groups)
  rmse_array3[,,j] <- resj$rmse
  numcoef_mat3[,j] <- resj$numcoef
}

boxplot(rmse_array3[,6,1:5]/rmse_array3[,2,1:5])
boxplot(rmse_array3[,7,1:5]/rmse_array3[,2,1:5])
boxplot(rmse_array3[,8,1:5]/rmse_array3[,2,1:5])

rmse_array4 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat4 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam, mux = 1,sigx = sigxv[j],sige = 1,M = 100,popn = popn, srh = 0.2,srl = 0.05,group = groups)
  rmse_array4[,,j] <- resj$rmse
  numcoef_mat4[,j] <- resj$numcoef
}

boxplot(rmse_array4[,6,1:5]/rmse_array4[,2,1:5])
boxplot(rmse_array4[,7,1:5]/rmse_array4[,2,1:5])
boxplot(rmse_array4[,8,1:5]/rmse_array4[,2,1:5])


save(rmse_array, rmse_array2, rmse_array3, rmse_array4, file = "result/rmse_reg0.RData")
#### another betam

betam1 <- matrix(c(1,0.2,1,1),2)

rmse_array11 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat11 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam1, mux = 1,sigx = sigxv[j],sige = 1,M = 100,sr = 0.2,popn = popn)
  rmse_array11[,,j] <- resj$rmse
  numcoef_mat11[,j] <- resj$numcoef
}

boxplot(rmse_array11[,6,1:5]/rmse_array11[,2,1:5])
boxplot(rmse_array11[,7,1:5]/rmse_array11[,2,1:5])
boxplot(rmse_array11[,8,1:5]/rmse_array11[,2,1:5])


rmse_array12 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat12 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam1, mux = 1,sigx = sigxv[j],sige = 1,M = 100,sr = 0.1,popn = popn)
  rmse_array12[,,j] <- resj$rmse
  numcoef_mat12[,j] <- resj$numcoef
}

boxplot(rmse_array12[,6,1:5]/rmse_array12[,2,1:5])
boxplot(rmse_array12[,7,1:5]/rmse_array12[,2,1:5])
boxplot(rmse_array12[,8,1:5]/rmse_array12[,2,1:5])


rmse_array13 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat13 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam1, mux = 1,sigx = sigxv[j],sige = 1,M = 100,popn = popn,srh = 0.2,srl = 0.1,group = groups)
  rmse_array13[,,j] <- resj$rmse
  numcoef_mat13[,j] <- resj$numcoef
}

boxplot(rmse_array13[,6,1:5]/rmse_array13[,2,1:5])
boxplot(rmse_array13[,7,1:5]/rmse_array13[,2,1:5])
boxplot(rmse_array13[,8,1:5]/rmse_array13[,2,1:5])

rmse_array14 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat14 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam1, mux = 1,sigx = sigxv[j],sige = 1,M = 100,popn = popn,srh = 0.2,srl = 0.05,group = groups)
  rmse_array14[,,j] <- resj$rmse
  numcoef_mat14[,j] <- resj$numcoef
}

boxplot(rmse_array14[,6,1:5]/rmse_array14[,2,1:5])
abline(h = 1)
boxplot(rmse_array14[,7,1:5]/rmse_array14[,2,1:5])
abline(h = 1)
boxplot(rmse_array14[,8,1:5]/rmse_array14[,2,1:5])

save(rmse_array11,rmse_array12, rmse_array12, rmse_array14, file = "result/rmse_reg1.RData")

##### three groups #####

betam2 <- matrix(c(1,0.5,1,1,1,1.5),nrow = 2)

rmse_array21 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat21 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam2, mux = 1,sigx = sigxv[j],sige = 1,M = 100,sr = 0.2,popn = popn)
  rmse_array21[,,j] <- resj$rmse
  numcoef_mat21[,j] <- resj$numcoef
}


boxplot(rmse_array21[,6,1:5]/rmse_array21[,2,1:5])
abline(h = 1)
boxplot(rmse_array21[,7,1:5]/rmse_array21[,2,1:5])
abline(h = 1)
boxplot(rmse_array21[,8,1:5]/rmse_array21[,2,1:5])
abline(h = 1)


rmse_array22 <- array(0, dim = c(M, 8, length(sigxv)))
numcoef_mat22 <-matrix(0, M, length(sigxv))
for(j in 1:length(sigxv))
{
  resj <- comp_epla_reg(betam = betam2, mux = 1,sigx = sigxv[j],sige = 1,M = 100,sr = 0.1,popn = popn)
  rmse_array22[,,j] <- resj$rmse
  numcoef_mat22[,j] <- resj$numcoef
}


boxplot(rmse_array22[,6,1:5]/rmse_array22[,2,1:5])
abline(h = 1)
boxplot(rmse_array22[,7,1:5]/rmse_array22[,2,1:5])
abline(h = 1)
boxplot(rmse_array22[,8,1:5]/rmse_array22[,2,1:5])
abline(h = 1)

save(rmse_array21,rmse_array22, file = "result/rmse_reg2.RData")
