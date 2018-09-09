### simluate simple model to compare eblup and lasso 
setwd("/Users/Xin/Research/survey/")
nseg <- read.csv("data/nseg_urban.csv",stringsAsFactors = FALSE)
popn <- nseg[,c("COUNTY","nseg")]
#sr <- 0.2
#gamma <- 1
# sig2 <- 1
# srh <- 0.2
# srl <- 0.05
set.seed(1040)
group <- sample.int(2, size = nrow(popn),replace = TRUE)

###### constant mean model ####
# group 1 for h 2 for l
comp_epla <- function(mu, gamma, sig2, M, sr = 0.2, popn,
                      srh = NULL, srl = NULL, group = NULL)
{
 
  sig2v <- gamma * sig2 
  N <- nrow(popn)
  doms <- popn[,1]
  rmse_mat1 <- matrix(0, M, 5)
  colnames(rmse_mat1) <- c("direct","eblup","reg","lasso","lasso1")
  numcoef1 <- rep(0,M)
  
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
  
  datp <- as.data.frame(matrix(0, Nt, 2))
  indp <- rep(popn[,1], Nv)
  colnames(datp) <- c("dom","y")
  datp$dom <- indp
  
  meanx <- cbind(popn[,1], 1)
  
  dats <- as.data.frame(matrix(0, n0, 3))
  colnames(dats) <- c("dom", "y", "weights")
  dats$dom <- rep(popn[,1], nv)
  
  fn <- nv/Nv
  
  for(m in 1:M)
  {
    set.seed(m)
    
    y <- mu + rep(rnorm(N)*sqrt(sig2v), Nv) + rnorm(Nt)*sqrt(sig2)
    datp$y <- y
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      ypi <- y[indpi]
      
      indi <- dats$dom == doms[i]
      
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$weights[indi] <- Nv[i]/nv[i]
    }
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    estbhf <- eblupBHF(y~1, dom = dom, meanxpop = meanx ,popnsize = popn, 
                       method = "REML",data = dats)
    
    res0 <- lm(y~1, weights = dats$weights, data = dats)
    betaest0 <- coef(res0)
    #estreg <- betaest0
  
    estreg <- meansy$meansy*fn + (1- fn)*betaest0
    
    xmp <- matrix(0, N , N)
    xmp[,1] <- 1
    xmp[1,-1] <- -1
    
  
    xm <- matrix(0, nrow(dats), N)
    xm[,1] <- 1
    xm[1:nv[1],-1] <- -1
    
    for(i in 2:N)
    {
      index1 <- sum(nv[1:(i-1)]) + 1
      index2 <- sum(nv[1:i])
      xm[index1:index2,i] <- 1
      
      xmp[i, i]  <- 1
    }
    
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    estlasso1 <- meansy$meansy*fn + (1- fn)*estlasso
    
    rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg)^2),
                   mean((meany[,2] - estlasso)^2),
                   mean((meany[,2] - estlasso1)^2)))
    rmse_mat1[m,]<- rmse
  }
  
  return(list(rmse = rmse_mat1, numcoef = numcoef1))
}

res1 <- comp_epla(mu = 5,gamma = 1,sig2 = 1,M = 100,sr = 0.2,popn = popn)
res2 <- comp_epla(mu = 5,gamma = 0.5,sig2 = 1,M = 100,sr = 0.2,popn = popn)
res3 <- comp_epla(mu = 5,gamma = 4,sig2 = 1,M = 100,sr = 0.2,popn = popn)

gammav <- c(seq(0.2,1,by = 0.2), seq(1.5,5,by = 0.5))

ratiom <- ratiome <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla(mu = 5,gamma = gammav[g],sig2 = 1,M = 100,sr = 0.2,popn = popn)
  ratiom[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom, names = gammav, main = "srs: 0.2")
abline(h = 1)

boxplot(ratiome, names = gammav, main = "srs: 0.2")
abline(h = 1)


ratiom2 <- matrix(0, 100, length(gammav))
ratiome2 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla(mu = 5,gamma = gammav[g],sig2 = 1,M = 100,sr = 0.1,popn = popn)
  ratiom2[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome2[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom2, names = gammav)
abline(h = 1)

boxplot(ratiome2, names = gammav)
abline(h = 1)


#### high and low 
ratiom3 <- ratiome3 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla(mu = 5,gamma = gammav[g],sig2 = 1,M = 100,popn = popn,srh = 0.2,srl = 0.05,group = group)
  ratiom3[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome3[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom3,names = gammav)
abline(h = 1)
boxplot(ratiome3,names = gammav)
abline(h = 1)

ratiom4 <- ratiome4 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla(mu = 5,gamma = gammav[g],sig2 = 1,M = 100,popn = popn,srh = 0.2,srl = 0.025,group = group)
  ratiom4[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome4[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom4,names = gammav)
abline(h = 1)

boxplot(ratiome4,names = gammav)
abline(h = 1)

save(ratiom, ratiom2, ratiom3, ratiom4, file = "result/meanmodel_comp_epla.RData")

save(ratiome, ratiome2, ratiome3, ratiome4, file = "result/meanmodel_comp_epla_e.RData")

par(mfrow = c(2,2))
boxplot(ratiom, names = gammav, main = "srs: 0.2")
abline(h = 1)

boxplot(ratiom2,names = gammav, main = "srs: 0.1")
abline(h = 1)

boxplot(ratiom3,names = gammav, main = "srs: 0.2 and 0.05")
abline(h = 1)

boxplot(ratiom4,names = gammav, main = "srs: 0.2 and 0.025")
abline(h = 1)

par(mfrow = c(2,2))
boxplot(ratiome, names = gammav, main = "srs: 0.2")
abline(h = 1)

boxplot(ratiome2,names = gammav, main = "srs: 0.1")
abline(h = 1)

boxplot(ratiome3,names = gammav, main = "srs: 0.2 and 0.05")
abline(h = 1)

boxplot(ratiome4,names = gammav, main = "srs: 0.2 and 0.025")
abline(h = 1)




########### two groups ############

comp_epla_v2 <- function(muv, gamma, sig2, M, sr = 0.2, popn,
                      srh = NULL, srl = NULL, group = NULL)
{
  
  ng <- length(muv)
  sig2v <- gamma * sig2 
  nc <- nrow(popn)
  doms <- popn[,1]
  rmse_mat1 <- matrix(0, M, 5)
  colnames(rmse_mat1) <- c("direct","eblup","reg","lasso","lasso1")
  numcoef1 <- rep(0,M)
  
  Nv <- popn[,2]
  Nt <- sum(Nv)
  
  N <- length(Nv) ## total number of counties
  
  if(!is.null(group)) ## two groups
  {
    srv <- c(srh, srl)[group]
    nv <- round(Nv*srv)
  }else{
    nv <- round(Nv * sr) ## sample size 
  }
  
  fn <- nv/Nv
  
  n0 <- sum(nv) ## total sample size
  
  datp <- as.data.frame(matrix(0, Nt, 2))
  indp <- rep(popn[,1], Nv)
  colnames(datp) <- c("dom","y")
  datp$dom <- indp
  
  meanx <- cbind(popn[,1], 1)
  
  dats <- as.data.frame(matrix(0, n0, 3))
  colnames(dats) <- c("dom", "y", "weights")
  dats$dom <- rep(popn[,1], nv)
  
  for(m in 1:M)
  {
    set.seed(m)
    mugroup <- sample.int(ng, size = nc,replace = TRUE)
    
    y <- rep(muv[mugroup], Nv) + rep(rnorm(nc)*sqrt(sig2v), Nv) + rnorm(Nt)*sqrt(sig2)
    datp$y <- y
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      ypi <- y[indpi]
      
      indi <- dats$dom == doms[i]
      
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$weights[indi] <- Nv[i]/nv[i]
    }
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    estbhf <- eblupBHF(y~1, dom = dom, meanxpop = meanx ,popnsize = popn, 
                       method = "REML",data = dats)
    
    res0 <- lm(y~1, weights = dats$weights, data = dats)
    betaest0 <- coef(res0)
    #estreg <- betaest0
    
    
    estreg <- meansy$meansy*fn + (1- fn)*betaest0
    
    
    
    xmp <- matrix(0, nc , nc)
    xmp[,1] <- 1
    xmp[1,-1] <- -1
    
    
    xm <- matrix(0, nrow(dats), nc)
    xm[,1] <- 1
    xm[1:nv[1],-1] <- -1
    
    for(i in 2:nc)
    {
      index1 <- sum(nv[1:(i-1)]) + 1
      index2 <- sum(nv[1:i])
      xm[index1:index2,i] <- 1
      
      xmp[i, i]  <- 1
    }
    
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    estlasso1 <- meansy$meansy*fn + (1- fn)*estlasso
    
    rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg)^2),
                   mean((meany[,2] - estlasso)^2),
                   mean((meany[,2] - estlasso1)^2)))
    rmse_mat1[m,]<- rmse
  }
  
  return(list(rmse = rmse_mat1, numcoef = numcoef1))
}


resv2 <- comp_epla_v2(muv = c(5,4),gamma = 1,sig2 = 1,M = 100,sr = 0.2,popn = popn)
boxplot(resv2$rmse[,4]/resv2$rmse[,2])
boxplot(resv2$rmse[,5]/resv2$rmse[,2])

ratiom21 <- matrix(0, 100, length(gammav))
ratiome21 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla_v2(muv = c(5,4),gamma = gammav[g],sig2 = 1,M = 100,sr = 0.2,popn = popn)
  ratiom21[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome21[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom21)
abline(h = 1)

boxplot(ratiome21)
abline(h = 1)

ratiom22 <- matrix(0, 100, length(gammav))
ratiome22 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla_v2(muv = c(5,4),gamma = gammav[g],sig2 = 1,M = 100,sr = 0.1,popn = popn)
  ratiom22[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome22[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom22)
abline(h = 1)
boxplot(ratiome22)
abline(h = 1)

ratiom23 <- matrix(0, 100, length(gammav))
ratiome23 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla_v2(muv = c(5,4),gamma = gammav[g],sig2 = 1,M = 100,popn = popn,srh = 0.2,srl = 0.05,group = group)
  ratiom23[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome23[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom23)
abline(h = 1)

boxplot(ratiome23)
abline(h = 1)

save(ratiom21, ratiom22, ratiom23,ratiome21, ratiome22, ratiome23, file = "result/meanmodel_comp_epla2.RData")


ratiom31 <- ratiome31 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla_v2(muv = c(5,1),gamma = gammav[g],sig2 = 1,M = 100,sr = 0.2,popn = popn)
  ratiom31[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome31[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom31)
abline(h = 1)

boxplot(ratiome31)
abline(h = 1)

ratiom32 <- ratiome32 <- matrix(0, 100, length(gammav))
for(g in 1:length(gammav))
{
  resg <- comp_epla_v2(muv = c(5,1),gamma = gammav[g],sig2 = 1,M = 100,sr = 0.1,popn = popn)
  ratiom32[,g] <- resg$rmse[,4]/resg$rmse[,2]
  ratiome32[,g] <- resg$rmse[,5]/resg$rmse[,2]
}

boxplot(ratiom32)
abline(h = 1)

boxplot(ratiome32)
abline(h = 1)

save(ratiom31, ratiom32, ratiome31, ratiome32, file = "result/meanmodel_comp_epla3.RData")
