library(sae)
library(plyr)
library(Spgr)
library(Matrix)
setwd("/Users/Xin/Research/NRI_urban/")

source('/Users/Xin/Research/cluster/code/refit.R')
source("/Users/Xin/Research/NRI_urban/code/Spgr_bhf.R")

nsegtotal <- read.csv("data/nsegtotal.csv",stringsAsFactors = FALSE)

load("data/adjMat.RData")
Cmatia <- adjMat[["19"]]

countyfips <- as.character(19000+ nsegtotal$COUNTYFP)
Cmatia <- Cmatia[countyfips, countyfips]
ordermat <- getorder(Matrix(Cmatia))



M <- 100
popn <- nsegtotal[,c("COUNTYFP","nseg")]
mux <- 1
sdx <- 1
sde <- 1 

#### model-based evaluation ####
# rate or ni specify 1 
# M is the total number of simulations
sim_pop <- function(beta, mux, sdx, sde, rate = 0.1, ni = NULL, M, popn, group,
                    alp = c(0.5,1,1.5), seed = 7231625)
{
  nobs <- nrow(nsegtotal)
  msemat <- matrix(0, nobs, 8)
  colnames(msemat) <- c("direct","bhf","reg1","reg2","equal1","equal2","sp1","sp2")
  
  
  doms <- popn[,1]
  Nv <- popn[,2]
  Nt <- sum(Nv)
  N <- nrow(popn)
  
  ng <- nrow(beta)
  betamdom <- beta[group,]
  
  
  #### population matrix
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(doms, Nv)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  
  ## sample 
  
  if(is.null(rate) & is.null(ni))
  {
    print("specify rate or sample size")
  }
  
  if(!is.null(rate) & is.null(ni))
  {
    nv <- round(Nv * rate)
  }
  
  if(is.null(rate) & !is.null(ni))
  {
    nv <- rep(ni, N)
  }
  
  
  n0 <- sum(nv) ## total sample size
  
  fn <- nv/Nv
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x","weights")
  dats$dom <- rep(doms, nv)
  dats$weights <-  rep(Nv/nv, nv)
  
  nrnzero <- 0 
  groupmat1 <- groupmat2 <- matrix(0, nobs, M)

  for(m in 1:M)
  {
    set.seed(seed + m)
    x <- rnorm(Nt)*sdx + mux
    datp$x <- x
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      xpi <- x[indpi]
      ypi <- cbind(1,xpi)%*%betamdom[i,] + rnorm(Nv[i])*sde
      datp$y[indpi] <- ypi
      
      indi <- dats$dom == doms[i]
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
    }
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    meansx <- ddply(dats,.(dom),summarize, meansx = mean(x)) ## sample mean x
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    

    estbhf <- 0
    tryCatch({estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                       method = "REML",data = dats)},
             error=function(e){cat("ERROR :",error(e), "\n")})
    
    ## regression 
    betaest_reg <- coef(lm(y~x, dats))
    
    reg01 <- cbind(1, meanx$meanx)%*%betaest_reg
    reg02 <- meansy$meansy*fn + reg01 - fn*cbind(1, meansx$meansx)%*%betaest_reg
             
    
    # equal weight ###
    c1 <- rep(1,nobs*(nobs-1)/2)
    lam1 <- seq(0.05,1,by = 0.05)
    bic1 <- rep(0,length(lam1))
    beta_array1 <- array(0, dim = c(nobs,2,length(lam1)))
    group1 <- matrix(0,nobs,length(lam1))
    
    betam01 <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
    for(l in 1:length(lam1))
    {
      res1 <- Spgrrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam01,lam = lam1[l],maxiter = 2000)
      betam01 <- res1$beta
      beta_array1[,,l] <- res1$beta
      bic1[l] <- BICcrx(obj = res1,indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),c0 = 0.2)
      group1[,l] <- res1$group
    }
    betaest1 <- beta_array1[,,which.min(bic1)]
    
    groupmat1[,m] <- group1[,which.min(bic1)]
    
  
    reg11 <- rowSums(cbind(1, meanx$meanx)*betaest1)
    reg12 <- meansy$meansy*fn + reg11 - fn*rowSums(cbind(1, meansx$meansx)*betaest1)
    ### spatial weight ####
    
    lam2 <- seq(0.05,1.5,by =0.05)
    bic2 <- rep(0, length(alp))
    beta_array2 <- array(0, dim = c(nobs,2,length(alp)))
    group2 <- matrix(0, nobs, length(alp))
    
    for(l1 in 1:length(alp))
    {
      bic2l <- rep(0,length(lam2))
      beta_array2l <- array(0, dim = c(nobs,2,length(lam2)))
      c2l <- exp(alp[l1]*(1 - ordermat))
      group2l <- matrix(0,nobs, length(lam2))
      
      betam02l <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
      for(l2 in 1:length(lam2))
      {
        res2l <- Spgrrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c2l,betam0 = betam02l,lam = lam2[l2],maxiter = 2000)
        betam02l <- res2l$beta
        beta_array2l[,,l2] <- res2l$beta
        bic2l[l2] <- BICcrx(obj = res2l,indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),c0 = 0.2)
        BICcrx(obj = res1,indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),c0 = 0.2)
        group2l[,l2] <- res2l$group
      }
      
    bic2[l1] <- min(bic2l)
    beta_array2[,,l1] <- beta_array2l[,,which.min(bic2l)]
    group2[,l1] <- group2l[,which.min(bic2l)]
      
    }

    betaest2 <- beta_array2[,,which.min(bic2)]
    groupmat2[,m] <- group2[,which.min(bic2)]
    
    reg21 <- rowSums(cbind(1, meanx$meanx)*betaest2)
    reg22 <-  meansy$meansy*fn + reg21 - fn*rowSums(cbind(1, meansx$meansx)*betaest2) 
    
    if(!is.na(estbhf$eblup)[1])
    {
      nrnzero <- nrnzero + 1
      msemat <-  msemat + cbind((estdir$Direct - meany[,2])^2,
                                (estbhf$eblup$eblup - meany[,2])^2,
                                (reg01 - meany[,2])^2,
                                (reg02 - meany[,2])^2,
                                (reg11 - meany[,2])^2,
                                (reg12 - meany[,2])^2,
                                (reg21 - meany[,2])^2,
                                (reg22 - meany[,2])^2)
      boxplot(msemat[,3:6])
      print(m)
    }
    

    
  }
  
  msemat <- msemat/nrnzero
  return(list(mse = msemat, nrnzero = nrnzero, group1 = groupmat1, group2 = groupmat2))
}


beta <- matrix(c(0.5,0.5,1.5,1.5,2.5,2.5),ncol=2,byrow = TRUE)
sim_d1_r05_sde1 <- sim_pop(beta = beta, mux = 1,sdx = 1,sde = 1,rate = 0.05,M = 100,popn = popn,group = nsegtotal$newgroup)

sim_d1_r05_sde2 <- sim_pop(beta = beta, mux = 1,sdx = 1,sde = 2,rate = 0.05,M = 100,popn = popn,group = nsegtotal$newgroup)

sim_d1_r01_sde1 <- sim_pop(beta = beta, mux = 1,sdx = 1,sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)

sim_d1_r01_sde2 <- sim_pop(beta = beta, mux = 1,sdx = 1,sde = 2,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)


pdf("doc/figures/fig_sae_unit_d1.pdf",width = 9,height = 9)
par(mfrow = c(2,2))
boxplot(sqrt(sim_d1_r05_sde1$mse[,-c(1,3,4)]),ylab = "RMSE", 
        main = expression(paste("r=0.05, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_d1_r05_sde2$mse[,-c(1,3,4)]),ylab = "RMSE",
        main = expression(paste("r=0.05, ",sigma[epsilon],"=2")))
boxplot(sqrt(sim_d1_r01_sde1$mse[,-c(1,3,4)]),ylab = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_d1_r01_sde2$mse[,-c(1,3,4)]), ylab  = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=2")))
dev.off()

# small 
beta2 <- matrix(c(0.5,0.5,1,1,1.5,1.5),ncol=2,byrow = TRUE)
sim_d05_r05_sde1 <- sim_pop(beta = beta2, mux = 1,sdx = 1,sde = 1,rate = 0.05,M = 100,popn = popn,group = nsegtotal$newgroup)
boxplot(sqrt(sim_d05_r05_sde1$mse[,-1]))

sim_d05_r05_sde2 <- sim_pop(beta = beta2, mux = 1,sdx = 1,sde = 2,rate = 0.05,M = 100,popn = popn,group = nsegtotal$newgroup)
boxplot(sqrt(sim_d05_r05_sde2$mse[,-1]))

# smaller sampling rate 
sim_d05_r01_sde1 <- sim_pop(beta = beta2, mux = 1,sdx = 1,sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup,alp = c(0.5,1,1.5))
boxplot(sqrt(sim_d05_r01_sde1$mse[,-1]))

sim_d05_r01_sde2 <- sim_pop(beta = beta2, mux = 1,sdx = 1,sde = 2,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup,alp = c(0.5,1,1.5))
boxplot(sqrt(sim_d05_r01_sde2$mse[,-1]))

pdf("doc/figures/fig_sae_unit_d05.pdf",width = 9,height = 9)
par(mfrow = c(2,2))
boxplot(sqrt(sim_d05_r05_sde1$mse[,-c(1,3,4)]),ylab = "RMSE", 
        main = expression(paste("r=0.05, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_d05_r05_sde2$mse[,-c(1,3,4)]),ylab = "RMSE",
        main = expression(paste("r=0.05, ",sigma[epsilon],"=2")))
boxplot(sqrt(sim_d05_r01_sde1$mse[,-c(1,3,4)]),ylab = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_d05_r01_sde2$mse[,-c(1,3,4)]), ylab  = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=2")))
dev.off()


boxplot(sqrt(sim_d05_r05_sde1$mse[,-c(1)]),ylab = "RMSE", 
        main = expression(paste("r=0.05, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_d05_r05_sde2$mse[,-c(1)]),ylab = "RMSE",
        main = expression(paste("r=0.05, ",sigma[epsilon],"=2")))
boxplot(sqrt(sim_d05_r01_sde1$mse[,-c(1)]),ylab = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_d05_r01_sde2$mse[,-c(1)]), ylab  = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=2")))

save(sim_d1_r05_sde1,sim_d1_r05_sde2,sim_d1_r01_sde1,sim_d1_r01_sde2,
     sim_d05_r05_sde1,sim_d05_r05_sde2,sim_d05_r01_sde1,sim_d05_r01_sde2, 
     file = "result/sim_subgroup3_sae_unit.RData")



load("result/sim_subgroup3_sae_unit.RData")





### one group 

beta1 <- matrix(c(1,1),ncol=2)
sim_g1_r01_sde1 <- sim_pop(beta = beta1, mux = 1,sdx = 1,sde = 1,rate = 0.01,M = 100,popn = popn,group = rep(1,nobs),alp = c(0.5,1,1.5))


sim_g1_r01_sde2 <- sim_pop(beta = beta1, mux = 1,sdx = 1,sde = 2,rate = 0.01,M = 100,popn = popn,group = rep(1,nobs),alp = c(0.5,1,1.5))


save(sim_g1_r01_sde1, sim_g1_r01_sde2, file = "result/sim_subgroup1_sae_unit.RData")
load("result/sim_subgroup1_sae_unit.RData")


pdf("doc/figures/fig_sae_unit_group1.pdf",width = 6,height = 6)
par(mfrow = c(2,1))
boxplot(sqrt(sim_g1_r01_sde1$mse[,-1]),ylab = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=1")))
boxplot(sqrt(sim_g1_r01_sde2$mse[,-1]), ylab = "RMSE",
        main = expression(paste("r=0.01, ",sigma[epsilon],"=2")))
dev.off()


#### simulate data from bhf model ####
sim_pop_bhf <- function(beta, mux, sdx, sdv, sde, rate = 0.1, ni = NULL, M, popn, group,
                    alp = c(0.5,1,1.5), seed = 7231625)
{
  nobs <- nrow(nsegtotal)
  msemat <- matrix(0, nobs, 8)
  colnames(msemat) <- c("direct","bhf","reg1","reg2","equal1","equal2","sp1","sp2")
  
  doms <- popn[,1]
  Nv <- popn[,2]
  Nt <- sum(Nv)
  N <- nrow(popn)
  
  ng <- nrow(beta)
  betamdom <- beta[group,]
  
  #### population matrix
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(doms, Nv)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  
  ## sample 
  
  if(is.null(rate) & is.null(ni))
  {
    print("specify rate or sample size")
  }
  
  if(!is.null(rate) & is.null(ni))
  {
    nv <- round(Nv * rate)
  }
  
  if(is.null(rate) & !is.null(ni))
  {
    nv <- rep(ni, N)
  }
  
  n0 <- sum(nv) ## total sample size
  fn <- nv/Nv
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x","weights")
  dats$dom <- rep(doms, nv)
  dats$weights <-  rep(Nv/nv, nv)
  
  nrnzero <- 0 
  groupmat1 <- groupmat2 <- matrix(0, nobs, M)
  
  for(m in 1:M)
  {
    set.seed(seed + m)
    x <- rnorm(Nt)*sdx + mux
    vrand <- rnorm(N)*sdv
    datp$x <- x
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      xpi <- x[indpi]
      ypi <- cbind(1,xpi)%*%betamdom[i,] + vrand[i] + rnorm(Nv[i])*sde
      datp$y[indpi] <- ypi
      
      indi <- dats$dom == doms[i]
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
    }
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    meansx <- ddply(dats,.(dom),summarize, meansx = mean(x)) ## sample mean x
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    
    estbhf <- 0
    tryCatch({estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                                 method = "REML",data = dats)},
             error=function(e){cat("ERROR :",error(e), "\n")})
    
    ## regression 
    betaest_reg <- coef(lm(y~x, dats))
    
    reg01 <- cbind(1, meanx$meanx)%*%betaest_reg
    reg02 <- meansy$meansy*fn + reg01 - fn*cbind(1, meansx$meansx)%*%betaest_reg
    
    
    # equal weight ###
    c1 <- rep(1,nobs*(nobs-1)/2)
    lam1 <- seq(0.05,1,by = 0.05)
    bic1 <- rep(0,length(lam1))
    beta_array1 <- array(0, dim = c(nobs,2,length(lam1)))
    group1 <- matrix(0,nobs,length(lam1))
    
    betam01 <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
    for(l in 1:length(lam1))
    {
      res1 <- Spgrrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam01,lam = lam1[l],maxiter = 2000)
      betam01 <- res1$beta
      beta_array1[,,l] <- res1$beta
      bic1[l] <- BICcrx(obj = res1,indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),c0 = 0.2)
      group1[,l] <- res1$group
    }
    betaest1 <- beta_array1[,,which.min(bic1)]
    
    groupmat1[,m] <- group1[,which.min(bic1)]
    
    
    reg11 <- rowSums(cbind(1, meanx$meanx)*betaest1)
    reg12 <- meansy$meansy*fn + reg11 - fn*rowSums(cbind(1, meansx$meansx)*betaest1)
    ### spatial weight ####
    
    lam2 <- seq(0.05,1.5,by =0.05)
    bic2 <- rep(0, length(alp))
    beta_array2 <- array(0, dim = c(nobs,2,length(alp)))
    group2 <- matrix(0, nobs, length(alp))
    
    for(l1 in 1:length(alp))
    {
      bic2l <- rep(0,length(lam2))
      beta_array2l <- array(0, dim = c(nobs,2,length(lam2)))
      c2l <- exp(alp[l1]*(1 - ordermat))
      group2l <- matrix(0,nobs, length(lam2))
      
      betam02l <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
      for(l2 in 1:length(lam2))
      {
        res2l <- Spgrrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c2l,betam0 = betam02l,lam = lam2[l2],maxiter = 2000)
        betam02l <- res2l$beta
        beta_array2l[,,l2] <- res2l$beta
        bic2l[l2] <- BICcrx(obj = res2l,indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),c0 = 0.2)
        BICcrx(obj = res1,indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),c0 = 0.2)
        group2l[,l2] <- res2l$group
      }
      
      bic2[l1] <- min(bic2l)
      beta_array2[,,l1] <- beta_array2l[,,which.min(bic2l)]
      group2[,l1] <- group2l[,which.min(bic2l)]
      
    }
    
    betaest2 <- beta_array2[,,which.min(bic2)]
    groupmat2[,m] <- group2[,which.min(bic2)]
    
    reg21 <- rowSums(cbind(1, meanx$meanx)*betaest2)
    reg22 <-  meansy$meansy*fn + reg21 - fn*rowSums(cbind(1, meansx$meansx)*betaest2) 
    
    if(!is.na(estbhf$eblup)[1])
    {
      nrnzero <- nrnzero + 1
      msemat <-  msemat + cbind((estdir$Direct - meany[,2])^2,
                                (estbhf$eblup$eblup - meany[,2])^2,
                                (reg01 - meany[,2])^2,
                                (reg02 - meany[,2])^2,
                                (reg11 - meany[,2])^2,
                                (reg12 - meany[,2])^2,
                                (reg21 - meany[,2])^2,
                                (reg22 - meany[,2])^2)
      boxplot(msemat[,3:6])
      print(m)
    }
    
    
    
  }
  
  msemat <- msemat/nrnzero
  return(list(mse = msemat, nrnzero = nrnzero, group1 = groupmat1, group2 = groupmat2))
}

beta2 <- matrix(c(0.5,0.5,1,1,1.5,1.5),ncol=2,byrow = TRUE)
sim_bhf_d05_r05_sde1 <- sim_pop_bhf(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.05,M = 100,popn = popn,group = nsegtotal$newgroup)
boxplot(sim_bhf_d05_r05_sde1$mse[,-c(1,3,4)])

sim_bhf_d05_r05_sde2 <- sim_pop_bhf(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 2,rate = 0.05,M = 100,popn = popn,group = nsegtotal$newgroup)
boxplot(sim_bhf_d05_r05_sde2$mse[,-c(1,3,4)])


#### simulate data from bhf model with subgroups estimation ####
library(Rcpp)
sourceCpp("code/Spgr_bhf.cpp")
sourceCpp("code/Spgr_bhf3.cpp")
sim_pop_bhf2 <- function(beta, mux, sdx, sdv, sde, rate = 0.1, ni = NULL, M, popn, group,
                        alp = c(0.5,1,1.5), c0=0.2,seed = 7231625)
{
  nobs <- nrow(nsegtotal)
  msemat <- matrix(0, nobs, 6)
  colnames(msemat) <- c("direct","bhf","equal1","equal2","sp1","sp2")
  
  doms <- popn[,1]
  Nv <- popn[,2]
  Nt <- sum(Nv)
  N <- nrow(popn)
  
  ng <- nrow(beta)
  betamdom <- beta[group,]
  
  #### population matrix
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(doms, Nv)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  
  ## sample 
  
  if(is.null(rate) & is.null(ni))
  {
    print("specify rate or sample size")
  }
  
  if(!is.null(rate) & is.null(ni))
  {
    nv <- round(Nv * rate)
  }
  
  if(is.null(rate) & !is.null(ni))
  {
    nv <- rep(ni, N)
  }
  
  n0 <- sum(nv) ## total sample size
  fn <- nv/Nv
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x","weights")
  dats$dom <- rep(doms, nv)
  dats$weights <-  rep(Nv/nv, nv)
  
  nrnzero <- 0 
  groupmat1 <- groupmat2 <- matrix(0, nobs, M)
  
  for(m in 1:M)
  {
    set.seed(seed + m)
    x <- rnorm(Nt)*sdx + mux
    vrand <- rnorm(N)*sdv
    datp$x <- x
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      xpi <- x[indpi]
      ypi <- cbind(1,xpi)%*%betamdom[i,] + vrand[i] + rnorm(Nv[i])*sde
      datp$y[indpi] <- ypi
      
      indi <- dats$dom == doms[i]
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
    }
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    meansx <- ddply(dats,.(dom),summarize, meansx = mean(x)) ## sample mean x
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    
    estbhf <- 0
    tryCatch({estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                                 method = "REML",data = dats)},
             error=function(e){cat("ERROR :",error(e), "\n")})
    
    # equal weight ###
    c1 <- rep(1,nobs*(nobs-1)/2)
    lam1 <- seq(0.1,1.5,by = 0.05)
    bic1 <- rep(0,length(lam1))
    beta_array1 <- array(0, dim = c(nobs,2,length(lam1)))
    group1 <- matrix(0,nobs,length(lam1))
    sig2mat <- matrix(0, length(lam1),2)
    loglikv <- rep(0,length(lam1))
    
    betam01 <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
    for(l in 1:length(lam1))
    {
      res1 <- Spgr_bhf3(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam01,lam = lam1[l],maxiter = 2000)
      betam01 <- res1$beta
      beta_array1[,,l] <- res1$beta
      bic1[l] <- BICc_bhf3(obj = res1,c0 = c0)
      group1[,l] <- res1$group
      sig2mat[l,] <- res1$sig2est
      loglikv[l] <- res1$loglikvalue
    }
    betaest1 <- beta_array1[,,which.min(bic1)]
    sig2est1 <- sig2mat[which.min(bic1),]
    groupmat1[,m] <- group1[,which.min(bic1)]
    
    gammav = sig2est1[1]/(sig2est1[1] + sig2est1[2]/nv)
    vhat <- gammav*(meansy$meansy - rowSums(cbind(1, meansx$meansx)*betaest1))
    reg11 <- rowSums(cbind(1, meanx$meanx)*betaest1) + vhat
    reg12 <- meansy$meansy*fn + reg11 - fn*rowSums(cbind(1, meansx$meansx)*betaest1)-fn*vhat
   
    ### spatial weight ####
    
    lam2 <- seq(0.1,1.5,by =0.05)
    alp <- c(0.5,1,1.5)
    bic2 <- rep(0, length(alp))
    beta_array2 <- array(0, dim = c(nobs,2,length(alp)))
    group2 <- matrix(0, nobs, length(alp))
    sig2mat2 <- matrix(0, length(alp), 2)
    
    for(l1 in 1:length(alp))
    {
      bic2l <- rep(0,length(lam2))
      beta_array2l <- array(0, dim = c(nobs,2,length(lam2)))
      c2l <- exp(alp[l1]*(1 - ordermat))
      group2l <- matrix(0,nobs, length(lam2))
      sig2mat2l <- matrix(0, length(lam2),2)
      loglikv <- rep(0, length(lam2))
      
      betam02l <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
      betam02l <- cal_initialrx2(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),K0 = 20)
      for(l2 in 1:length(lam2))
      {
        res2l <- Spgr_bhf3(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights  = c2l,betam0 = betam02l,lam = lam2[l2],maxiter = 1000)
        betam02l <- res2l$beta
        beta_array2l[,,l2] <- res2l$beta
        bic2l[l2] <- BICc_bhf3(obj = res2l,c0 = c0)
        group2l[,l2] <- res2l$group
        sig2mat2l[l2,] <- res2l$sig2est
        loglikv[l2] <- res2l$loglikvalue
      }
      
      bic2[l1] <- min(bic2l)
      beta_array2[,,l1] <- beta_array2l[,,which.min(bic2l)]
      group2[,l1] <- group2l[,which.min(bic2l)]
      sig2mat2[l1,] <- sig2mat2l[which.min(bic2l),]
      
    }
    
    betaest2 <- beta_array2[,,which.min(bic2)]
    sig2est2 <- sig2mat2[which.min(bic2),]
    groupmat2[,m] <- group2[,which.min(bic2)]
    
    gammav2 = sig2est2[1]/(sig2est2[1] + sig2est2[2]/nv)
    vhat2 <- gammav2*(meansy$meansy - rowSums(cbind(1, meansx$meansx)*betaest2))
    reg21 <- rowSums(cbind(1, meanx$meanx)*betaest2) + vhat2
    reg22 <- meansy$meansy*fn + reg21 - fn*rowSums(cbind(1, meansx$meansx)*betaest2)-fn*vhat2
    
    if(!is.na(estbhf$eblup)[1])
    {
      nrnzero <- nrnzero + 1
      msemat <-  msemat + cbind((estdir$Direct - meany[,2])^2,
                                (estbhf$eblup$eblup - meany[,2])^2,
                                (reg11 - meany[,2])^2,
                                (reg12 - meany[,2])^2,
                                (reg21 - meany[,2])^2,
                                (reg22 - meany[,2])^2)

    }
    
    
    print(m) 
  }
  
  msemat <- msemat/nrnzero
  return(list(mse = msemat, nrnzero = nrnzero, group1 = groupmat1, group2 = groupmat2))
}

# BIC 
sim_pop_bhf3 <- function(beta, mux, sdx, sdv, sde, rate = 0.1, ni = NULL, M, popn, group,
                         alp = c(0.5,1,1.5),seed = 7231625)
{
  nobs <- nrow(nsegtotal)
  msemat <- matrix(0, nobs, 6)
  colnames(msemat) <- c("direct","bhf","equal1","equal2","sp1","sp2")
  
  doms <- popn[,1]
  Nv <- popn[,2]
  Nt <- sum(Nv)
  N <- nrow(popn)
  
  ng <- nrow(beta)
  betamdom <- beta[group,]
  
  #### population matrix
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(doms, Nv)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  
  ## sample 
  
  if(is.null(rate) & is.null(ni))
  {
    print("specify rate or sample size")
  }
  
  if(!is.null(rate) & is.null(ni))
  {
    nv <- round(Nv * rate)
  }
  
  if(is.null(rate) & !is.null(ni))
  {
    nv <- rep(ni, N)
  }
  
  n0 <- sum(nv) ## total sample size
  fn <- nv/Nv
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x","weights")
  dats$dom <- rep(doms, nv)
  dats$weights <-  rep(Nv/nv, nv)
  
  nrnzero <- 0 
  groupmat1 <- groupmat2 <- matrix(0, nobs, M)
  
  for(m in 1:M)
  {
    set.seed(seed + m)
    x <- rnorm(Nt)*sdx + mux
    vrand <- rnorm(N)*sdv
    datp$x <- x
    
    for(i in 1:N)
    {
      indpi <- indp == doms[i]
      xpi <- x[indpi]
      ypi <- cbind(1,xpi)%*%betamdom[i,] + vrand[i] + rnorm(Nv[i])*sde
      datp$y[indpi] <- ypi
      
      indi <- dats$dom == doms[i]
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
    }
    
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
    
    meansy <- ddply(dats,.(dom),summarize, meansy = mean(y)) ## sample mean y
    meansx <- ddply(dats,.(dom),summarize, meansx = mean(x)) ## sample mean x
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    
    estbhf <- 0
    tryCatch({estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                                 method = "REML",data = dats)},
             error=function(e){cat("ERROR :",error(e), "\n")})
    
    # equal weight ###
    c1 <- rep(1,nobs*(nobs-1)/2)
    lam1 <- seq(0.1,1.5,by = 0.05)
    bic1 <- rep(0,length(lam1))
    beta_array1 <- array(0, dim = c(nobs,2,length(lam1)))
    group1 <- matrix(0,nobs,length(lam1))
    sig2mat <- matrix(0, length(lam1),2)
    loglikv <- rep(0,length(lam1))
    
    betam01 <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
    for(l in 1:length(lam1))
    {
      res1 <- Spgr_bhf3(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam01,lam = lam1[l],maxiter = 2000)
      betam01 <- res1$beta
      beta_array1[,,l] <- res1$beta
      bic1[l] <- BIC_bhf3(obj = res1)
      group1[,l] <- res1$group
      sig2mat[l,] <- res1$sig2est
      loglikv[l] <- res1$loglikvalue
    }
    betaest1 <- beta_array1[,,which.min(bic1)]
    sig2est1 <- sig2mat[which.min(bic1),]
    groupmat1[,m] <- group1[,which.min(bic1)]
    
    gammav = sig2est1[1]/(sig2est1[1] + sig2est1[2]/nv)
    vhat <- gammav*(meansy$meansy - rowSums(cbind(1, meansx$meansx)*betaest1))
    reg11 <- rowSums(cbind(1, meanx$meanx)*betaest1) + vhat
    reg12 <- meansy$meansy*fn + reg11 - fn*rowSums(cbind(1, meansx$meansx)*betaest1)-fn*vhat
    
    ### spatial weight ####
    
    lam2 <- seq(0.1,1.5,by =0.05)
    alp <- c(0.5,1,1.5)
    bic2 <- rep(0, length(alp))
    beta_array2 <- array(0, dim = c(nobs,2,length(alp)))
    group2 <- matrix(0, nobs, length(alp))
    sig2mat2 <- matrix(0, length(alp), 2)
    
    for(l1 in 1:length(alp))
    {
      bic2l <- rep(0,length(lam2))
      beta_array2l <- array(0, dim = c(nobs,2,length(lam2)))
      c2l <- exp(alp[l1]*(1 - ordermat))
      group2l <- matrix(0,nobs, length(lam2))
      sig2mat2l <- matrix(0, length(lam2),2)
      loglikv <- rep(0, length(lam2))
      
      betam02l <- cal_initialrx(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x))
      betam02l <- cal_initialrx2(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),K0 = 20)
      for(l2 in 1:length(lam2))
      {
        res2l <- Spgr_bhf3(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights  = c2l,betam0 = betam02l,lam = lam2[l2],maxiter = 1000)
        betam02l <- res2l$beta
        beta_array2l[,,l2] <- res2l$beta
        bic2l[l2] <- BIC_bhf3(obj = res2l)
        group2l[,l2] <- res2l$group
        sig2mat2l[l2,] <- res2l$sig2est
        loglikv[l2] <- res2l$loglikvalue
      }
      
      bic2[l1] <- min(bic2l)
      beta_array2[,,l1] <- beta_array2l[,,which.min(bic2l)]
      group2[,l1] <- group2l[,which.min(bic2l)]
      sig2mat2[l1,] <- sig2mat2l[which.min(bic2l),]
      
    }
    
    betaest2 <- beta_array2[,,which.min(bic2)]
    sig2est2 <- sig2mat2[which.min(bic2),]
    groupmat2[,m] <- group2[,which.min(bic2)]
    
    gammav2 = sig2est2[1]/(sig2est2[1] + sig2est2[2]/nv)
    vhat2 <- gammav2*(meansy$meansy - rowSums(cbind(1, meansx$meansx)*betaest2))
    reg21 <- rowSums(cbind(1, meanx$meanx)*betaest2) + vhat2
    reg22 <- meansy$meansy*fn + reg21 - fn*rowSums(cbind(1, meansx$meansx)*betaest2)-fn*vhat2
    
    if(!is.na(estbhf$eblup)[1])
    {
      nrnzero <- nrnzero + 1
      msemat <-  msemat + cbind((estdir$Direct - meany[,2])^2,
                                (estbhf$eblup$eblup - meany[,2])^2,
                                (reg11 - meany[,2])^2,
                                (reg12 - meany[,2])^2,
                                (reg21 - meany[,2])^2,
                                (reg22 - meany[,2])^2)
      
    }
    
    
    print(m) 
  }
  
  msemat <- msemat/nrnzero
  return(list(mse = msemat, nrnzero = nrnzero, group1 = groupmat1, group2 = groupmat2))
}

betad15 <-  matrix(c(0.5,0.5,2,2,3.5,3.5),ncol=2,byrow = TRUE)
beta <- matrix(c(0.5,0.5,1.5,1.5,2.5,2.5),ncol=2,byrow = TRUE)
beta1 <- matrix(c(0.5,0.5,1.25,1.25,2,2),ncol=2,byrow = TRUE)
beta2 <- matrix(c(0.5,0.5,1,1,1.5,1.5),ncol=2,byrow = TRUE)
mux <- 1
sdx <- sdv <- sde <- 1
rate <- 0.01

sim_bhf2_d1_r01_sde1_v2 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup,c0 = 1)

sim_bhf2_d1_r01_sde1_v3 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup,c0 = 0.5)

sim_bhf2_d1_r01_sde1_v4 <- sim_pop_bhf3(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)



save(sim_bhf2_d1_r01_sde1_v2,sim_bhf2_d1_r01_sde1_v3,sim_bhf2_d1_r01_sde1_v4,
     file = "result/sim_bhf2_d1_r01_sde1_bic.RData" )

apply(sim_bhf2_d1_r01_sde1_v2$group2,2,fun1)
apply(sim_bhf2_d1_r01_sde1_v3$group1,2,fun1)


sim_bhf2_d1_r01_sde05_v2 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup,c0 = 1)
sim_bhf2_d1_r01_sde05_v3 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup,c0 = 0.5)

sim_bhf2_d1_r01_sde05_v4 <- sim_pop_bhf3(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)

save(sim_bhf2_d1_r01_sde05_v2,sim_bhf2_d1_r01_sde05_v3,sim_bhf2_d1_r01_sde05_v4,
     file = "result/sim_bhf2_d1_r01_sde05_bic.RData" )


apply(sim_bhf2_d1_r01_sde05_v2$group2,2,fun1)
apply(sim_bhf2_d1_r01_sde05_v3$group2,2,fun1)


sim_bhf2_d1_r01_sde1 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d1_r01_sde2 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 2,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d1_r01_sde05 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)

save(sim_bhf2_d1_r01_sde1, file = "result/sim_bhf2_d1_r01_sde1.RData")
save(sim_bhf2_d1_r01_sde2, file = "result/sim_bhf2_d1_r01_sde2.RData")
save(sim_bhf2_d1_r01_sde05 , file = "result/sim_bhf2_d1_r01_sde05.RData")

sim_bhf2_d075_r01_sde05 <- sim_pop_bhf2(beta = beta1, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d075_r01_sde1 <- sim_pop_bhf2(beta = beta1, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d075_r01_sde2 <- sim_pop_bhf2(beta = beta1, mux = 1,sdx = 1,sdv = 1, sde = 2,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)

save(sim_bhf2_d075_r01_sde05 , file = "result/sim_bhf2_d075_r01_sde05.RData")
save(sim_bhf2_d075_r01_sde1 , file = "result/sim_bhf2_d075_r01_sde1.RData")
save(sim_bhf2_d075_r01_sde2 , file = "result/sim_bhf2_d075_r01_sde2.RData")


boxplot(sqrt(sim_bhf2_d075_r01_sde05$mse))
boxplot(sqrt(sim_bhf2_d075_r01_sde1$mse))
boxplot(sqrt(sim_bhf2_d075_r01_sde2$mse))


sim_bhf2_d05_r01_sde1 <- sim_pop_bhf2(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d05_r01_sde05 <- sim_pop_bhf2(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d05_r01_sde2 <- sim_pop_bhf2(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 2,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)

save(sim_bhf2_d05_r01_sde1 , file = "result/sim_bhf2_d05_r01_sde1.RData")
save(sim_bhf2_d05_r01_sde05 , file = "result/sim_bhf2_d05_r01_sde05.RData")
save(sim_bhf2_d05_r01_sde2 , file = "result/sim_bhf2_d05_r01_sde2.RData")

boxplot(sqrt(sim_bhf2_d05_r01_sde05$mse))
boxplot(sqrt(sim_bhf2_d05_r01_sde1$mse))
boxplot(sqrt(sim_bhf2_d05_r01_sde2$mse))


sim_bhf2_d15_r01_sde05 <- sim_pop_bhf3(beta = betad15, mux = 1,sdx = 1,sdv = 1, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d15_r01_sde05 <- sim_pop_bhf3(beta = betad15, mux = 1,sdx = 1,sdv = 0.5, sde = 0.5,rate = 0.01,M = 100,popn = popn,group = nsegtotal$newgroup)




#### draw figures ####

pdf("doc/figures/rmse_sae_bhf_b1.pdf",width = 10,height = 4)
par(mfrow = c(1,3))
boxplot(sqrt(sim_bhf2_d1_r01_sde05$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=0.5")),las = 2)
boxplot(sqrt(sim_bhf2_d1_r01_sde1$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=1")),las = 2)
boxplot(sqrt(sim_bhf2_d1_r01_sde2$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=2")),las = 2)
dev.off()

pdf("doc/figures/rmse_sae_bhf_b075.pdf",width = 10,height = 4)
par(mfrow = c(1,3))
boxplot(sqrt(sim_bhf2_d075_r01_sde05$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=0.5")),las = 2)
boxplot(sqrt(sim_bhf2_d075_r01_sde1$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=1")),las = 2)
boxplot(sqrt(sim_bhf2_d075_r01_sde2$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=2")),las = 2)
dev.off()

pdf("doc/figures/rmse_sae_bhf_b05.pdf",width = 10,height = 4)
par(mfrow = c(1,3))
boxplot(sqrt(sim_bhf2_d05_r01_sde05$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=0.5")),las = 2)
boxplot(sqrt(sim_bhf2_d05_r01_sde1$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=1")),las = 2)
boxplot(sqrt(sim_bhf2_d05_r01_sde2$mse),ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=2")),las = 2)
dev.off()


fun1 <- function(x){return(length(unique(x)))}


apply(sim_bhf2_d1_r01_sde05$group1,2,fun1)
apply(sim_bhf2_d1_r01_sde05$group2,2,fun1)
apply(sim_bhf2_d1_r01_sde1$group1,2,fun1)
apply(sim_bhf2_d1_r01_sde1$group2,2,fun1)
apply(sim_bhf2_d1_r01_sde2$group1,2,fun1)
apply(sim_bhf2_d1_r01_sde2$group2,2,fun1)



sim_bhf2_d05_r025_sde1 <- sim_pop_bhf2(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.025,M = 100,popn = popn,group = nsegtotal$newgroup)
sim_bhf2_d05_r025_sde2 <- sim_pop_bhf2(beta = beta2, mux = 1,sdx = 1,sdv = 1, sde = 2,rate = 0.025,M = 100,popn = popn,group = nsegtotal$newgroup)

boxplot(sqrt(sim_bhf2_d05_r025_sde2$mse))
save(sim_bhf2_d05_r025_sde1, file = "result/sim_bhf2_d05_r025_sde1.RData")

sim_bhf2_d075_r025_sde1 <- sim_pop_bhf2(beta = beta1, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.025,M = 100,popn = popn,group = nsegtotal$newgroup)

save(sim_bhf2_d075_r025_sde1 , file = "result/sim_bhf2_d075_r025_sde1.RData")

sim_bhf2_d1_r025_sde1 <- sim_pop_bhf2(beta = beta, mux = 1,sdx = 1,sdv = 1, sde = 1,rate = 0.025,M = 100,popn = popn,group = nsegtotal$newgroup)

save(sim_bhf2_d1_r025_sde1 , file = "result/sim_bhf2_d1_r025_sde1.RData")

boxplot(sqrt(sim_bhf2_d1_r025_sde1$mse))
boxplot(sqrt(sim_bhf2_d1_r01_sde1$mse))


pdf("doc/figures/rmse_sae_bhf_b05_p2.pdf",height = 4,width = 8)
par(mfrow = c(1,2))
boxplot(sim_bhf2_d05_r01_sde1$mse,ylim = c(0,0.15), ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=1")),las = 2)
boxplot(sim_bhf2_d05_r025_sde1$mse,ylim = c(0,0.15),ylab = "RMSE",main = expression(paste("r=0.025, ",sigma[epsilon],"=1")),las = 2)
dev.off()

pdf("doc/figures/rmse_sae_bhf_b075_p2.pdf",height = 4,width = 8)
par(mfrow = c(1,2))
boxplot(sim_bhf2_d075_r01_sde1$mse,ylim = c(0,0.25), ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=1")),las = 2)
boxplot(sim_bhf2_d075_r025_sde1$mse,ylim = c(0,0.25),ylab = "RMSE",main = expression(paste("r=0.025, ",sigma[epsilon],"=1")),las = 2)
dev.off()

pdf("doc/figures/rmse_sae_bhf_b1_p2.pdf",height = 4,width = 8)
par(mfrow = c(1,2))
boxplot(sim_bhf2_d1_r01_sde1$mse,ylim = c(0,0.25), ylab = "RMSE",main = expression(paste("r=0.01, ",sigma[epsilon],"=1")),las = 2)
boxplot(sim_bhf2_d1_r025_sde1$mse,ylim = c(0,0.25),ylab = "RMSE",main = expression(paste("r=0.025, ",sigma[epsilon],"=1")),las = 2)
dev.off()


##

par(mfrow= c(3,1))
plot(dats$x, dats$y, cex = 0.5,main = expression(paste(sigma[epsilon],"=0.5")))
plot(dats$x, dats$y, cex=0.5,main = expression(paste(sigma[epsilon],"=1")))
plot(dats$x, dats$y, cex=0.5,main = expression(paste(sigma[epsilon],"=2")))
