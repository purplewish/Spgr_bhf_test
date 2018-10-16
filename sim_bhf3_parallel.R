#### subgroup_sae_unit #### parallel with tradtional BIC##
library(sae)
library(Spgr)
library(Rcpp)
library(plyr)
library(Matrix)
library(doParallel) 
library(foreach)
library(SpgrBHF)
#sourceCpp("code/Spgr_bhf3.cpp")

load("data/Cmat_popn.RData")
ordermat <- getorder(Matrix(Cmatia))

# BIC 
sim_bhf3_parallel <- function(beta, mux, sdx, sdv, sde, rate = 0.1, ni = NULL, M = 1,
                              popn, group,alp = c(0.5,1,1.5),seed = 7231625)
{
  nobs <- nrow(popn)
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
    
  }
  
  msemat <- msemat/nrnzero
  return(list(mse = msemat, nrnzero = nrnzero, group1 = groupmat1, group2 = groupmat2))
}

subfun <- function(mm)
{
  sim_bhf3_parallel(beta = betad15,sdx = 1,mux = 1,sdv = 1,sde = 1,rate = 0.01,M = 1,popn = popn[,1:2],group = popn[,3],seed = mm + 82245)
}

betad15 <-  matrix(c(0.5,0.5,2,2,3.5,3.5),ncol=2,byrow = TRUE)

cl <- makeCluster(2)  
registerDoParallel(cl)  
sim3 <- foreach(mm=1:2,.packages=c("plyr","sae","Spgr","SpgrBHF")) %dopar% subfun(mm)
stopCluster(cl) 
save(sim3, file = "test.RData")



