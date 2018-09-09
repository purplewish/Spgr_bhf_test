###### subgroup analysis with area random effect #####
# scad function is used ####
source("/Users/Xin/Research/cluster/code/mcp_scad.R")

Spgr_bhf <- function(indexy, y, x,cvec, nu = 1, gam = 3, lam = 0.5,
                     betam0,maxiter = 500, tolabs = 1e-4, tolrel = 1e-2)
{
  n0 <- length(y)
  ncx <- ncol(x)
  
  uniqxy <- unique(indexy)
  nobs <- length(uniqxy)
  npair <- nobs*(nobs-1)/2
  
  Ip <- diag(1,ncx,ncx)
  
  D <- matrix(0,nobs*(nobs-1)/2,nobs)
  for(j in 1:(nobs-1))
  {
    indexj <- (nobs-1 + nobs-j+1)*(j-1)/2
    indexvj <- indexj + (1:(nobs-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nobs)] <- -1
  }
  
  AtA <- t(D)%*%D %x% Ip
  
  Xm <- matrix(0, n0, nobs*ncx)
  nJ <- rep(0, nobs)
  for(i in 1:nobs)
  {
    nJ[i] <- sum(indexy == uniqxy[i])
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }
  
  matls <- list()
  diagls <- list()
  
  for(i in 1:nobs)
  {
    matls[[i]] <-  matrix(1,ncol=1,nrow = nJ[i])
    diagls[[i]] <- diag(1, nJ[i],nJ[i])
  }
  
  vm  <-  matrix(0, ncx, nobs*(nobs-1)/2)
  deltam.old <- t(D %*% betam0)
  betam <- betam0
  betanew <- c(t(betam))
  
  ## two functions used to update 
  dfun <- function(par, betam)
  {
    sigmav2 <- par[1]
    sigmae2 <- par[2]
    
    d11 <- d12 <- 0
    for(i in 1:nobs)
    {
      betami <- betam[i,]
      n1 <- nJ[i]
      mat1 <- matls[[i]]
      diag1 <- diagls[[i]]
      
      M1 <- - mat1%*%t(mat1)*(1/(sigmae2 + n1*sigmav2)^2)
      M2 <- (sigmav2*(2*sigmae2 + n1*sigmav2)/(sigmae2+ n1*sigmav2)^2* mat1%*%t(mat1) - diag1)/(sigmae2)^2
      
      indexi = indexy == uniqxy[i]
      
      d11 <- d11  -0.5*n1/(sigmae2+ n1*sigmav2) - 0.5*t(y[indexi] - x[indexi,]%*%betami)%*%M1%*%(y[indexi] - x[indexi,]%*%betami)
      d12 <- d12 - 0.5*n1*(1- sigmav2/(sigmae2 +n1*sigmav2))/sigmae2 - 0.5*t(y[indexi] - x[indexi,]%*%betami)%*%M2%*%(y[indexi] - x[indexi,]%*%betami)
    }
    return(c(d11,d12))
  }
  
  Ifun <- function(par)
  {
    sigmav2 <- par[1]
    sigmae2 <- par[2]
    
    I11 <- 0.5*sum(nJ^2/(sigmae2 + nJ*sigmav2)^2)
    I12 <- 0.5*sum(nJ/(sigmae2 + nJ*sigmav2)^2)
    I22 <- 0.5*sum(nJ*(1- sigmav2*(2*sigmae2 + nJ*sigmav2)/(sigmae2 + nJ*sigmav2)^2)/sigmae2^2)
    Imat1 <- matrix(0,2,2)
    Imat1[1,1] <- I11
    Imat1[1,2] <- Imat1[2,1] <- I12
    Imat1[2,2] <- I22
    return(Imat1)
  }
  
  sige2 <- mean((y - Xm %*% betanew)^2)
  siga2 <- 0 
  
  sig2est <- c(siga2, sige2)
  
  flag <- 0
  for(m in 1:maxiter)
  {
    XtX <- matrix(0, nobs*ncx, nobs*ncx)
    Xty <- rep(0, nobs*ncx)
    for(i in 1:nobs)
    {
      n1 <- nJ[i]
      mat1 <- matls[[i]]
      diag1 <- diagls[[i]]
      S1 <- mat1%*%t(mat1)*siga2 + diag1*sige2
      S1inv <- solve(S1)
      
      xi <- x[indexy == uniqxy[i],]
      
      index1 <- (i-1)*ncx + 1
      index2 <- i*ncx
      
      XtX[index1:index2,index1:index2] <- t(xi)%*%S1inv%*%xi
      Xty[index1:index2] <- t(xi)%*%S1inv%*%y[indexy == uniqxy[i]]
    }
    
    betanew <- solve(XtX + nu*AtA)%*%(Xty + nu * c((deltam.old -  vm/nu) %*% D))
    betam <- matrix(betanew, nobs, ncx, byrow = TRUE)
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    
    groupest <- getgroup(deltam = deltam,n = nobs)
    ngest <- length(unique(groupest))
    
    
    # sig2est <- sig2est + solve(Ifun(sig2est))%*%dfun(sig2est,betam = betam)
    # sig2est[sig2est <0] <- 0 
    

    for(j in 1:10)
    {
      sig2est <- sig2est + solve(Ifun(sig2est))%*%dfun(sig2est, betam = betam)
      sig2est[sig2est <0] <- 0 
    }
    
    #sig2est[sig2est <0] <- 0 
    
    sige2 <- sig2est[2]
    siga2 <- sig2est[1] 
  
    tolpri <- tolabs*sqrt(npair*ncx) + tolrel*max(sqrt(sum(betadiff^2)),sqrt(sum(deltam^2)))
    toldual <- tolabs*sqrt(nobs*ncx) + tolrel* sqrt(sum((vm %*% D)^2))
    
    rm <- sqrt(sum((betadiff - deltam)^2))
    sm <- nu*sqrt(sum(((deltam - deltam.old)%*%D)^2))
    
    deltam.old <- deltam
    
    if(rm <= tolpri & sm <= toldual){break}
  }
  
  if(m == maxiter) {flag <- 1}

  rownames(sig2est) <- c("random","error")
  
  alpest <- do.call("cbind",by(betam, groupest, colMeans,simplify = TRUE))
  betaest <- t(alpest[,groupest])
  
  ## likelihood value ###
  loglikvalue <- 0
  for(i in 1:nobs)
  {
    n1 <- nJ[i]
    mat1 <- matls[[i]]
    diag1 <- diagls[[i]]
    
    S1 <- mat1%*%t(mat1)*siga2 + diag1*sige2
    S1inv <- solve(S1)
    
    indexi = indexy == uniqxy[i]
    loglikvalue <- loglikvalue - 0.5*sum(log(eigen(S1)$values)) - 0.5*t(y[indexi] - x[indexi,]%*%betaest[i,])%*%S1inv%*%(y[indexi] - x[indexi,]%*%betaest[i,])
  }
  
  
  outls <- list(betaest = betaest, betam = betam, sig2est = sig2est, group = groupest,
                deltam = deltam, flag = flag,loglikvalue = loglikvalue,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                niteration = m)
  return(outls)
  
}

# indexy <- rep(1:length(nv), nv)
# betam0 <- cal_initialrx(indexy = indexy,y = yv,x = xm)
# 
# res1 <- Spgr_bhf(indexy = indexy, y = yv,x = xm,cvec = rep(1,length(nv)*(length(nv)-1)/2),lam = 0.85,betam0 = betam0,maxiter = 2000)
# 
# res2 <- Spgr_bhf(indexy = indexy, y = yv,x = xm,cvec = rep(1,length(nv)*(length(nv)-1)/2),lam = 1,betam0 = betam0,maxiter = 2000)
# 
# res3 <- Spgr_bhf(indexy = indexy, y = yv,x = xm,cvec = rep(1,length(nv)*(length(nv)-1)/2),lam = 3,betam0 = betam0,maxiter = 2000)

#BIC is defined as -2loglikevalue + Cn log(nobs)(K*ncx), Cn  = c0*log(log(nobs*ncx))

BICc_bhf <- function(obj, c0=0.2)
{
  ncx <- ncol(obj$betam)
  nobs <- nrow(obj$betam)
  ngest <- length(unique(obj$group))
  Cn <- c0*log(log(nobs*ncx))
  
  bicvalue <- -2*obj$loglikvalue + Cn*log(nobs)*(ngest*ncx)
  return(bicvalue)
  
}

# BICc_bhf(res1,0.1)
# BICc_bhf(res2,0.1)
# BICc_bhf(res3,0.1)
# 
# 
#
#### weighted 
Spgr_bhfw <- function(indexy, y, x,cvec, nu = 1, gam = 3, lam = 0.5,
betam0,maxiter = 500, tolabs = 1e-4, tolrel = 1e-2)
{
  n0 <- length(y)
  ncx <- ncol(x)
  
  uniqxy <- unique(indexy)
  nobs <- length(uniqxy)
  npair <- nobs*(nobs-1)/2
  
  Ip <- diag(1,ncx,ncx)
  
  D <- matrix(0,nobs*(nobs-1)/2,nobs)
  for(j in 1:(nobs-1))
  {
    indexj <- (nobs-1 + nobs-j+1)*(j-1)/2
    indexvj <- indexj + (1:(nobs-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nobs)] <- -1
  }
  
  AtA <- t(D)%*%D %x% Ip
  
  Xm <- matrix(0, n0, nobs*ncx)
  nJ <- rep(0, nobs)
  for(i in 1:nobs)
  {
    nJ[i] <- sum(indexy == uniqxy[i])
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }
  
  matls <- list()
  diagls <- list()
  
  for(i in 1:nobs)
  {
    matls[[i]] <-  matrix(1,ncol=1,nrow = nJ[i])
    diagls[[i]] <- diag(1, nJ[i],nJ[i])
  }
  
  vm  <-  matrix(0, ncx, nobs*(nobs-1)/2)
  deltam.old <- t(D %*% betam0)
  betam <- betam0
  betanew <- c(t(betam))
  
  ## two functions used to update 
  dfun <- function(par, betam)
  {
    sigmav2 <- par[1]
    sigmae2 <- par[2]
    
    d11 <- d12 <- 0
    for(i in 1:nobs)
    {
      betami <- betam[i,]
      n1 <- nJ[i]
      mat1 <- matls[[i]]
      diag1 <- diagls[[i]]
      
      M1 <- - mat1%*%t(mat1)*(1/(sigmae2 + n1*sigmav2)^2)
      M2 <- (sigmav2*(2*sigmae2 + n1*sigmav2)/(sigmae2+ n1*sigmav2)^2* mat1%*%t(mat1) - diag1)/(sigmae2)^2
      
      indexi = indexy == uniqxy[i]
      
      d11 <- d11  -0.5/(sigmae2+ n1*sigmav2) - 0.5*t(y[indexi] - x[indexi,]%*%betami)%*%M1%*%(y[indexi] - x[indexi,]%*%betami)/n1
      d12 <- d12 - 0.5*(1- sigmav2/(sigmae2 +n1*sigmav2))/sigmae2 - 0.5*t(y[indexi] - x[indexi,]%*%betami)%*%M2%*%(y[indexi] - x[indexi,]%*%betami)/n1
    }
    return(c(d11,d12))
  }
  
  Ifun <- function(par)
  {
    sigmav2 <- par[1]
    sigmae2 <- par[2]
    
    I11 <- 0.5*sum(nJ/(sigmae2 + nJ*sigmav2)^2)
    I12 <- 0.5*sum(1/(sigmae2 + nJ*sigmav2)^2)
    I22 <- 0.5*sum((1- sigmav2*(2*sigmae2 + nJ*sigmav2)/(sigmae2 + nJ*sigmav2)^2)/sigmae2^2)
    Imat1 <- matrix(0,2,2)
    Imat1[1,1] <- I11
    Imat1[1,2] <- Imat1[2,1] <- I12
    Imat1[2,2] <- I22
    return(Imat1)
  }
  
  sige2 <- mean((y - Xm %*% betanew)^2)
  siga2 <- 0 
  
  sig2est <- c(siga2, sige2)
  
  flag <- 0
  for(m in 1:maxiter)
  {
    XtX <- matrix(0, nobs*ncx, nobs*ncx)
    Xty <- rep(0, nobs*ncx)
    for(i in 1:nobs)
    {
      n1 <- nJ[i]
      mat1 <- matls[[i]]
      diag1 <- diagls[[i]]
      S1 <- mat1%*%t(mat1)*siga2 + diag1*sige2
      S1inv <- solve(S1)/n1
      
      xi <- x[indexy == uniqxy[i],]
      
      index1 <- (i-1)*ncx + 1
      index2 <- i*ncx
      
      XtX[index1:index2,index1:index2] <- t(xi)%*%S1inv%*%xi
      Xty[index1:index2] <- t(xi)%*%S1inv%*%y[indexy == uniqxy[i]]
    }
    
    betanew <- solve(XtX + nu*AtA)%*%(Xty + nu * c((deltam.old -  vm/nu) %*% D))
    betam <- matrix(betanew, nobs, ncx, byrow = TRUE)
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    
    groupest <- getgroup(deltam = deltam,n = nobs)
    ngest <- length(unique(groupest))
    
    
    # sig2est <- sig2est + solve(Ifun(sig2est))%*%dfun(sig2est,betam = betam)
    # sig2est[sig2est <0] <- 0 
    
    
    for(j in 1:10)
    {
      sig2est <- sig2est + solve(Ifun(sig2est))%*%dfun(sig2est, betam = betam)
      sig2est[sig2est <0] <- 0 
    }
    
    #sig2est[sig2est <0] <- 0 
    
    sige2 <- sig2est[2]
    siga2 <- sig2est[1] 
    
    tolpri <- tolabs*sqrt(npair*ncx) + tolrel*max(sqrt(sum(betadiff^2)),sqrt(sum(deltam^2)))
    toldual <- tolabs*sqrt(nobs*ncx) + tolrel* sqrt(sum((vm %*% D)^2))
    
    rm <- sqrt(sum((betadiff - deltam)^2))
    sm <- nu*sqrt(sum(((deltam - deltam.old)%*%D)^2))
    
    deltam.old <- deltam
    
    if(rm <= tolpri & sm <= toldual){break}
  }
  
  if(m == maxiter) {flag <- 1}
  
  rownames(sig2est) <- c("random","error")
  
  alpest <- do.call("cbind",by(betam, groupest, colMeans,simplify = TRUE))
  betaest <- t(alpest[,groupest])
  
  ## likelihood value ###
  loglikvalue <- 0
  for(i in 1:nobs)
  {
    n1 <- nJ[i]
    mat1 <- matls[[i]]
    diag1 <- diagls[[i]]
    
    S1 <- mat1%*%t(mat1)*siga2 + diag1*sige2
    S1inv <- solve(S1)/n1
    
    indexi = indexy == uniqxy[i]
    loglikvalue <- loglikvalue - 0.5*sum(log(eigen(S1)$values))/n1 - 0.5*t(y[indexi] - x[indexi,]%*%betaest[i,])%*%S1inv%*%(y[indexi] - x[indexi,]%*%betaest[i,])
  }
  
  
  outls <- list(betaest = betaest, betam = betam, sig2est = sig2est, group = groupest,
                deltam = deltam, flag = flag,loglikvalue = loglikvalue,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                niteration = m)
  return(outls)
  
}
