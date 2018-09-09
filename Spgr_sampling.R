# without spatial group structure  only x included different intercept and different slope 
# weighted least square 

# maxiter <- 1000
# nu <- 1
# gam <- 3
# lam <- 0.4
# tolabs = 1e-4
# tolrel = 1e-2
# cvec <- rep(1, n*(n-1)/2)

Spgr_sampling <- function(yhat, x, sde, cvec, betam0,
                          nu = 1, gam = 3, lam = 0.5,
                          maxiter = 1000, tolabs = 1e-4, tolrel = 1e-2)
{
  nobs <- length(yhat)
  n1 <- nobs*(nobs - 1)/2
  ncx <- ncol(x)
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
  
  Xm <- matrix(0, nobs, nobs*ncx) ## x to large matrix 
  for(i in 1:nobs)
  {
    Xm[i,(ncx*(i-1) + 1) : (ncx*i)] <- x[i,]
  }
  
  #### initial values 
  deltae <- sde^2

  #Omega <- Omega0
  vm  <-  matrix(0, ncx, nobs*(nobs-1)/2)
  deltam.old <- t(D %*% betam0)
  ###########
  sig2 <- mean((yhat - rowSums(x*betam0))^2)
  
  Omega <- 1/(deltae + sig2)
  
  flag <- 0
  for(m in 1:maxiter)
  {
    XtX <- t(Omega*Xm)%*%Xm
    Xinv <- solve(XtX + nu*AtA)
    beta.new <-  Xinv %*% (t(Xm)%*%(Omega*yhat) + nu * c((deltam.old -  vm/nu) %*% D))
    betam <- matrix(beta.new, nobs, ncx, byrow = TRUE)
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    
    groupest <- getgroup(deltam = deltam,n = nobs)
    ngest <- length(unique(groupest))
    
    Ivalue<- sum(1/(sig2 + deltae)^2)
    svalue <- sum((yhat - rowSums(x*betam))^2/(sig2 + deltae)^2) - sum(1/(sig2 + deltae))
    sig2 <- sig2 + svalue/Ivalue
    if(sig2 < 0) {sig2 <- 0}
    
    
    ### moment ##
    # avalue <- sum((yhat - Xm%*%beta.new)^2/(sig2 + deltae))
    # aprime <- -sum((yhat - Xm%*%beta.new)^2/(sig2 + deltae)^2) 
    # sig2 <- sig2 + (n - ncx*ngest - avalue)/aprime
    # if(sig2 < 0) {sig2 <- 0}
    
    Omega <- 1/(deltae + sig2)
    
    
    tolpri <- tolabs*sqrt(n1*ncx) + tolrel*max(sqrt(sum(betadiff^2)),sqrt(sum(deltam^2)))
    toldual <- tolabs*sqrt(nobs*ncx) + tolrel* sqrt(sum((vm %*% D)^2))
    
    rm <- sqrt(sum((betadiff - deltam)^2))
    sm <- nu*sqrt(sum(((deltam - deltam.old)%*%D)^2))
    
    deltam.old <- deltam
    
    if(rm <= tolpri & sm <= toldual){break}
    
  }
  
  if(m == maxiter) {flag <- 1}
  
  alpest <- do.call("cbind",by(betam, groupest, colMeans,simplify = TRUE))
  betaest <- t(alpest[,groupest])
  
  loglikvalue <- -0.5*sum(log(deltae + sig2)) - 0.5* sum((yhat - rowSums(x*betaest))^2/(deltae + sig2))
  
  outls <- list(betaest = betaest, betam = betam, sig2 = sig2, group = groupest,
                deltam = deltam, flag = flag,loglikvalue = loglikvalue,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                niteration = m)
  return(outls)
}
