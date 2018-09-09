# em with phi 
###### simulation to two level model ###
# without spatial group structure  only x included different intercept and different slope 
# EM algorithm 
# maxiter <- 1000
# nu <- 1
# gam <- 3
# lam <- 0.4
# tolabs = 1e-4
# tolrel = 1e-2
# cvec <- rep(1, n*(n-1)/2)



Spgr_sampling_emp <- function(yhat, x, sde, cvec, betam0, phi = 1,
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
  vm  <-  matrix(0, ncx, nobs*(nobs-1)/2)
  deltam.old <- t(D %*% betam0)
  betam <- betam0
  ###########
  sig2 <- mean((yhat - rowSums(x*betam))^2)
  
  flag <- 0
  for(m in 1:maxiter)
  {
    # XtX <- t(Omega*Xm)%*%Xm
    # Xinv <- solve(XtX + nu*AtA)
    # beta.new <-  Xinv %*% (t(Xm)%*%(Omega*yhat) + nu * c((deltam.old -  vm/nu) %*% D))
    
    muest <- rowSums(x * betam)
    mhat <- (sig2*yhat + phi*deltae*muest)/(sig2 + phi*deltae)
    Vhat <- phi*deltae*sig2/(sig2 + phi*deltae)
    
    XtX <- t(Xm)%*%Xm
    Xinv <- solve(XtX + nu*AtA)
    beta.new <-  Xinv %*% (t(Xm)%*%mhat + nu * c((deltam.old -  vm/nu) %*% D))
    
    betam <- matrix(beta.new, nobs, ncx, byrow = TRUE)
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    
    groupest <- getgroup(deltam = deltam,n = nobs)
    ngest <- length(unique(groupest))
    
    sig2 <- mean((mhat - rowSums(x * betam))^2 + Vhat)
    phi <- mean(((yhat - mhat)^2 + Vhat)/deltae)
    
    
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
  
  # marginal likelihood 
  loglikvalue <- -0.5*sum(log(deltae + sig2)) - 0.5* sum((yhat - rowSums(x*betaest))^2/(deltae + sig2))
  
  outls <- list(betaest = betaest, betam = betam, sig2 = sig2, phi = phi,
                group = groupest, mhat = mhat,
                deltam = deltam, flag = flag,loglikvalue = loglikvalue,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                niteration = m)
  return(outls)
}