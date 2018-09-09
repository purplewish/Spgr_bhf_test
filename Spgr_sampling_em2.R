
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


Spgr_sampling_em2 <- function(yhat, z, x, sde, cvec, betam0,
                             nu = 1, gam = 3, lam = 0.5,
                             maxiter = 1000, tolabs = 1e-4, tolrel = 1e-2)
{
  nobs <- length(yhat)
  n1 <- nobs*(nobs - 1)/2
  ncx <- ncol(x)
  ncz <- ncol(z)
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
  ztz <- solve(t(z) %*% z)%*%t(z)
  Qz <- diag(1, nobs, nobs) - z%*% ztz
  
  etaest <- ztz %*% (yhat - rowSums(x*betam))
  sig2 <- mean((yhat - z%*%etaest - rowSums(x*betam))^2)
  
  flag <- 0
  for(m in 1:maxiter)
  {
    
    muest <- z%*%etaest + rowSums(x * betam)
    mhat <- (sig2*yhat + deltae*muest)/(sig2 + deltae)
    Vhat <- deltae*sig2/(sig2 + deltae)
    
    XtX <- t(Xm)%*%Qz%*%Xm
    Xinv <- solve(XtX + nu*AtA)
    beta.new <-  Xinv %*% (t(Xm)%*%Qz%*%mhat + nu * c((deltam.old -  vm/nu) %*% D))
    
    betam <- matrix(beta.new, nobs, ncx, byrow = TRUE)
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    
    if(is.numeric(deltam)){deltam = matrix(deltam, nrow=1)}
    
    etaest <- ztz %*% (yhat - rowSums(x*betam))
    
    sig2 <- mean((mhat - z%*%etaest - rowSums(x * betam))^2 + Vhat)
    
    
    tolpri <- tolabs*sqrt(n1*ncx) + tolrel*max(sqrt(sum(betadiff^2)),sqrt(sum(deltam^2)))
    toldual <- tolabs*sqrt(nobs*ncx) + tolrel* sqrt(sum((vm %*% D)^2))
    
    rm <- sqrt(sum((betadiff - deltam)^2))
    sm <- nu*sqrt(sum(((deltam - deltam.old)%*%D)^2))
    
    deltam.old <- deltam
    
    if(rm <= tolpri & sm <= toldual){break}
    
  }
  
  if(m == maxiter) {flag <- 1}
  
  groupest <- getgroup(deltam = deltam,n = nobs)
  ngest <- length(unique(groupest))
  
  if(ncx ==1){
    alpest <-  matrix(by(betam, groupest, colMeans,simplify = TRUE), nrow= 1)
  }else{  alpest <- do.call("cbind",by(betam, groupest, colMeans,simplify = TRUE))}
  

  betaest <- t(alpest[,groupest,drop= FALSE])
  
  etaest <- ztz %*% (yhat - rowSums(x*betaest))
  
  
  # marginal likelihood 
  loglikvalue <- -0.5*sum(log(deltae + sig2)) - 0.5* sum((yhat - z%*%etaest-rowSums(x*betaest))^2/(deltae + sig2))
  
  outls <- list(betaest = betaest, etaest = etaest, betam = betam, sig2 = sig2, 
                group = groupest, mhat = mhat,
                deltam = deltam, flag = flag,loglikvalue = loglikvalue,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                niteration = m)
  return(outls)
}