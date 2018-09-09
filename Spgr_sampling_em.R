
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


Spgr_sampling_em <- function(yhat, x, sde, cvec, betam0,
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
    mhat <- (sig2*yhat + deltae*muest)/(sig2 + deltae)
    Vhat <- deltae*sig2/(sig2 + deltae)
    
    XtX <- t(Xm)%*%Xm
    Xinv <- solve(XtX + nu*AtA)
    beta.new <-  Xinv %*% (t(Xm)%*%mhat + nu * c((deltam.old -  vm/nu) %*% D))
    
    betam <- matrix(beta.new, nobs, ncx, byrow = TRUE)
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    
    groupest <- getgroup(deltam = deltam,n = nobs,tol = 1e-2)
    ngest <- length(unique(groupest))
    
    sig2 <- mean((mhat - rowSums(x * betam))^2 + Vhat)
    
    
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
  
  muest <- rowSums(x * betam)
  mhat <- (sig2*yhat + deltae*muest)/(sig2 + deltae)
  Vhat <- deltae*sig2/(sig2 + deltae)
  
  sig2 <- mean((mhat - rowSums(x * betam))^2 + Vhat)
  
  # marginal likelihood 
  loglikvalue <- -0.5*sum(log(deltae + sig2)) - 0.5* sum((yhat - rowSums(x*betaest))^2/(deltae + sig2))
  
  outls <- list(betaest = betaest, betam = betam, sig2 = sig2, group = groupest, mhat = mhat, Vhat = Vhat,
                deltam = deltam, flag = flag,loglikvalue = loglikvalue,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                niteration = m)
  return(outls)
}


###### group information is known #####
Spgr_sampling_em_fixed <- function(yhat, x, sde, groupest, betav0, maxiter=1000, tol= 1e-4)
{
  nobs <- length(yhat)
  ncx <- ncol(x)
  Ip <- diag(1,ncx,ncx)
  
  ngest <- length(unique(groupest))

  W1 <- matrix(0,nobs, ngest*ncx)
  for(j in 1:ngest)
  {
    W1[groupest == j, (ncx*(j-1)+1):(ncx*j)] <- x[groupest == j]
  }

  
  #### initial values 
  deltae <- sde^2
  betav <- betav0
  ###########
  sig2 <- mean((yhat - W1%*%betav)^2)
  
  flag <- 0
  for(m in 1:maxiter)
  {
    muest <- W1%*%betav
    mhat <- (sig2*yhat + deltae*muest)/(sig2 + deltae)
    Vhat <- deltae*sig2/(sig2 + deltae)
    
    betavnew <- solve(t(W1)%*%W1)%*%t(W1)%*%mhat
    
    sig2 <- mean((mhat - W1%*%betavnew)^2 + Vhat)
    
    diffnorm <- sum((betav - betavnew)^2)
    
    betav <- betavnew
    
    if(diffnorm <= tol){break}
    
  }
  
  if(m == maxiter) {flag <- 1}
  
  muest <- W1%*%betav
  mhat <- (sig2*yhat + deltae*muest)/(sig2 + deltae)
  Vhat <- deltae*sig2/(sig2 + deltae)
  
  betaest <- matrix(betav, nrow = ncx)
  
  Omega <- diag(1/(sig2+ deltae),nobs,nobs)
  
  covm <- solve(t(W1)%*%Omega%*%W1)

  outls <- list(betaest = betaest,sig2 = sig2, mhat = mhat, Vhat = Vhat,
                flag = flag,niteration = m, covm = covm)
  return(outls)
}