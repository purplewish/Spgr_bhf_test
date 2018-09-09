library(numDeriv)
nv <- rep(c(20,30,40),c(20,20,20))
sigmae2 <- 1
sigmav2 <- 1

groupv <- rep(1:3, each = 20)
betag <- matrix(c(1,1,3,3,5,5),ncol=2,byrow = TRUE)[groupv,]



xm <- cbind(1,rnorm(sum(nv)))
#betamv <- c(1,1)
yv <- rep(rnorm(length(nv))*sqrt(sigmav2),nv) + rnorm(sum(nv))*sqrt(sigmae2)

for(i in 1:length(nv))
{
  if(i==1){index1=1}else{
    index1 <- sum(nv[1:(i-1)]) + 1
  }
  index2 <- sum(nv[1:i])
  
  yv[index1:index2] <- yv[index1:index2] +  xm[index1:index2,]%*%betag[i,]
}

yv = dats$y
xm = cbind(1,dats$x)
indexy = dats$dom
betag = betamdom

uindexy = unique(indexy)

loglik <- function(par)
{
  sigmav2 <- par[1]
  sigmae2 <- par[2]

  value <- 0
  for(i in 1:length(nv))
  {
    n1 <- nv[i]
    mat1 <- matrix(1,ncol=1,nrow = n1)
    diag1 <- diag(1, n1,n1)
    
    S1 <- mat1%*%t(mat1)*sigmav2 + diag1*sigmae2
    S1inv <- solve(S1)
    
    # if(i==1){index1=1}else{
    #   index1 <- sum(nv[1:(i-1)]) + 1
    # }
    # index2 <- sum(nv[1:i])
    # 
    indexi = indexy == uindexy[i]
    #value <- value  -0.5*sum(log(eigen(S1)$values)) - 0.5*t(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])%*%S1inv%*%(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])
    
    value <- value  -0.5*sum(log(eigen(S1)$values)) - 0.5*t(yv[indexi] - xm[indexi,]%*%betag[i,])%*%S1inv%*%(yv[indexi] - xm[indexi,]%*%betag[i,])
  }
  
  return(value)
}



d1 <- function(par)
{
  sigmav2 <- par[1]
  sigmae2 <- par[2]
  
  d11 <- d12 <- 0
  for(i in 1:length(nv))
  {
    n1 <- nv[i]
    mat1 <- matrix(1,ncol=1,nrow = n1)
    diag1 <- diag(1, n1,n1)
    
    M1 <- - mat1%*%t(mat1)*(1/(sigmae2 + n1*sigmav2)^2)
    M2 <- (sigmav2*(2*sigmae2 + n1*sigmav2)/(sigmae2+ n1*sigmav2)^2* mat1%*%t(mat1) - diag1)/(sigmae2)^2
    
    if(i==1){index1=1}else{
      index1 <- sum(nv[1:(i-1)]) + 1
    }
    index2 <- sum(nv[1:i])
    
    d11 <- d11  -0.5*n1/(sigmae2+ n1*sigmav2) - 0.5*t(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])%*%M1%*%(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])
    d12 <- d12 - 0.5*n1*(1- sigmav2/(sigmae2 +n1*sigmav2))/sigmae2 - 0.5*t(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])%*%M2%*%(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])
  }
  return(c(d11,d12))
}


optim(par = c(1,1),fn = loglik,method = "BFGS",control = list(fnscale=-1))

I1 <- function(par)
{
  sigmav2 <- par[1]
  sigmae2 <- par[2]
  
  I11 <- 0.5*sum(nv^2/(sigmae2 + nv*sigmav2)^2)
  I12 <- 0.5*sum(nv/(sigmae2 + nv*sigmav2)^2)
  I22 <- 0.5*sum(nv*(1- sigmav2*(2*sigmae2 + nv*sigmav2)/(sigmae2 + nv*sigmav2)^2)/sigmae2^2)
  Imat1 <- matrix(0,2,2)
  Imat1[1,1] <- I11
  Imat1[1,2] <- Imat1[2,1] <- I12
  Imat1[2,2] <- I22
  return(Imat1)
}


est <- c(1,1)
for(j in 1:10)
{
  est <- est + solve(I1(est))%*%d1(est)
}


d1(c(0.22,0.51))
grad(loglik, c(0.22,0.51))


loglik2 <- function(par)
{
  sigmav2 <- par[1]
  sigmae2 <- par[2]
  
  value <- 0
  for(i in 1:length(nv))
  {
    n1 <- nv[i]
    mat1 <- matrix(1,ncol=1,nrow = n1)
    diag1 <- diag(1, n1,n1)
    
    S1 <- mat1%*%t(mat1)*sigmav2 + diag1*sigmae2
    S1inv <- solve(S1)
    
    # if(i==1){index1=1}else{
    #   index1 <- sum(nv[1:(i-1)]) + 1
    # }
    # index2 <- sum(nv[1:i])
    # 
    indexi = indexy == uindexy[i]
    #value <- value  -0.5*sum(log(eigen(S1)$values)) - 0.5*t(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])%*%S1inv%*%(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])
    
    value <- value  -0.5*sum(log(eigen(S1)$values))/n1 - 0.5*t(yv[indexi] - xm[indexi,]%*%betag[i,])%*%S1inv%*%(yv[indexi] - xm[indexi,]%*%betag[i,])/n1
  }
  
  return(value)
}

optim(par = c(1,1),fn = loglik2,method = "BFGS",control = list(fnscale=-1))

d2 <- function(par)
{
  sigmav2 <- par[1]
  sigmae2 <- par[2]
  
  d11 <- d12 <- 0
  for(i in 1:length(nv))
  {
    n1 <- nv[i]
    mat1 <- matrix(1,ncol=1,nrow = n1)
    diag1 <- diag(1, n1,n1)
    
    M1 <- - mat1%*%t(mat1)*(1/(sigmae2 + n1*sigmav2)^2)
    M2 <- (sigmav2*(2*sigmae2 + n1*sigmav2)/(sigmae2+ n1*sigmav2)^2* mat1%*%t(mat1) - diag1)/(sigmae2)^2
    
    if(i==1){index1=1}else{
      index1 <- sum(nv[1:(i-1)]) + 1
    }
    index2 <- sum(nv[1:i])
    
    d11 <- d11  -0.5/(sigmae2+ n1*sigmav2) - 0.5*t(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])%*%M1%*%(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])/n1
    d12 <- d12 - 0.5*(1- sigmav2/(sigmae2 +n1*sigmav2))/sigmae2 - 0.5*t(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])%*%M2%*%(yv[index1:index2] - xm[index1:index2,]%*%betag[i,])/n1
  }
  return(c(d11,d12))
}


I2 <- function(par)
{
  sigmav2 <- par[1]
  sigmae2 <- par[2]
  
  I11 <- 0.5*sum(nv/(sigmae2 + nv*sigmav2)^2)
  I12 <- 0.5*sum(1/(sigmae2 + nv*sigmav2)^2)
  I22 <- 0.5*sum((1- sigmav2*(2*sigmae2 + nv*sigmav2)/(sigmae2 + nv*sigmav2)^2)/sigmae2^2)
  Imat1 <- matrix(0,2,2)
  Imat1[1,1] <- I11
  Imat1[1,2] <- Imat1[2,1] <- I12
  Imat1[2,2] <- I22
  return(Imat1)
}


est <- c(1,1)
for(j in 1:10)
{
  est <- est + solve(I2(est))%*%d2(est)
}

S1 <- mat1%*%t(mat1)*sigmav2 + diag1*sigmae2
S1inv <- solve(S1)

sum(diag(S1inv  %*% S1inv))

n1^2/(sigmae2 + n1*sigmav2)^2

sum(diag(S1inv %*% S1inv))

n1*(1- sigmav2*(2*sigmae2 + n1*sigmav2)/(sigmae2 + n1*sigmav2)^2)/sigmae2^2

library(Rcpp)
sourceCpp("code/Spgr_bhf.cpp")
sourceCpp("code/Spgr_bhf3.cpp")
t1 <- Sys.time()
res3 <- Spgr_bhf2(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam01,maxiter = 2000,lam = 0.5)
t2 <- Sys.time()


res31 <- Spgr_bhf(indexy = dats$dom, y = dats$y, x = cbind(1, dats$x),cvec = c1,betam0 = betam01,maxiter = 2000,lam = 0.5)
t3 <- Sys.time()
res4 <- dfun3(par = c(1,1),betam = betag,n = 99,nJ = nv,matls = matls,diagls = diagls,indexy = indexy,uindexy = uindexy,y = y,x = x)

res5 <- Ifun3(par = sig2est,nJ = nv)


res32 <- Spgr_bhf2(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam,maxiter = 1,lam = 0.8)


res33 <- Spgr_bhfw(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),cvec = c1,betam0 = betam,maxiter = 1,lam = 0.5)

res34 <- Spgr_bhf3(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam,maxiter = 1000,lam = 0.75)
res35 <- Spgr_bhf2(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = c1,betam0 = betam,maxiter = 1000,lam = 0.8)

cl = exp(1*(1-ordermat))
res36 <- Spgr_bhf3(indexy = dats$dom,y = dats$y,x = cbind(1,dats$x),weights = cl,betam0 = betam,maxiter = 1000,lam = 1)

betag = res32$betaest[res32$group,]
loglik2(res34$sig2est)

loglik3(par = as.numeric(res33$sig2est),indexy = dats$dom,y = dats$y,x = cbind(1, dats$x),betam = res33$betaest,matls = matls,diagls = diagls,nJ = nv)


