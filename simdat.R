##### simulate data #####
# N: number of domains
# Nv: set of population size
# srh: higher sampling rate
# srl: lower sampling rate
# Ax: area of domain 
# ax bx beta distribution parameters
# meanx: mean of x
# sigx: sd of x
# sige: sd of error
# sigv: sd of domain random effect
# sigb: sd of beta 
# model1: linear regression with common parameters


simdat <- function(N = 99, Nv, srh = 0.3, srl = 0.1, Ax, ax = 1, bx = 10,
                   model = 1, beta = NULL, betam = NULL, 
                   sde = 4, sdv = NULL, sdb0 = NULL, sdb1 = NULL, seed = 1019)
{
  set.seed(seed)
  
  # group1 has higher sampling rate and group 2 has lower sampling rate
  group <- sample(1:2,size = N,replace = TRUE)
  
  Nv0 <- sample(Nv, size = N, replace = TRUE) ## population domain size 
  
  Nt <- sum(Nv0) 
  
  x <- rbeta(Nt, ax, bx)* rep(Ax,  Nv0)
  
  if(model == 1)
  {
    y <- beta[1] + beta[2]*x + rnorm(Nt)*sde
  }
  
  if(model == 2)
  {
    domeff <- rep(rnorm(N)*sdv, Nv0)
    y = beta[1] + beta[2]*x + domeff + rnorm(Nt)*sde
  }
  
  if(model == 3)
  {
    b0 <- beta[1] + rnorm(N)*sdb0
    b1 <- beta[2] + rnorm(N)*sdb1
    b0 <- rep(b0, Nv0)
    b1 <- rep(b1, Nv0)
    y <- b0 + b1*x + rnorm(Nt)*sde
  }
  
  ### sample procedure ### 
  # srs 
  
  nv <- rep(0, N) # sample size in each domain
  nv[group==1] <- round(Nv0[group==1]*srh)
  nv[group==2] <- round(Nv0[group==2]*srl)
  n0 <- sum(nv)
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x", "weights")
  dats$dom <- rep(1:N, nv)
  
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(1:N, Nv0)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  datp$y <- y
  datp$x <- x
  
  for(i in 1:N)
  {
    indpi <- indp == i
    ypi <- y[indpi]
    xpi <- x[indpi]
    
    indi <- dats$dom == i
   
    inds <- sample.int(Nv0[i], nv[i])
    
    dats$y[indi] <- ypi[inds]
    dats$x[indi] <- xpi[inds]
    dats$weights[indi] <- Nv0[i]/nv[i]
  }
  
  domsize <- data.frame(dom = 1:N, size = Nv0, sampsize = nv)
  out <- list(sample = dats, pop = datp, domsize = domsize)
  return(out)
}


simpop <- function(N = 99, Nv, Ax, ax = 1, bx = 10,
                   model = 1, beta = NULL, 
                   sde = 4, sdv = NULL, sdb0 = NULL, sdb1 = NULL, seed = 1019)
{
  set.seed(seed)
  
  # group1 has higher sampling rate and group 2 has lower sampling rate
  group <- sample(1:2,size = N,replace = TRUE)
  
  Nv0 <- sample(Nv, size = N, replace = TRUE) ## population domain size 
  
  Nt <- sum(Nv0) 
  
  x <- rbeta(Nt, ax, bx)* rep(Ax,  Nv0)
  
  if(model == 1)
  {
    y <- beta[1] + beta[2]*x + rnorm(Nt)*sde
  }
  
  if(model == 2)
  {
    domeff <- rep(rnorm(N)*sdv, Nv0)
    y = beta[1] + beta[2]*x + domeff + rnorm(Nt)*sde
  }
  
  if(model == 3)
  {
    b0 <- beta[1] + rnorm(N)*sdb0
    b1 <- beta[2] + rnorm(N)*sdb1
    b0 <- rep(b0, Nv0)
    b1 <- rep(b1, Nv0)
    y <- b0 + b1*x + rnorm(Nt)*sde
  }
  
  datp <- as.data.frame(matrix(0, Nt, 3))
  indp <- rep(1:N, Nv0)
  colnames(datp) <- c("dom","y","x")
  datp$dom <- indp
  datp$y <- y
  datp$x <- x
  
  domsize <- data.frame(dom = 1:N, size = Nv0)
  out <- list(pop = datp, domsize = domsize, group = group)
  return(out)
}


simsrs <- function(N, group, Nv0, datp, srh = 0.3, srl = 0.1, seed = 2233)
{
  set.seed(seed)
  nv <- rep(0, N) # sample size in each domain
  nv[group==1] <- round(Nv0[group==1]*srh)
  nv[group==2] <- round(Nv0[group==2]*srl)
  n0 <- sum(nv)
  
  dats <- as.data.frame(matrix(0, n0, 4))
  colnames(dats) <- c("dom", "y", "x", "weights")
  dats$dom <- rep(1:N, nv)
  
  
  y <- datp$y
  x <- datp$x
  indp <- datp$dom
  
  for(i in 1:N)
  {
    indpi <- indp == i
    ypi <- y[indpi]
    xpi <- x[indpi]
    
    indi <- dats$dom == i
    
    inds <- sample.int(Nv0[i], nv[i])
    
    dats$y[indi] <- ypi[inds]
    dats$x[indi] <- xpi[inds]
    dats$weights[indi] <- Nv0[i]/nv[i]
  }
  
  samplesize <- data.frame(dom = 1:N, sampsize = nv)
  out <- list(sample = dats, samplesize = samplesize)
  return(out)
}



### simulation based on model 

sim_model <- function(M, N, Nv,srh, srl, Ax, ax, bx, model, beta, 
                      sde = NULL, sdv = NULL, sdb0 = NULL, sdb1 = NULL,
                      seed0 = 2330)
{
  
  rmse_mat1 <- matrix(0, M, 6)
  colnames(rmse_mat1) <- c("direct","eblup","reg","ridge","lasso","Alasso")
  numcoef1 <- rep(0,M)
  
  for(m in 1: M)
  {
    sim1 <- simdat(N = N,Nv = Nv,srh = srh, srl = srl,
                   Ax = Ax,ax = ax,bx = bx,model = model,beta = beta,
                   sde = sde, sdv = sdv, sdb0 = sdb0, sdb1 = sdb1, seed = seed0 + m)
    
    dats <- sim1$sample
    datp <- sim1$pop
    popn <- sim1$domsize
    meany <- ddply(datp,.(dom),summarize, meany = mean(y))
    
    set.seed(2*m + seed0)
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
    
    estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx,popnsize = popn, 
                       method = "REML",data = dats)
    
    xmp <- matrix(0, nrow(meanx), 2*N)
    xmp[,1] <- 1
    xmp[,2] <- meanx$meanx
    xmp[1,seq(3,2*N,by = 2)] <- -1
    xmp[1,seq(4,2*N,by = 2)] <- - meanx$meanx[1]
    
    xm <- matrix(0, nrow(dats), 2*(N))
    xm[,1] <- 1
    xm[,2] <- dats$x
    xm[1:popn$sampsize[1],seq(3,ncol(xm),by = 2)] <- -1
    xm[1:popn$sampsize[1],seq(4,ncol(xm),by = 2)] <- -dats$x[1:popn$sampsize[1]]
    
    for(i in 2:N)
    {
      index1 <- sum(popn$sampsize[1:(i-1)]) + 1
      index2 <- sum(popn$sampsize[1:i])
      xm[index1:index2,(i-1)*2 +1] <- 1
      xm[index1:index2,2*i] <- dats$x[dats$dom==i]
      
      xmp[i, (i-1)*2 +1]  <- 1
      xmp[i,2*i]  <- meanx$meanx[i]
    }
    
    
    ## regression 
    
    res0 <- lm(y~x, weights = dats$weights, data = dats)
    betaest0 <- coef(res0)
    estreg <- cbind(1,meanx$meanx) %*% betaest0
    
    
    ## ridge regression 
    cv1 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights, alpha = 0, standardize = TRUE)
    betaest1 <- coef(cv1,s = cv1$lambda.min)
    estridge <- as.matrix(xmp %*% betaest1)[,1]
    
    
    ## lasso
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    ## adaptive lasso 
    ww <- 1/abs(betaest1[-1])
    cv3 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1,
                     penalty.factor = ww)
    betaest3 <- coef(cv3,s = cv3$lambda.min)
    estalasso <- as.matrix(xmp%*% betaest3)[,1]
    
    
    rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg)^2),
                   mean((meany[,2] - estridge)^2),
                   mean((meany[,2] - estlasso)^2),
                   mean((meany[,2] - estalasso)^2)))
    rmse_mat1[m,] <- rmse
    
  }
  
  out <- list(rmse = rmse_mat1, numcoef = numcoef1)
  return(out)
}


sim_design <- function(M, N, Nv, srh, srl, Ax, ax, bx, model, beta, 
                       sde = NULL, sdv = NULL, sdb0 = NULL, sdb1 = NULL,
                       seed0 = 2330)
{
  sim1 <- simpop(N = N,Nv = Nv,Ax = Ax,ax = ax,bx = bx,model = model,beta = beta,
                 sde = sde, sdv = sdv, sdb0 = sdb0, sdb1 = sdb1, seed = seed0)
  datp <- sim1$pop
  popn <- sim1$domsize
  meanx <- ddply(datp,.(dom),summarize, meanx = mean(x))
  meany <- ddply(datp,.(dom),summarize, meany = mean(y))
  
  rmse_mat_d1 <- matrix(0, M, 6)
  colnames(rmse_mat_d1) <- c("direct","eblup","reg","ridge","lasso","Alasso")
  numcoef_d1 <- rep(0,M)
  
  for(m in 1: M)
  {
    
    simm <- simsrs(N = N,group = sim1$group, Nv0 = sim1$domsize[,2],datp = sim1$pop,
                   srh = srh,srl = srl,seed = seed0 + m)
    
    dats <- simm$sample
    ssizem <- simm$samplesize[,2]
    
    
    set.seed(2*m + seed0)
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx,popnsize = popn, 
                       method = "REML",data = dats)
    
    xmp <- matrix(0, nrow(meanx), 2*N)
    xmp[,1] <- 1
    xmp[,2] <- meanx$meanx
    xmp[1,seq(3,2*N,by = 2)] <- -1
    xmp[1,seq(4,2*N,by = 2)] <- - meanx$meanx[1]
    
    xm <- matrix(0, nrow(dats), 2*(N))
    xm[,1] <- 1
    xm[,2] <- dats$x
    xm[1:ssizem[1],seq(3,ncol(xm),by = 2)] <- -1
    xm[1:ssizem[1],seq(4,ncol(xm),by = 2)] <- -dats$x[1:ssizem[1]]
    
    for(i in 2:N)
    {
      index1 <- sum(ssizem[1:(i-1)]) + 1
      index2 <- sum(ssizem[1:i])
      xm[index1:index2,(i-1)*2 +1] <- 1
      xm[index1:index2,2*i] <- dats$x[dats$dom==i]
      
      xmp[i, (i-1)*2 +1]  <- 1
      xmp[i,2*i]  <- meanx$meanx[i]
    }
    
    
    ## regression 
    
    res0 <- lm(y~x, weights = dats$weights, data = dats)
    betaest0 <- coef(res0)
    estreg <- cbind(1,meanx$meanx) %*% betaest0
    
    
    ## ridge regression 
    cv1 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights, alpha = 0, standardize = TRUE)
    betaest1 <- coef(cv1,s = cv1$lambda.min)
    estridge <- as.matrix(xmp %*% betaest1)[,1]
    
    
    ## lasso
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef_d1[m] <- sum(betaest2!=0)
    
    ## adaptive lasso 
    ww <- 1/abs(betaest1[-1])
    cv3 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1,
                     penalty.factor = ww)
    betaest3 <- coef(cv3,s = cv3$lambda.min)
    estalasso <- as.matrix(xmp%*% betaest3)[,1]
    
    
    rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg)^2),
                   mean((meany[,2] - estridge)^2),
                   mean((meany[,2] - estlasso)^2),
                   mean((meany[,2] - estalasso)^2)))
    rmse_mat_d1[m,] <- rmse
    
  }
  
  out <- list(rmse = rmse_mat_d1, numcoef = numcoef_d1)
  return(out)
}
