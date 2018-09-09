library(plyr)
library(sae)
library(glmnet)
dat_nriia <- read.csv("/Volumes/cssm_groups$/Grads/XinWang/NRI/data/Pgen/PgenAllWithLonLatIDIA.csv")
urbants <- read.csv("/Users/Xin/Research/NRI_urban/data/ia_urban_area_800m_ts.csv",stringsAsFactors = FALSE)

### nri data segment level area 
datia = ddply(dat_nriia,.(COUNTY, LONID, LATID),summarize,
              nri97 = sum(WEIGHT*(BU1997==7))/sum(WEIGHT), nri01 = sum(WEIGHT*(BU2001==7))/sum(WEIGHT), 
              nri06 = sum(WEIGHT*(BU2006==7))/sum(WEIGHT),nri11 = sum(WEIGHT*(BU2011==7))/sum(WEIGHT))

datia <- datia[!is.nan(datia$nri97),]

datia$ts11 <- datia$ts06 <- datia$ts01 <- datia$ts97 <- rep(0, nrow(datia))
counties <- unique(datia$COUNTY)

distmin <- rep(0, nrow(datia))

## match with ts data
for(j in counties)
{
  datj <- datia[datia$COUNTY==j,] # county j 
  longj <- range(datj$LONID)
  latj <- range(datj$LATID)
  
  distminj <- rep(0, nrow(datj))
  
  urbantsj <- urbants[urbants$x >= longj[1]- 0.1 & urbants$x <= longj[2] + 0.1 & urbants$y >=latj[1] - 0.1 & urbants$y <= latj[2] + 0.1,]
  
  for(k in 1:nrow(datj))
  {
    distk <- (urbantsj$x - datj$LONID[k])^2 + (urbantsj$y - datj$LATID[k])^2
    indexjk <- which.min(distk)
    
    datj$ts97[k] <- urbantsj[indexjk,"A1997"]
    datj$ts01[k] <- urbantsj[indexjk,"A2001"]
    datj$ts06[k] <- urbantsj[indexjk,"A2006"]
    datj$ts11[k] <- urbantsj[indexjk,"A2011"]
    
    distminj[k] <- sqrt(min(distk))
  }
  
  datj$ts97 <- datj$ts97*247.105
  datj$ts01 <- datj$ts01*247.105
  datj$ts06 <- datj$ts06*247.105
  datj$ts11 <- datj$ts11*247.105
  
  #datj$area <- datj$area*160
  
  
  datia[datia$COUNTY==j,"ts97"] <- datj$ts97
  datia[datia$COUNTY==j,"ts01"] <- datj$ts01
  datia[datia$COUNTY==j,"ts06"] <- datj$ts06
  datia[datia$COUNTY==j,"ts11"] <- datj$ts11
  
  distmin[datia$COUNTY==j] <- distminj
  
}

datia$nri97 <- datia$nri97*160
datia$nri01 <- datia$nri01*160
datia$nri06 <- datia$nri06*160
datia$nri11 <- datia$nri11*160

temp1 <- ddply(datia,.(COUNTY),summarize, nri = sum(nri), ts = sum(ts),nseg = length(PSU))
temp2 <- ddply(datia,.(COUNTY),summarize, nri = mean(nri), ts = mean(ts),nseg = length(PSU))

write.csv(temp2, "/Users/Xin/Research/NRI_urban/data/nseg_urban.csv",row.names = FALSE) # county mean area 
write.csv(datia, "/Users/Xin/Research/NRI_urban/data/datia_urban97.csv",row.names = FALSE) # nseg in NRI

write.csv(datia, "/Users/Xin/Research/NRI_urban/data/datia_urban.csv",row.names = FALSE) 


datia97 <- read.csv("data/datia_urban97.csv",stringsAsFactors = FALSE)
############ relationship at segment level ####
# some figures #

# some examples 


#counties1 <- counties[sample(99,12)]
counties1 <- c(13,39,67,141,113,15,185,19,21,151,171,57)

par(mfrow = c(3,4))
for(i in 1:length(counties1))
{
  plot(datia[datia$COUNTY==counties1[i],c("ts97","nri97")],main = paste("county",counties1[i]))
}

par(mfrow = c(3,4))
for(i in 1:length(counties1))
{
  plot(datia[datia$COUNTY==counties1[i],c("ts11","nri11")],main = paste("county",counties1[i]))
}

tabseg <- rbind(dplyr::filter(datia, COUNTY==39, PSU == "00010PZ"),

dplyr::filter(datia, COUNTY==185, PSU == "00053PZ"),

dplyr::filter(datia, COUNTY==21, PSU == "00005PZ"),

dplyr::filter(datia, COUNTY==141, PSU == "00038PZ"))


dplyr::filter(dat_nriia, COUNTY==39, PSU == "00010PZ")

dplyr::filter(dat_nriia, COUNTY==185, PSU == "00053PZ")

dplyr::filter(dat_nriia, COUNTY==21, PSU == "00005PZ")

dplyr::filter(dat_nriia, COUNTY==141, PSU == "00038PZ")


dat_nriia[dat_nriia$LONID==-93.57128,1:10]
dat_nriia[dat_nriia$LONID==-93.10216,1:10]
dat_nriia[dat_nriia$LONID==-94.92861,1:10]
dat_nriia[dat_nriia$LONID==-95.39317,1:10]

dat_nriia[dat_nriia$LONID==-91.23699,1:10]

head(dplyr::filter(dat_nriia, POINT==99))


tabseg1 <- rbind(dplyr::filter(datia, COUNTY==39, LONID == -93.57128),
                
                dplyr::filter(datia, COUNTY==185, LONID == -93.10216),
                
                dplyr::filter(datia, COUNTY==21, LONID == -94.92861),
                
                dplyr::filter(datia, COUNTY==141, LONID == -95.39317))


dplyr::filter(datia, nri97==160, ts97<1)


dat_nriia[dat_nriia$LONID==-91.50716,1:12]



##### srs in each 
# datp population include x and y
# popn includes number of population size
group <- rep(1, 99)
sr <- 0.2
popn <- temp2[,c("COUNTY","nseg")]
meanx <- temp2[,c("COUNTY","ts")]
meany <- temp2[,c("COUNTY","nri")]
M <- 100
datp <- datia
xname <- "ts"
yname <- "nri"
srs <- function(datp, popn, xname, yname, sr, srh, srl, group, meanx, meany, seed = 2345)
{
  numcoef1 <- rep(0, M)
  rmse_mat <- rmse_mat1 <- matrix(0, M, 6)
  
  colnames(rmse_mat1) <- c("direct","eblup","reg","ridge","lasso","Alasso")
  Nv <- popn[,2]
  
  N <- length(Nv) ## total number of counties
  nv <- round(Nv * sr) ## sample size 
  n0 <- sum(nv) ## total sample size
  
  dats <- as.data.frame(matrix(0, n0, 4)) # sampled data
  colnames(dats) <- c("dom", "y", "x", "weights")
  dats$dom <- rep(popn[,1],nv)
  doms <- popn[,1]
  
  fn <- nv/Nv ## observed proportion
  
  for(m in 1:M)
  {
    set.seed(seed + m)
    for(i in 1:N)
    {
      indpi <- datp[,1] == doms[i]
      ypi <- datp[indpi,yname] ## population y
      xpi <- datp[indpi,xname] ## population x
      
      indi <- dats$dom == doms[i]
      
      inds <- sample.int(Nv[i], nv[i])
      
      dats$y[indi] <- ypi[inds]
      dats$x[indi] <- xpi[inds]
      dats$weights[indi] <- Nv[i]/nv[i]
    }
    
    meansy <- ddply(dats,.(dom), summarize, meansy = mean(y))
    meansx <- ddply(dats,.(dom), summarize, meansx = mean(x))
    meancx <- (Nv* meanx[,2] - nv* meansx[,2])/(Nv - nv)
    
    estdir <- direct(dats$y,dom = dats$dom, sweight = dats$weights,domsize = popn)
    
    estbhf <- eblupBHF(y~x, dom = dom, meanxpop = meanx ,popnsize = popn, 
                       method = "REML",data = dats)
    
    res0 <- lm(y~x, weights = dats$weights, data = dats)
    
    betaest0 <- coef(res0)
    
    estreg <- cbind(1,meanx[,2]) %*% betaest0
    
    estreg1 <- fn*meansy[,2] + (1-fn)*cbind(1,meancx) %*% betaest0
    
    ########### expand matrix ##############
    
    xmp <-  xmc <- matrix(0, nrow(meanx), 2*N)
    xmp[,1] <- 1
    xmp[,2] <- meanx[,2]
    xmp[1,seq(3,2*N,by = 2)] <- -1
    xmp[1,seq(4,2*N,by = 2)] <- - meanx[1,2]
    
    xmc[,1] <- 1
    xmc[,2] <- meancx
    xmc[1,seq(3,2*N,by = 2)] <- -1
    xmc[1,seq(4,2*N,by = 2)] <- - meancx[1]
    
    
    xm <- matrix(0, nrow(dats), 2*(N))
    xm[,1] <- 1
    xm[,2] <- dats$x
    xm[1:nv[1],seq(3,ncol(xm),by = 2)] <- -1
    xm[1:nv[1],seq(4,ncol(xm),by = 2)] <- -dats$x[1:nv[1]]
    
    for(i in 2:N)
    {
      index1 <- sum(nv[1:(i-1)]) + 1
      index2 <- sum(nv[1:i])
      xm[index1:index2,(i-1)*2 +1] <- 1
      xm[index1:index2,2*i] <- dats$x[dats$dom==doms[i]]
      
      xmp[i, (i-1)*2 +1]  <- 1
      xmp[i,2*i]  <- meanx[i,2]
      
      xmc[i, (i-1)*2 +1]  <- 1
      xmc[i,2*i]  <- meancx[i]
    }
    
    cv1 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights, alpha = 0, standardize = TRUE)
    betaest1 <- coef(cv1,s = cv1$lambda.min)
    estridge <- as.matrix(xmp %*% betaest1)[,1]
    
    estridge1 <- fn*meansy[,2] + (1-fn)*estridge
    
    cv2 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1)
    betaest2 <- coef(cv2,s = cv2$lambda.min)
    estlasso <- as.matrix(xmp %*% betaest2)[,1]
    numcoef1[m] <- sum(betaest2!=0)
    
    estlasso1 <- fn*meansy[,2]  +  (1-fn)*as.matrix(xmp %*% betaest2)[,1]
    
    ww <- 1/abs(betaest1[-1])
    cv3 <- cv.glmnet(x = xm[,-1],y = dats$y,weights = dats$weights,alpha = 1,
                     penalty.factor = ww)
    betaest3 <- coef(cv3,s = cv3$lambda.min)
    estalasso <- as.matrix(xmp%*% betaest3)[,1]
    
    estalasso1 <- fn*meansy[,2] + (1-fn)*estalasso
    
    ##### subgroups ###
    load('/Volumes/cssm_groups$/Grads/XinWang/county_adjacency/data/county_adjacency.RData')
    
    Cmatia <- adjMat[['19']]
    rm(adjMat)

    neworder <- as.character(19000+ counties)
    Cmatia <- Cmatia[neworder, neworder]
    order1 <- getorder(Matrix(Cmatia))
    cvec1 <- exp(0.7*(1- order1))


    dats1 <- dats
    dats1$y <- scale(dats1$y)
    dats1$x <- scale(dats1$x)
    betam01 <- matrix(0, length(counties),2)
    for(i in 1:length(counties))
    {
      indexi <- which(dats$dom== counties[i])
      betam01[i, ] <- coef(lm(dats1$y[indexi] ~ dats1$x[indexi]))
    }

    z0 <- matrix(1,ncol = 1, nrow(dats1))
    betam001 <- initialrx(indexy = dats$dom,y = dats1$y,x = as.matrix(cbind(1,dats1$x)),lam0 = 0.0001)
    
    betam002 <- initialr(indexy = dats$dom,y = dats1$y,x = as.matrix(dats1$x),z = z0,lam0 = 0.0001)
    #
    cvec0 <- rep(1, 99*98/2)
    lamv <- seq(0.05,2,0.05)
   #betam1 <- betam01[,2]
    betam1 <- betam001
    bicv0 <- rep(0, length(lamv))
    group0 <- matrix(0, length(counties), length(lamv))
    beta.array0 <- array(0, dim = c(length(counties),2,length(lamv)))
   
    estv <- rep(0, length(lamv))

   # 
   # 
    for(l in 1:length(lamv))
    {
      # res0l<- concavefusionr(indexy = dats1$dom,y = dats1$y,x = as.matrix(dats1$x),z = z0,weights = cvec1,betam0 = as.matrix(betam1),lam = lamv[l])
      # 
      # bicv0[l] <- BICcr(obj = res0l,indexy = dats1$dom,y = dats1$y,x = as.matrix(dats1$x),z = z0,c0 = 0.3)
      # 


      res0l<- concavefusionrx(indexy = dats1$dom,y = dats1$y,x = as.matrix(cbind(1,dats1$x)),weights = cvec1,betam0 = as.matrix(betam1),lam = lamv[l])

      bicv0[l] <- BICcrx(obj = res0l,indexy = dats1$dom,y = dats1$y,x = as.matrix(cbind(1,dats1$x)),c0 = 0.1)



      beta.array0[,,l] <- res0l$beta
      betam1 <- res0l$beta
      group0[,l] <- getgroup(deltam = res0l$deltam,n = 99)
     #estv[l] <- res0l$eta
    }
   #  
    
    plot(bicv0)
    
   indexm <- which.min(bicv0)
   betaest1 <- beta.array0[,,which.min(bicv0)]
   plot(betaest1)
   
  #  groupall <- rep(group0[,indexm],nv)
  #  ng <- length(unique(group0[,indexm]))
  #  
  #  ncx = 1
  #  X0 <- matrix(0, n0,ng*ncx)
  #  for(k in 1:ng)
  #  {
  #    X0[groupall==k,(ncx*(k-1) + 1) : (ncx*k)] <- dats1$x[groupall==k,]
  #  }
  #  
  #  W0 <- cbind(z0, X0)
  # est0 <- solve(t(W0) %*% W0) %*% t(W0) %*% dats1$y
  # eta0 <- est0[1:ncz]
  # alpest0 <- t(matrix(est0[(ncz+1) : length(est0)],ncx))
  # betam <- alpest0[group0,]
   
   # est1 <- (meancx - attr(dats1$x, "scaled:center"))/attr(dats1$x,"scaled:scale") *betaest1 + estv[indexm]
   # esty1 <- est1 * attr(dats1$y,"scaled:scale") + attr(dats1$y, "scaled:center")
   #   
   # estg1 <- fn*meansy[,2] + (1-fn)* esty1
   
   est1 <- rowSums(cbind(1,(meancx - attr(dats1$x, "scaled:center"))/attr(dats1$x,"scaled:scale")) *betaest1)
 
    esty1 <- est1 * attr(dats1$y,"scaled:scale") + attr(dats1$y, "scaled:center")

    estg1 <- fn*meansy[,2] + (1-fn)* esty1
   

    rmse <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg)^2),
                   mean((meany[,2] - estridge)^2),
                   mean((meany[,2] - estlasso)^2),
                   mean((meany[,2] - estalasso)^2)))
    
    rmse1 <- sqrt(c(mean((meany[,2] - estdir$Direct)^2),
                   mean((meany[,2] - estbhf$eblup$eblup)^2),
                   mean((meany[,2] - estreg1)^2),
                   mean((meany[,2] - estridge1)^2),
                   mean((meany[,2] - estlasso1)^2),
                   mean((meany[,2] - estalasso1)^2)))
    
    mean((meany[,2] - estg1)^2)
    
    rmse_mat[m,] <- rmse
    rmse_mat1[m,] <- rmse1
  }
}

save(rmse_mat, rmse_mat1,file = "/Users/Xin/Research/survey/result/urban_nri_ts_02.RData")

boxplot(rmse_mat1)
boxplot(rmse_mat)


boxplot(rmse_mat1)


betam0 <- matrix(0, length(counties),2)
for(i in 1:length(counties))
{
  indexi <- which(datp$COUNTY== counties[i])
  betam0[i, ] <- coef(lm(datp$nri[indexi] ~ datp$ts[indexi]))
}


betam01 <- matrix(0, length(counties),2)
for(i in 1:length(counties))
{
  indexi <- which(dats$dom== counties[i])
  betam01[i, ] <- coef(lm(dats$y[indexi] ~ dats$x[indexi]))
}


residvalue <- dats$y - lm(y~x,data = dats)$fitted.values
plot(density(residvalue[abs(residvalue)<=5]))


residvalue0 <- datp$nri - lm(nri~ts,data = datp)$fitted.values
plot(density(residvalue0[abs(residvalue0)<=5]))


library(concavefusion)





library(shapefiles)
dat1 <- read.dbf("data/IA_TS_800m/IA_UrbanTS_800m.dbf")
dat2 <- read.shx("data/IA_TS_800m/IA_UrbanTS_800m.shx")
dat3 <- read.shp("data/IA_TS_800m/IA_UrbanTS_800m.shp")
dat4 <- shapefile("data/IA_TS_800m/IA_UrbanTS_800m.shp")

urbanall <- read.csv("Research/survey/data/IA_TS_800m/IA_UrbanArea_TS.csv",stringsAsFactors = FALSE)


dat_nriia[dat_nriia$PSU =="020302R" & dat_nriia$COUNTY == 1,c("BU1997","WEIGHT")]

dat_nriia[dat_nriia$LONID==-93.57128,1:10]
dat_nriia[dat_nriia$LONID==-93.10216,1:10]
dat_nriia[dat_nriia$LONID==-94.92861,1:10]
dat_nriia[dat_nriia$LONID==-95.39317,1:10]


