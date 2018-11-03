### parallel ###
library(doParallel) 
library(foreach)
source("sim_bhf3_parallel.R")

load("Cmat_popn.RData")
ordermat <- getorder(Matrix(Cmatia))

betad05 <-  matrix(c(0.5,0.5,1,1,1.5,1.5),ncol=2,byrow = TRUE)

subfun1 <- function(mm)
{
  sim_bhf3_parallel(beta = betad05,sdx = 1,mux = 1,sdv = 1,sde = 2,rate = 0.01,popn = popn[,1:2],maxiter =2000, group = popn[,3],seed = mm + 82245)
}

subfun2 <- function(mm)
{
  sim_bhf3_parallel(beta = betad05,sdx = 1,mux = 1,sdv = 1,sde = 1,rate = 0.01,popn = popn[,1:2],maxiter =2000, group = popn[,3],seed = mm + 82245)
}

subfun3 <- function(mm)
{
  sim_bhf3_parallel(beta = betad05,sdx = 1,mux = 1,sdv = 1,sde = 0.5,rate = 0.01,popn = popn[,1:2],maxiter =2000, group = popn[,3],seed = mm + 82245)
}

subfun4 <- function(mm)
{
  sim_bhf3_parallel(beta = betad05,sdx = 1,mux = 1,sdv = 0.5,sde = 0.5,rate = 0.01,popn = popn[,1:2],maxiter =2000, group = popn[,3],seed = mm + 82245)
}

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
result_beta05_x1v1e2 <- foreach(mm=1:100,.packages=c("plyr","sae","Spgr","SpgrBHF")) %dopar%{subfun1(mm)}
result_beta05_x1v1e1 <- foreach(mm=1:100,.packages=c("plyr","sae","Spgr","SpgrBHF")) %dopar%{subfun2(mm)}
result_beta05_x1v1e05 <- foreach(mm=1:100,.packages=c("plyr","sae","Spgr","SpgrBHF")) %dopar%{subfun3(mm)}
result_beta05_x1v05e05 <- foreach(mm=1:100,.packages=c("plyr","sae","Spgr","SpgrBHF")) %dopar%{subfun4(mm)}
stopCluster(cl)

save(result_beta05_x1v1e2, file = "result_beta05_x1v1e2.RData")
save(result_beta05_x1v1e1, file = "result_beta05_x1v1e1.RData")
save(result_beta05_x1v1e05, file = "result_beta05_x1v1e05.RData")
save(result_beta05_x1v05e05, file = "result_beta05_x1v05e05.RData")
