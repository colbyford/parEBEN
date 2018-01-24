library(EBEN)
library(parEBEN)
data(BASISbinomial)
data(yBinomial)

## Register your local cluster.
library(doParallel)  
no_cores <- detectCores() 
cl <- makeCluster(no_cores)
#clusterExport(cl, c("parEBEN.cv.doParallel"))
registerDoParallel(cl)

## Generate your test matrix
tests <- data.frame( n = c(25,25,50,50),
                     k = c(100,200,100,200),
                     nFolds = c(3,3,3,3),
                     ser_user_time = c(0,0,0,0),
                     ser_system_time = c(0,0,0,0),
                     ser_elapsed_time = c(0,0,0,0),
                     par_user_time = c(0,0,0,0),
                     par_system_time = c(0,0,0,0),
                     par_elapsed_time = c(0,0,0,0)
                     )

## Loop throught the tests in the matrix and report the time
for (i in 1:nrow(tests)){
  n <- tests[i,1]
  k <- tests[i,2]
  nFolds <- tests[i,3]
  N <- length(yBinomial)
  set.seed(1)
  set <- sample(N,n)
  BASIS <- BASISbinomial[set,1:k]
  y  <- yBinomial[set]
  sertime <- system.time(EBelasticNet.BinomialCV(BASIS,y, nFolds, Epis = "no"))
  tests[i,4] <- sertime[1]
  tests[i,5] <- sertime[2]
  tests[i,6] <- sertime[3]
  partime <- system.time(parEBEN.cv(BASIS, y, nFolds, Epis = "no", parMethod = "doParallel", prior = "binomial"))
  tests[i,7] <- partime[1]
  tests[i,8] <- partime[2]
  tests[i,9] <- partime[3]
}
