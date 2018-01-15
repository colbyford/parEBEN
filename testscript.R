#Test Script for parEBEN

library(EBEN)
data(BASISbinomial)
data(yBinomial)
#reduce sample size to speed up the running time
n <- 50;
k <- 200;
N <- length(yBinomial);
set.seed(1)
set  <- sample(N,n);
BASIS <- BASISbinomial[set,1:k];
y <- yBinomial[set];
nFolds <- 3

output <- EBelasticNet.Binomial(BASIS, y,lambda = serCV$Lambda_optimal,alpha = serCV$Alpha_optimal, Epis = "no",verbose = 5)

library(doParallel)  
no_cores <- detectCores() 
cl <- makeCluster(no_cores)
clusterExport(cl, c("parEBEN.cv.doParallel"))
registerDoParallel(cl) 
#stopCluster(cl)

#Other Variables for Testing
Target <- y
Epis <- "no"
foldId <- 0
i <- 1

serCV <- EBelasticNet.BinomialCV(BASIS, y, nFolds = 3,Epis = "no")
parCV <- parEBEN.cv.doParallel(BASIS, y, nFolds = 3,Epis = "no", prior = "binomial")


system.time(serCV <- EBelasticNet.BinomialCV(BASIS, y, nFolds = 3,Epis = "no"), gcFirst = TRUE)
system.time(parCV <- parEBEN.cv.doParallel(BASIS, y, nFolds = 3,Epis = "no", prior = "binomial"), gcFirst = TRUE)


serCV <- EBelasticNet.GaussianCV(BASIS, y, nFolds = 3,Epis = "no")
parCV <- parEBEN.cv.doParallel(BASIS, y, nFolds = 3,Epis = "no", prior = "gaussian")
