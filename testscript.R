#Test Script for parEBEN

library(EBEN)
data(BASISbinomial)
data(yBinomial)
#reduce sample size to speed up the running time
n = 50;
k = 100;
N = length(yBinomial);
set.seed(1)
set  = sample(N,n);
BASIS = BASISbinomial[set,1:k];
y = yBinomial[set];
nFolds = 3
CV = EBelasticNet.BinomialCV(BASIS, y, nFolds = 3,Epis = "yes")
output = EBelasticNet.Binomial(BASIS, y,lambda = CV$Lambda_optimal,alpha = CV$Alpha_optimal, Epis = "yes",verbose = 5)

library(doParallel)  
no_cores <- detectCores() 
cl <- makeCluster(no_cores)
clusterExport(cl, c("parEBEN.cv.doParallel"))
registerDoParallel(cl) 
#stopCluster(cl)

parcv <- parEBEN.cv.doParallel(BASIS, y, nFolds = 3,Epis = "yes", prior = "binomial")
