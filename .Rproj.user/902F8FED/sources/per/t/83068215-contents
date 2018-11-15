## To use doParallel or doMPI funtionality on a multicore/multi-machine environment, you must register your local cluster.
library(doParallel)  #or library(doMPI)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#clusterExport(cl, c("CrossValidate"))
registerDoParallel(cl)

## Load in data and required EBEN and parEBEN packages
library(EBEN)
library(parEBEN)

## Create small sample matrix for testing
data(BASIS)
data(y)
n = 50
k = 100
BASIS = BASIS[1:n,1:k]
y  = y[1:n]

parEBENcv <- CrossValidate(BASIS, y, nFolds = 3, Epis = "no", prior = "gaussian", search = "global")

#Use the optimal values in the EBEN model
EBENoutput <- EBelasticNet.Gaussian(BASIS, y, lambda = parEBENcv$lambda.optimal, alpha = parEBENcv$alpha.optimal, Epis = "no", verbose = 1)
