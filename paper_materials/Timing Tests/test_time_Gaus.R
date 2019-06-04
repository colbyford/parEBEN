setwd("~/Downloads/par_EBEN/")
# library(devtools)
# install_github("colbyford/parEBEN")
# library(parEBEN)
rm(list=ls())
library("EBEN")
library(foreach)
library(iterators)
library(parallel)
library(parEBEN)
library(readr)

genotype <- read.table("./genotype_full.txt",header = T)
pheno <- read.table("./pheno_left.txt",header=F)
geno1 <- t(genotype)
pheno <- as.matrix(pheno)
pheno1 <- matrix(pheno[,2],ncol=1)
geno2 <- geno1[,which(geno1[1,] %in% pheno[,1])]
geno3 <- t(geno2[2:nrow(geno2),])

tests <- read.csv("./test_time_template_gaussian.csv")
## Register your local cluster.
library(doParallel)
#no_cores <- detectCores()
no_cores <- 8
cl <- makeCluster(no_cores)
#clusterExport(cl, c("parEBEN.cv.doParallel"))
registerDoParallel(cl)
#stopImplicitCluster()

for (i in 1:nrow(tests)){
#for (i in c(60)){
  n <- tests[i,1]
  k <- tests[i,2]
  nFolds <- tests[i,3]
  set.seed(1)
  sample_n <- sample(1:nrow(geno3),n,replace = F)
  sample_k <- sample(1:ncol(geno3),k,replace = F)
  BASISset <- geno3[sample_n,sample_k]
  class(BASISset) <- "numeric"
  yset  <- matrix(pheno1[1:n,],ncol=1)
  class(yset) <- "numeric"
  cat("Iteration:",i ," (n:", n,"k:", k,"nFolds:", nFolds,")\n")
  sertime <- system.time(EBelasticNet.GaussianCV(BASISset,yset, nFolds, Epis = "yes"))
  tests[i,4] <- sertime[1]
  tests[i,5] <- sertime[2]
  tests[i,6] <- sertime[3]
  partime <- system.time(CrossValidate(BASISset, yset, nFolds, Epis = "yes", prior = "gaussian", search = "local"))
  tests[i,7] <- partime[1]
  tests[i,8] <- partime[2]
  tests[i,9] <- partime[3]
  gc()
  write_csv(tests[i,], "testoutput_time_4-12-2018_gaussian_EPIS_cf.csv", append = TRUE)
}



