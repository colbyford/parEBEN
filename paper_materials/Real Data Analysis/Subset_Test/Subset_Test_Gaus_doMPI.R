rm(list=ls())
library("EBEN")
library(parEBEN)
library(readr)
library(dplyr)

phenotype <- read_delim("pheno1", "\t", col_names = c("y"))
genotype <- read_delim("filter_matrix","\t", col_names = c("SAMID",paste0("x",seq(1:5356))))


## Register your local cluster.
#library(doParallel)
#no_cores <- detectCores()
#cl <- makeCluster(no_cores-1, outfile="Log.txt")
#cl <- makeCluster(no_cores)
#clusterExport(cl, c("CrossValidate"))
#registerDoParallel(cl)
#stopImplicitCluster()
#paste("Using", no_cores, "cores.")

library(doMPI)

# create and register a doMPI cluster if necessary
if (!identical(getDoParName(), 'doMPI')) {
  # set count to (cores_requested-1)
  cl <- startMPIcluster(count=255,verbose=TRUE)
# clusterExport(cl, c("CrossValidate","AssignToFolds","BuildGrid","GetModelError","LocalSearch","TestModel"))
  registerDoMPI(cl)
}


output <- CrossValidate(as.matrix(genotype[,2:ncol(genotype)]),
                        phenotype$y,
                        nFolds = 3,
                        Epis = "no",
                        prior = "gaussian",
                        search = "global")

saveRDS(output,"SubsetParCV_5-2-2018.RDS")

closeCluster(cl)
