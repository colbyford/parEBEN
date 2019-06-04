#rm(list=ls())
library("EBEN")
library(parEBEN)
library(readr)
library(dplyr)

phenotype <- read_delim("pheno1", "\t", col_names = c("y"))
genotype <- read_delim("filter_matrix_looser_0.02_main_0.15_epi-19872","\t", col_names = c("SAMID",paste0("x",seq(1:19871))))


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
#library(Rmpi)

# create and register a doMPI cluster if necessary
if (!identical(getDoParName(), 'doMPI')) {
  # set count to (cores_requested-1)
  cl <- startMPIcluster(count=215,verbose=TRUE)
# clusterExport(cl, c("CrossValidate","AssignToFolds","BuildGrid","GetModelError","LocalSearch","TestModel"))
  registerDoMPI(cl)
}

#library(doSNOW)
#cl <- makeMPIcluster(215)
#registerDoSNOW(cl)

output <- CrossValidate(as.matrix(genotype[,2:ncol(genotype)]),
                        phenotype$y,
                        nFolds = 3,
                        Epis = "no",
                        prior = "gaussian",
                        search = "global")

saveRDS(output,"FullTest_19872_ParCV_8-7-2018.RDS")

closeCluster(cl)
mpi.quit()
