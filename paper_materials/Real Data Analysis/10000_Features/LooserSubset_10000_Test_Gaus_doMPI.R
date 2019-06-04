rm(list=ls())
library("EBEN")
library(parEBEN)
library(readr)
library(dplyr)

print("Loaded Packages")

phenotype <- read_delim("pheno1", "\t", col_names = c("y"))

print("Data loaded: phenotype")

genotype <- read_delim("filter_matrix_looser","\t", col_names = c("SAMID",paste0("x",seq(1:53703))))

print("Data loaded: genotype")

library(doMPI)
# create and register a doMPI cluster if necessary
if (!identical(getDoParName(), 'doMPI')) {
  # set count to (cores_requested-1)
  cl <- startMPIcluster(count=255,verbose=TRUE)
  registerDoMPI(cl)
}

print("Created MPI cluster")

output <- CrossValidate(as.matrix(genotype[,2:10001]),
#output <- CrossValidate(as.matrix(genotype[,2:ncol(genotype)]),
                        phenotype$y,
                        nFolds = 3,
                        Epis = "no",
                        prior = "gaussian",
                        search = "global")

print("Finished CV")

saveRDS(output,"LooserSubsetParCV_5-3-2018.RDS")

print("Wrote output to disk")

closeCluster(cl)

print("Closed cluster")
