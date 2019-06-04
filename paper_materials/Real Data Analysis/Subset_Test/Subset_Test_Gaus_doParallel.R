rm(list=ls())
library("EBEN")
library(parEBEN)
library(readr)
library(dplyr)

phenotype <- read_delim("pheno1", "\t", col_names = c("y"))
genotype <- read_delim("filter_matrix","\t", col_names = c("SAMID",paste0("x",seq(1:5356))))


## Register your local cluster.
library(doParallel)
no_cores <- detectCores()
#cl <- makeCluster(no_cores-1, outfile="Log.txt")
cl <- makeCluster(no_cores)
clusterExport(cl, c("CrossValidate"))
registerDoParallel(cl)
#stopImplicitCluster()
paste("Using", no_cores, "cores.")

output <- CrossValidate(as.matrix(genotype[,2:ncol(genotype)]),
                        phenotype$y,
                        nFolds = 3,
                        Epis = "no",
                        prior = "gaussian",
                        search = "global")

saveRDS(output,"SubsetParCV_4-18-2018_cvoutput.RDS")

model <- EBelasticNet.Gaussian(as.matrix(genotype[,2:ncol(genotype)]),
                               phenotype$y,
                               lambda = output$lambda.optimal,
                               alpha = output$alpha.optimal,
                               Epis = "no")

saveRDS(model,"SubsetParCV_4-18-2018_model.RDS")


