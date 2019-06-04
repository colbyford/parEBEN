rm(list=ls())
library("EBEN")
library(readr)
library(dplyr)

phenotype <- read_delim("pheno1", "\t", col_names = c("y"))
genotype <- read_delim("filter_matrix","\t", col_names = c("SAMID",paste0("x",seq(1:5356))))


CVoutput <- readRDS("Subset_4-15-2018_parCV.RDS")

model <- EBelasticNet.Gaussian(as.matrix(genotype[,2:ncol(genotype)]),
                               phenotype$y,
                               lambda = CVoutput$lambda.optimal,
                               alpha = CVoutput$alpha.optimal,
                               Epis = "no")

saveRDS(model,"Subset_4-15-2018_model.RDS")


