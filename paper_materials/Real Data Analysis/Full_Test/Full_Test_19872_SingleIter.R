#rm(list=ls())
library("EBEN")
#library(parEBEN)
library(readr)
library(dplyr)

phenotype <- read_delim("pheno1", "\t", col_names = c("y"))
genotype <- read_delim("filter_matrix_looser_0.02_main_0.15_epi-19872","\t", col_names = c("SAMID",paste0("x",seq(1:19871))))

system.time(
output <- EBelasticNet.Gaussian(as.matrix(genotype[,2:ncol(genotype)]),
                                phenotype$y,
                                lambda = 0.5,
                                alpha = 0.5,
                                Epis = "no",
                                verbose = TRUE)
)

saveRDS(output,"FullTest_19872_SingleIter_8-7-2018.RDS")
