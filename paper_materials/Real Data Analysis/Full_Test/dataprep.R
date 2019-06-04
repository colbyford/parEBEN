## MUST LOAD IN DATA AS MATRICES

BASIS <- as.matrix(read.delim("filter_matrix_looser_0.02_main_0.15_epi-19872.tsv", sep = "\t")[,-1])
y <- as.matrix(read.delim("pheno1.tsv", sep = "\t"))

## THEN USE parEBEN or EBEN

parEBENcv <- CrossValidate(BASIS, y, nFolds = 3, Epis = "no", prior = "gaussian", search = "global")
EBENoutput <- EBelasticNet.Gaussian(BASIS, y, lambda = parEBENcv$lambda.optimal, alpha = parEBENcv$alpha.optimal, Epis = "no", verbose = 1)
