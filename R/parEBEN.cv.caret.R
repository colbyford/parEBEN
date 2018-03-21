library(caret)

EBENgrid <- function(x, y, len = NULL, search = "grid"){
  BuildGrid(x, y, Epis = "no")
}

EBENfit_Gaussian <- function(x, y, param, wts, lev, last, weights, classProbs, ...){
  x <- as.matrix(x)
  y <- y
  out <- parEBEN.Gaussian(x, y, param$lambda, param$alpha, Epis = "no")
  out
}

EBENpredict_Gaussian <- function(modelFit, x, y, submodels = NULL){
  x <- as.matrix(x)
  y <- y
  out <- EBelasticNet.Gaussian(x, y, modelFit$hyperparameters[1], modelFit$hyperparameters[2])
  out
}

EBENmodel_Gaussian  <- list(label = "Gaussian EBEN Model",
                            type = "Regression",
                            loop = NULL,
                            library = "EBEN",
                            parameters = data.frame(parameter = c("lambda", "alpha"),
                                                    class = c("numeric", "numeric"),
                                                    label = c("lambda", "alpha")),
                            grid = EBENgrid,
                            fit = EBENfit_Gaussian,
                            prob = NULL,
                            predict = EBENpredict_Gaussian)