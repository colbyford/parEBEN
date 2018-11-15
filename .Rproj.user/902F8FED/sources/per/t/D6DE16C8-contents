#' @title Perform Hyperparameter Sweep and Cross Validate Empirical Bayesian Elastic Net Model
#'
#' @description This function serves as the front end for the internal parallelized cross-validation functions of the parEBEN package. `CrossValidate` allows for the user to define the input data, number of cross-validation folds, the prior, and all other parameters from a single function call.
#' In the EBEN algorithm, the hyperparameters control the degree of shrinkage, and are obtained via cross-validation. This program calculates the optimal alpha and lambda values, which minimize the mean square error of the prediction.
#'
#' @usage CrossValidate(BASIS, Target, nFolds, foldId = 0, Epis = "no", prior = "gaussian", search = "global")
#'
#' @param BASIS Sample matrix; rows correspond to samples, columns correspond to features. Can be a matrix or data frame.
#' @param Target Dependent variable of each individual. Usually a vector
#' @param nFolds Number of folds n-fold cross validation.
#' @param Epis Epistasis - \code{"yes"} or \code{"no"} for including two-way interactions
#' @param foldId Randomly assign samples to different folds
#' @param prior Model prior - \code{"binomial"} or \code{"gaussian"}
#' @param search Search type - \code{"global"} or \code{"local"}
#'
#' @details   If Epis = \code{"yes"}, the program adds two-way interaction K*(K-1)/2 more columns to BASIS
#' If search = \code{"global"}, the program will search through all 400 combinations of alpha and lambda. Otherwise, the \code{"local"} search is performed using the logic from the original EBEN package, which skips certain iterations based on error changes.
#' The empirical Bayesian elastic net is fully parallelized when performing a \code{"global"} search, which may take longer than the serial EBEN package as more iterations may be computed.
#' For the \code{"local"} search, only the cross-validation steps are parallelized whereas the hyperparameter sweeping is still serial.
#'
#' @return Four objects are returned from the CV exercise:
#' \item{Results.Detail}{Dataframe of CV results: foldId, alpha, lambda, and MSE.}
#' \item{Results.Fold}{Summary of CV by alpha and lambda. Reports SSE.}
#' \item{lambda.optimal}{The optimal lambda hyperparameter value as computed.}
#' \item{alpha.optimal}{The optimal alpha hyperparameter value as computed.}
#'
#' @references Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79
#' @author Colby T. Ford, Ph.D.; Dept. of Bioinformatics and Genomics, The University of North Carolina at Charlotte
#'
#' @examples ## To use doParallel or doMPI funtionality on a multicore/multi-machine environment, you must register your local cluster.
#' library(doParallel)  #or library(doMPI)
#' no_cores <- detectCores() - 1
#' cl <- makeCluster(no_cores)
#' #clusterExport(cl, c("CrossValidate"))
#' registerDoParallel(cl)
#'
#' ## Load in data and required EBEN and parEBEN packages
#' library(EBEN)
#' library(parEBEN)
#'
#' ## Create small sample matrix for testing
#' data(BASIS)
#' data(y)
#' n = 50
#' k = 100
#' BASIS = BASIS[1:n,1:k]
#' y  = y[1:n]
#'
#' parEBENcv <- CrossValidate(BASIS, y, nFolds = 3, Epis = "no", prior = "gaussian", search = "global")
#'
#' #Use the optimal values in the EBEN model
#' EBENoutput <- EBelasticNet.Gaussian(BASIS, y, lambda = parEBENcv$lambda.optimal, alpha = parEBENcv$alpha.optimal, Epis = "no", verbose = 1)
#'
#' @import EBEN
#' @import foreach
#' @import dplyr
#' @import magrittr
#' @export CrossValidate
###############################

CrossValidate <- function(BASIS, Target, nFolds, foldId = 0, Epis = "no", prior = "gaussian", search = "global"){
  if(search == "global"){
    ParameterGrid <- BuildGrid(BASIS, Target, nFolds, Epis)

    if(prior == "gaussian"){
      MSE_temp <- foreach (row = 1:nrow(ParameterGrid), .combine = rbind, .packages = "parEBEN", .verbose = TRUE) %dopar% {
        alpha <- ParameterGrid$alpha[row]
        lambda <- ParameterGrid$lambda[row]
        TestModel(BASIS, Target, lambda, alpha, nFolds, foldId, Epis, prior = "gaussian")
      }

      Error <- MSE_temp %>%
        group_by(alpha, lambda) %>%
        summarise(SE = sd(MSE)/sqrt(max(foldId)),
                  MSE = mean(MSE)
        )

      index <- which.min(Error$MSE)
      lambda.optimal <- Error[index,]$lambda
      alpha.optimal <- Error[index,]$alpha

      out <- list(Results.Detail = MSE_temp,
                  Results.Summary = Error,
                  lambda.optimal = lambda.optimal,
                  alpha.optimal = alpha.optimal)

    }else{
      logL_temp <- foreach (row = 1:nrow(ParameterGrid), .combine = rbind, .package = "parEBEN", .verbose = FALSE) %dopar% {
        alpha <- ParameterGrid$alpha[row]
        lambda <- ParameterGrid$lambda[row]
        TestModel(BASIS, Target, lambda, alpha, nFolds, foldId, Epis, prior = "binomial")
      }
      Error <- logL_temp %>%
        group_by(alpha, lambda) %>%
        summarise(SE = sd(logL)/sqrt(max(foldId)),
                  Likelihood = -mean(logL)
        )

      index <- which.min(Error$MSE)
      lambda.optimal <- Error[index,]$lambda
      alpha.optimal <- Error[index,]$alpha

      out <- list(Results.Detail = logL_temp,
                  Results.Summary = Error,
                  lambda.optimal = lambda.optimal,
                  alpha.optimal = alpha.optimal)
    }

    return(out)

  }else{
    out <- LocalSearch(BASIS, Target, nFolds, Epis, foldId, prior)

    return(out)
  }

}
