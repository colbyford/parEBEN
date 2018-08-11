#' @title Test Each Model Using Hyperparameter Combinations
#'
#' @export TestModel
#############################

TestModel <- function(BASIS, Target, lambda, alpha, nFolds, foldId = 0, Epis = "no", prior = "gaussian"){
  out <- foreach (i = 1:nFolds, .combine = rbind, .packages = c("EBEN", "parEBEN")) %dopar% {
    #cat("Testing Fold", j, "\n")
    foldId <- AssignToFolds(BASIS, nFolds)

    index <- which(foldId!=i)
    Basis.Train <- BASIS[index,]
    Target.Train <- Target[index]

    index <- which(foldId == i)
    Basis.Test <- BASIS[index,]
    Target.Test <- Target[index]

    if(prior == "gaussian"){
      fit <- EBEN::EBelasticNet.Gaussian(Basis.Train, Target.Train, lambda, alpha, Epis)
	  #fit <- parEBEN.Gaussian(Basis.Train, Target.Train, lambda, alpha, Epis)
      FoldError <- GetFoldError(Basis.Test, Target.Test, fit, prior = "gaussian")
      out <- data.frame(foldId = i,
                        alpha = alpha,
                        lambda = lambda,
                        MSE = FoldError)
    }else{
      fit <- EBEN::EBelasticNet.Binomial(Basis.Train, Target.Train, lambda, alpha, Epis)
	  #fit <- parEBEN.Binomial(Basis.Train, Target.Train, lambda, alpha, Epis)
      FoldError <- GetFoldError(Basis.Test, Target.Test, fit, prior = "binomial")
      out <- data.frame(foldId = i,
                        alpha = alpha,
                        lambda = lambda,
                        logL = FoldError)
    }
    return(out)
  }
  return(out)
}
