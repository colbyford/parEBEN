TestModel <- function(BASIS, Target, lambda, alpha, nFolds, foldId = 0, Epis = "no", prior = "gaussian"){
  out <- foreach (i = 1:nFolds, .combine = rbind, .packages = c("parEBEN","EBEN")) %dopar% {
    #cat("Testing Fold", j, "\n")
    foldId <- AssignToFolds(BASIS, nFolds)
    
    index <- which(foldId!=i)
    Basis.Train <- BASIS[index,]
    Target.Train <- Target[index]
    
    index <- which(foldId == i)
    Basis.Test <- BASIS[index,]
    Target.Test <- Target[index]
    
    if(prior == "gaussian"){
      fit <- parEBEN.Gaussian(Basis.Train, Target.Train, lambda, alpha, Epis)
      FoldError <- GetFoldError(fit,prior)
      out <- data.frame(foldId = i,
                        alpha = alpha,
                        lambda = lambda,
                        MSE = FoldError)
    }else{
      fit <- parEBEN.Binomial(Basis.Train, Target.Train, lambda, alpha, Epis)
      FoldError <- GetFoldError(fit,prior)
      out <- data.frame(foldId = i,
                        alpha = alpha,
                        lambda = lambda,
                        logL = FoldError)
    }
    return(out)
  }
  return(out)
}