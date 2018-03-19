parEBEN.cv <- function (BASIS, Target, nFolds, Epis = "no", foldId = 0, parMethod = "doParallel", prior = "binomial"){
  if(parMethod == "doParallel"){
    result <- parEBEN.cv.doParallel(BASIS, Target, nFolds, Epis, foldId, prior)
    return(result)
  } else if(parMethod == "doMPI"){
    result <- parEBEN.cv.doMPI(BASIS, Target, nFolds, Epis, foldId, prior)
    return(result)
  }
  else {
    cat("Currently, only the `doParallel` and `doMPI` function sets for parallelization are available.")
  }
} 