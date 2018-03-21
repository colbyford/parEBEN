CrossValidate <- function(BASIS, Target, nFolds, foldId = 0, Epis = "no", prior = "gaussian"){
  ParameterGrid <- BuildGrid(BASIS, Target, nFolds, Epis)
  
  if(prior == "gaussian"){
    MSE <- foreach (row = 1:nrow(ParameterGrid), .combine = rbind) %dopar% {
      alpha <- ParameterGrid$alpha[row]
      lambda <- ParameterGrid$lambda[row]
      TestModel(BASIS, Target, lambda, alpha, nFolds, foldId, Epis, prior = "gaussian")
      }
    
    SSE <- MSE %>%
      group_by(alpha, lambda) %>% 
      summarise(SSE = mean(MSE)/(sd(MSE)/sqrt(max(foldId))))
    
    index <- which.min(SSE$SSE)
    lambda.optimal <- SSE[index,]$lambda
    alpha.optimal <- SSE[index,]$alpha
    
    out <- list(Results.Detail = MSE,
                Results.ByFold = SSE,
                lambda.optimal = lambda.optimal,
                alpha.optimal = alpha.optimal)
    
  }else{
    logL <- foreach (row = 1:nrow(ParameterGrid), .combine = rbind) %dopar% {
      alpha <- ParameterGrid$alpha[row]
      lambda <- ParameterGrid$lambda[row]
      TestModel(BASIS, Target, lambda, alpha, nFolds, foldId, Epis, prior = "binomial")
    }
    SSE <- logL %>% 
      group_by(alpha, lambda) %>% 
      summarise(SSE = mean(logL)/(sd(logL)/sqrt(max(foldId))))
    
    index <- which.min(SSE$SSE)
    lambda.optimal <- SSE[index,]$lambda
    alpha.optimal <- SSE[index,]$alpha
    
    out <- list(Results.Detail = logL,
                Results.ByFold = SSE,
                lambda.optimal = lambda.optimal,
                alpha.optimal = alpha.optimal)
  }
  
  return(out)
}