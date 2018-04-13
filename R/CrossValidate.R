CrossValidate <- function(BASIS, Target, nFolds, foldId = 0, Epis = "no", prior = "gaussian", search = "global"){
  if(search == "global"){
    ParameterGrid <- BuildGrid(BASIS, Target, nFolds, Epis)
    
    if(prior == "gaussian"){
      MSE_temp <- foreach (row = 1:nrow(ParameterGrid), .combine = rbind) %dopar% {
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
      logL_temp <- foreach (row = 1:nrow(ParameterGrid), .combine = rbind) %dopar% {
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