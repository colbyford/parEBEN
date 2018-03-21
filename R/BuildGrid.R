GetLambdaMax <- function (BASIS, Target, nFolds, Epis = "no", foldId = 0){
  N <- nrow(BASIS)
  K <- ncol(BASIS)
  set.seed(1)
  if(length(foldId)!=N){
    if(N%%nFolds!=0){
      foldId <- sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N)
    }else{
      foldId <- sample(rep(1:nFolds,floor(N/nFolds)),N)
    }
  }
  lambda_Max <- log(1.1)
  response <- Target-mean(Target)
  response <- response/sqrt(sum(response*response))
  
  for(i_b in 1:K){
    basis <- BASIS[,i_b]
    basis <- basis/sqrt(sum(basis*basis))
    corBy <- basis%*%response
    if(corBy>lambda_Max) lambda_Max = corBy
  }	
  
  if(Epis == "yes"){
    for(i_b in 1:(K-1)){
      for(i_bj in (i_b + 1):K){
        basis <- BASIS[,i_b]*BASIS[,i_bj]
        basis <- basis/sqrt(sum(basis*basis))
        corBy <- basis%*%(Target-mean(Target))
        if(corBy>lambda_Max) lambda_Max = corBy
      }
    }
  }
return(as.vector(lambda_Max))  
}

BuildGrid <- function(BASIS, Target, nFolds, Epis = "no", foldId = 0){
  lambda_Max <- GetLambdaMax(BASIS, Target, nFolds, Epis, foldId)
  lambda_Max <- lambda_Max * 10
  lambda_Min <- log(0.001 * lambda_Max)
  step <- (log(lambda_Max) - lambda_Min)/19
  Lambda <- exp(seq(from = log(lambda_Max), to = lambda_Min, by = -step))
  #lambda <- 1 #Initialize lambda
  N_step <- length(Lambda)
  
  step <- 1;
  Alpha <- seq(from = 1, to = 0.05, by = -0.05)
  #alpha <- 0.5 #Initialize alpha
  nAlpha <- length(Alpha);
  
  grid <- merge(Alpha,Lambda,all=TRUE)
  
  #MSEcv <- mat.or.vec((N_step*nAlpha),4)
  #MSEeachAlpha <- mat.or.vec(nAlpha,4) # minimum MSE for each alpha
  #MeanSqErr <- mat.or.vec(nFolds,1)
  #SSE1Alpha	<- matrix(1e10,N_step,2)# temp matrix to keep MSE + std in each step
  
  return(grid)
}