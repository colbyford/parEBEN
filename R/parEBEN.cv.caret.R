library(caret)

EBENgrid <- function(x, y, len = NULL, search = "grid"){
  N <- nrow(x)
  K <- ncol(x)
  set.seed(1)
  #if(length(foldId)!=N){
  #  if(N%%nFolds!=0){
  #    foldId <- sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N)
  #  }else{
  #    foldId <- sample(rep(1:nFolds,floor(N/nFolds)),N)
  #  }
  #}
  lambda_Max <- log(1.1)
  response <- y-mean(y)
  response <- response/sqrt(sum(response*response))
  
  for(i_b in 1:K){
    basis <- x[,i_b]
    basis <- basis/sqrt(sum(basis*basis))
    corBy <- basis%*%response
    if(corBy>lambda_Max) lambda_Max = corBy
  }	
  
  #if(Epis == "yes"){
  #  for(i_b in 1:(K-1)){
  #    for(i_bj in (i_b + 1):K){
  #      basis <- x[,i_b]*x[,i_bj]
  #      basis <- basis/sqrt(sum(basis*basis))
  #      corBy <- basis%*%(y-mean(y))
  #      if(corBy>lambda_Max) lambda_Max = corBy
  #    }
  #  }
  #}
  
  lambda_Max <- as.vector(lambda_Max) * 10
  lambda_Min <- log(0.001 * lambda_Max)
  step <- (log(lambda_Max) - lambda_Min)/19
  Lambda <- exp(seq(from = log(lambda_Max), to = lambda_Min, by = -step))
  #lambda <- 1 #Initialize lambda
  N_step <- length(Lambda)
  
  step <- 1
  Alpha <- seq(from = 1, to = 0.05, by = -0.05)
  #alpha <- 0.5 #Initialize alpha
  nAlpha <- length(Alpha)
  
  #grid <- merge(Alpha,Lambda,all=TRUE)
  #colnames(grid) <- c("alpha","lambda")
  grid <- expand.grid(alpha = Alpha, lambda = Lambda)
  
  grid
}

EBENfit_Gaussian <- function(x, y, param, wts, lev, last, weights, classProbs, ...){
  x <- as.matrix(x)
  y <- y
  EBelasticNet.Gaussian(x, y, param$lambda, param$alpha)
  
}

EBENpredict_Gaussian <- function(modelFit, x, y){
  #predict(modelFit, x)
  #print(modelFit)
  #EBelasticNet.Gaussian(x, y, modelFit$hyperparameters[1], modelFit$hyperparameters[2])
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