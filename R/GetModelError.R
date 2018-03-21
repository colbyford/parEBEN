GetModelError <- function(x, y, lambda, alpha, Epis = "no", prior = "gaussian"){
  if(prior == "gaussian"){
    fit <- parEBEN.Gaussian(x, y, lambda, alpha, Epis, verbose=0)
    M	<- length(fit$weight)/6
    Betas <- matrix(fit$weight, nrow = M, ncol = 6, byrow= FALSE)
    Mu <- Betas[,3]
    Mu0 <- fit$Intercept[1]
    if(is.na(Mu0)){break}
    ntest <- nrow(x)
    xTest <- matrix(rep(0,ntest*M),ntest,M)
    for(i_basis in 1:M){
      loc1 <- Betas[i_basis,1]
      loc2 <- Betas[i_basis,2]
      if(loc1 !=0){
        if(loc1==loc2){xTest[,i_basis] <- x[,loc1]}
        else{xTest[,i_basis] <- x[,loc1] * x[,loc2]}
        }else{
          xTest <- rep(0,length(y))
        }						
    }
    #Compute MSE:
    temp <- y - (Mu0 + xTest%*%Mu)				
    MeanSqErr[i] <- t(temp)%*%temp
  }else{
    fit <- parEBEN.Binomial(x, y, lambda, alpha, Epis, verbose=0)
    M	<- length(fit$weight)/6
    Betas <- matrix(fit$weight, nrow= M,ncol =6, byrow= FALSE)
    Mu <- Betas[,3]
    Mu0 <- fit$Intercept[1]
    
    rm(list="fit")
    ntest <- nrow(Basis.Test)
    if(M==1 && Betas[1,1]== 0){
      logL[i] <- 0
      }else{
      xTest <- matrix(rep(0,ntest*M),ntest,M)
      for(i_basis in 1:M){
        loc1 <- Betas[i_basis,1]
        loc2 <- Betas[i_basis,2]
        if(loc1==loc2){xTest[,i_basis] <- x[,loc1]
        }else{
        xTest[,i_basis] <-  x[,loc1]* x[,loc2]
        }
      }
      
      temp <- exp(Mu0 + xTest%*%Mu)
      if(max(temp)>1e10) temp[which(temp>1e10)] <- 1e5
      if(min(temp)<1e-10) temp[which(temp<1e-10)] <- 1e-5
      logL[i] <- mean(y*log(temp/(1+temp)) + (1-y)*log(1/(1+temp)))
      }
  return(logL)
    }
}