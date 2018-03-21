GetFoldError <- function(BASIS, Target, fit, prior = "gaussian"){
  if(prior == "gaussian"){
    M	<- length(fit$weight)/6
    Betas <- matrix(fit$weight,nrow= M,ncol =6, byrow= FALSE)
    Mu <- Betas[,3]
    Mu0 <- fit$Intercept[1]
    if(is.na(Mu0)){break}
    
    ntest <- nrow(BASIS)
    basisTest <- matrix(rep(0,ntest*M),ntest,M)
    
    for(i_basis in 1:M){
      loc1 <- Betas[i_basis,1]
      loc2 <- Betas[i_basis,2]
      if(loc1 !=0){
        if(loc1==loc2){basisTest[,i_basis] <- BASIS[,loc1]}
        else{basisTest[,i_basis] <- BASIS[,loc1] * BASIS[,loc2]}
      }else{
        basisTest<- rep(0,length(Target))
      }						
    }
    
    #compute mean square error:
    temp <- Target - (Mu0 + basisTest%*%Mu)				
    MeanSqErr <- t(temp)%*%temp
    return(MeanSqErr)
  }else{
    M	<- length(fit$weight)/6
    Betas <- matrix(fit$weight, nrow = M, ncol = 6, byrow = FALSE)
    Mu <- Betas[,3]
    Mu0 <- fit$Intercept[1]
    
    rm(list="fit")
    ntest <- nrow(BASIS)
    #M <- nrow(Betas)
    if(M==1 && Betas[1,1]== 0){
      logL <- 0
    }else{
      basisTest <- matrix(rep(0,ntest*M),ntest,M)
      for(i_basis in 1:M){
        loc1 <- Betas[i_basis,1]
        loc2 <- Betas[i_basis,2]
        if(loc1==loc2){basisTest[,i_basis] <- BASIS[,loc1]}
        else{basisTest[,i_basis] <- BASIS[,loc1] * BASIS[,loc2]}
      }
      temp <- exp(Mu0 + basisTest%*%Mu)
      if(max(temp) > 1e10) temp[which(temp>1e10)] <- 1e5
      if(min(temp) < 1e-10) temp[which(temp<1e-10)] <- 1e-5
      logL <- mean(Target*log(temp/(1+temp)) + (1-Target)*log(1/(1+temp)))
      }
  return(logL)
    }
}