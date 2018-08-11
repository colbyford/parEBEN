#' @title Assign Data to Folds
#'
#' @export AssignToFolds
#############################

AssignToFolds <- function(BASIS, nFolds = 0, foldId = 0){
  N <- nrow(BASIS)
  K <- ncol(BASIS)
  set.seed(1)
  if(length(foldId)!=N)
  {
    if(N%%nFolds!=0){
      foldId <- sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N)
    }else{
      foldId <- sample(rep(1:nFolds,floor(N/nFolds)),N)
    }
  }
  return(foldId)
}
