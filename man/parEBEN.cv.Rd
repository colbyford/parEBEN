\name{parEBEN.cv}
\alias{parEBEN.cv}

\title{Parallelized Cross-Validation Function to Determine Hyperparameter of the Empirical Bayesian Elastic Net Algorithm.}

\description{This function serves as the front end for the internal parallelized cross-validation functions of the parEBEN package. parEBEN.cv allows for the user to define the parallel environment, the prior, and all other parameters from a single function call.

In the EBEN algorithm, the hyperparameters control the degree of shrinkage, and are obtained via cross-validation. This program calculates the optimal alpha and lambda values, which minimize the mean square error of the prediction.}

\usage{
parEBEN.cv(BASIS, Target, nFolds, Epis = "no", foldId = 0, parMethod = "doParallel", prior = "binomial")
}

\arguments{
  \item{BASIS}{Sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{nFolds}{Number of n-fold cv}
  \item{Epis}{Epistasis - "yes" or "no" for including two-way interactions}
  \item{foldId}{Randomly assign samples to different folds}
  \item{parMethod}{Currently, only "doParallel" is supported}
  \item{prior}{Model prior - "binomial" of "gaussian"}
}

\details{
  If Epis="yes", the program adds two-way interaction K*(K-1)/2 more columns to BASIS
}

\value{
  \item{CrossValidation}{Dataframe of CV results: alpha, lambda, MSE, and StdErr values}
  \item{Optimal_Stats}{The optimal set of hyperparameters with minimum MSE}
  \item{Lambda_optimal}{The optimal lambda hyperparameter value as computed}
  \item{Alpha_optimal}{The optimal alpha hyperparameter value as computed}
}

\references{Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79}

\author{Colby Ford; Dept. of Bioinformatics and Genomics, The University of North Carolina at Charlotte}

\examples{
## To use doParallel funtionality on a multicore environment, you must register your local cluster.
library(doParallel)  
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#clusterExport(cl, c("parEBEN.cv.doParallel"))
registerDoParallel(cl)

## Load in data and required EBEN and parEBEN packages
library(EBEN)
library(parEBEN)
data(BASISbinomial)
data(yBinomial)

## Create small sample matrix for testing
n <- 50
k <- 100
N <- length(yBinomial)
set.seed(1)
set <- sample(N,n)
BASIS <- BASISbinomial[set,1:k]
y <- yBinomial[set]

cv <- parEBEN.cv(BASIS, y, nFolds = 3,Epis = "no", parMethod = "doParallel", prior = "binomial")

#Use the optimal values in the EBEN model
EBENoutput <- EBelasticNet.Binomial(BASIS, y, lambda = cv$Lambda_Optimal, alpha = cv$Alpha_Optimal, Epis = "no",verbose = 5)
}

