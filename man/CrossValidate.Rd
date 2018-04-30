\name{CrossValidate}
\alias{CrossValidate}

\title{Parallel Implementations of the Empirical Bayesian Elastic Net Cross-Validation}

\description{This function serves as the front end for the internal parallelized cross-validation functions of the parEBEN package. `CrossValidate` allows for the user to define the input data, number of cross-validation folds, the prior, and all other parameters from a single function call.

In the EBEN algorithm, the hyperparameters control the degree of shrinkage, and are obtained via cross-validation. This program calculates the optimal alpha and lambda values, which minimize the mean square error of the prediction.}

\usage{
CrossValidate(BASIS, Target, nFolds, foldId = 0, Epis = "no", prior = "gaussian", search = "global")
}

\arguments{
  \item{BASIS}{Sample matrix; rows correspond to samples, columns correspond to features. Can be a matrix or data frame.}
  \item{Target}{Dependent variable of each individual. Usually a vector}
  \item{nFolds}{Number of n-fold cv}
  \item{Epis}{Epistasis - "yes" or "no" for including two-way interactions}
  \item{foldId}{Randomly assign samples to different folds}
  \item{prior}{Model prior - "binomial" or "gaussian"}
  \item{search}{Search type - "global" or "local"}
}

\details{
  \item{1}{If Epis = "yes", the program adds two-way interaction K*(K-1)/2 more columns to BASIS}
  \item{2}{If search = "global", the program will search through all 400 combinations of alpha and lambda. Otherwise, the "local" search is performed using the logic from the original EBEN package, which skips certain iterations based on error changes. The empirical Bayesian elastic net is fully parallelized when performing a "global" search, which may take longer than the serial EBEN package as more iterations may be computed. For the "local" search, only the cross-validation steps are parallelized whereas the hyperparameter sweeping is still serial.}
}

\value{
  \item{Results.Detail}{Dataframe of CV results: foldId, alpha, lambda, and MSE.}
  \item{Results.Fold}{Summary of CV by alpha and lambda. Reports SSE.}
  \item{lambda.optimal}{The optimal lambda hyperparameter value as computed.}
  \item{alpha.optimal}{The optimal alpha hyperparameter value as computed.}
}

\references{Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79}

\author{Colby Ford; Dept. of Bioinformatics and Genomics, The University of North Carolina at Charlotte}

\examples{
## To use doParallel or doMPI funtionality on a multicore/multi-machine environment, you must register your local cluster.
library(doParallel)  #or library(doMPI)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#clusterExport(cl, c("CrossValidate"))
registerDoParallel(cl)

## Load in data and required EBEN and parEBEN packages
library(EBEN)
library(parEBEN)

## Create small sample matrix for testing
data(BASIS)
data(y)
n = 50
k = 100
BASIS = BASIS[1:n,1:k]
y  = y[1:n]

parEBENcv <- CrossValidate(BASIS, y, nFolds = 3, Epis = "no", prior = "gaussian", search = "global")

#Use the optimal values in the EBEN model
EBENoutput <- EBelasticNet.Gaussian(BASIS, y, lambda = parEBENcv$lambda.optimal, alpha = cv$alpha.optimal, Epis = "no", verbose = 1)
}

