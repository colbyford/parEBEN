\name{EBelasticNet.BinomialCV}
\alias{EBelasticNet.BinomialCV}

\title{ Cross Validation (CV) Function to Determine Hyperparameter of the EB_Elastic Net Algorithm for Binomial Model with Normal-Gamma (NG) Prior Distribution}
\description{Hyperparameter controls degree of shrinkage, and is obtained via Cross Validation (CV). This program calculates the maximum lambda
			that allows one non-zero basis; and performs a search down to 0.001*lambda_max at even steps. (20 steps)}
\usage{
EBelasticNet.BinomialCV(BASIS, Target, nFolds, Epis = "no",foldId)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{nFolds}{number of n-fold cv}
  \item{Epis}{"yes" or "no" for including two-way interactions}
  \item{foldId}{random assign samples to different folds}
}
\details{If Epis="yes", the program adds two-way interaction K*(K-1)/2 more columns to BASIS}
\value{
  \item{CrossValidation}{col1: hyperparameter; col2: loglikelihood mean; standard ERROR of nfold mean log likelihood}
  \item{Lmabda_optimal}{the optimal hyperparameter as computed}
    \item{Alpha_optimal}{the optimal hyperparameter as computed}
}
\references{Huang A, Xu S, Cai X: Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics 2013, 14(1):5.}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}



%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(EBEN)
data(BASISbinomial)
data(yBinomial)
#reduce sample size to speed up the running time
n = 50;
k = 100;
N = length(yBinomial);
set.seed(1)
set  = sample(N,n);
BASIS = BASISbinomial[set,1:k];
y  = yBinomial[set];
nFolds = 3
CV = EBelasticNet.BinomialCV(BASIS, y, nFolds = 3,Epis = "no")



}

