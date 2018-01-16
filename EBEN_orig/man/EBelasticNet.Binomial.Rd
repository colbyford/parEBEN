\name{EBelasticNet.Binomial}
\alias{EBelasticNet.Binomial}
\title{The EB Elastic Net Algorithm for Binomial Model with Normal-Gamma(NG) Prior Distribution}
\description{Generalized linear regression, normal-Gxponential (NG) hierarchical prior for regression coefficients}
\usage{EBelasticNet.Binomial(BASIS, Target, lambda, alpha,Epis = "no",verbose = 0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Class label of each individual, TAKES VALUES OF 0 OR 1}
  \item{lambda}{Hyperparameter controls degree of shrinkage; can be obtained via Cross Validation; lambda>0}
  \item{alpha}{Hyperparameter controls degree of shrinkage; can be obtained via Cross Validation; 0<alpha<1}
  \item{Epis}{"yes" or "no" for including two-way interactions}
  \item{verbose}{0 or 1; 1: display message; 0 no message}
}
\details{If Epis="yes", the program adds two-way interaction of K*(K-1)/2 more columns to BASIS}
\value{
	\item{weight}{the none-zero regression coefficients: \cr 
	col1,col2 are the indices of the bases(main if equal); \cr
	col3: coefficent value; \cr
	col4: posterior variance;\cr
	col5: t-value;\cr
	col6: p-value}
	\item{logLikelihood}{log likelihood from the final regression coefficients}
	\item{WaldScore}{Wald Score}
	\item{Intercept}{Intercept}
	\item{lambda}{the hyperparameter; same as input lambda}
	\item{alpha}{the hyperparameter; same as input alpha}	
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
set  = sample(N,n);
BASIS = BASISbinomial[set,1:k];
y  = yBinomial[set];
output = EBelasticNet.Binomial(BASIS, y,lambda = 0.1,alpha = 0.5, Epis = "no",verbose = 5)
}
