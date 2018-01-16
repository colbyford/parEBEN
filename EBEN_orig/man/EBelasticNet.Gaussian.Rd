\name{EBelasticNet.Gaussian}
\alias{EBelasticNet.Gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The EB Elastic Net Algorithm for Gaussian Model}
\description{General linear regression, normal-Gamma (NG) hierarchical prior for regression coefficients}
\usage{EBelasticNet.Gaussian(BASIS, Target, lambda, alpha,Epis = "no",verbose = 0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{sample matrix; rows correspond to samples, columns correspond to features}
  \item{Target}{Response each individual}
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
	\item{WaldScore}{Wald Score}
	\item{Intercept}{Intercept}
	\item{lambda}{the hyperparameter; same as input lambda}
	\item{alpha}{the hyperparameter; same as input alpha}	
}
\references{	Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79 \cr
}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(EBEN)
data(BASIS)
data(y)
n = 50;
k = 100;
BASIS = BASIS[1:n,1:k];
y  = y[1:n];
Blup = EBelasticNet.Gaussian(BASIS, y,lambda = 0.0072,alpha = 0.95, Epis = "no",verbose = 0)
betas 			= Blup$weight
betas
}
