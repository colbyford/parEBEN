\name{EBEN-package}
\alias{EBEN-package}
\alias{EBEN}
\docType{package}
\title{Empirical Bayesian Elastic Net (EBEN)}

\description{ Fast EBEN algorithms. \cr
			EBEN implements a normal and generalized gamma hierearchical priors. \cr
				(    **    )	Two parameters (alpha, lambda) are equivalent with elastic net priors. 	\cr
				(    **    )	When parameter alpha = 1, it is equivalent with EBlasso-NE (normal + exponential) \cr					
			Two models are available for both methods: 			\cr
				(    **    )	General linear regression model.			\cr
				(    **    )	Logistic regression model.			\cr
			Multi-collinearity: 	\cr
				(    **    )	for group of high correlated or collinear variables: EBEN identifies the group of variables estimates their effects together. \cr
				(    **    )	group of variables can be selected together. \cr
			*Epistasis (two-way interactions) can be included for all models/priors \cr
			*model implemented with memory efficient c code. \cr
			*LAPACK/BLAS are used for most linear algebra computations. 
}
\details{
\tabular{ll}{
Package: \tab EBEN\cr
Type: \tab Package\cr
Version: \tab 4.6\cr
Date: \tab 2015-10-06\cr
License: \tab gpl\cr
}

}
\author{
Anhui Huang
}
\references{
key algorithms: \cr
	Cai, X., Huang, A., and Xu, S. (2011). Fast empirical Bayesian LASSO for multiple quantitative trait locus mapping. BMC Bioinformatics 12, 211.\cr
	Huang A, Xu S, Cai X. (2013). Empirical Bayesian LASSO-logistic regression for multiple binary trait locus mapping. BMC genetics  14(1):5. \cr
	Huang, A., Xu, S., and Cai, X. (2014). Empirical Bayesian elastic net for multiple quantitative trait locus mapping. Heredity 10.1038/hdy.2014.79 \cr
Other publications: \cr
	Huang, A., E. Martin, et al. (2014). "Detecting genetic interactions in pathway-based genome-wide association studies." Genet Epidemiol 38(4): 300-309.\cr
	Huang, A., S. Xu, et al. (2014). "Whole-genome quantitative trait locus mapping reveals major role of epistasis on yield of rice." PLoS ONE 9(1): e87330.\cr
	Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).
}
\keyword{ package }
