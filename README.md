# parEBEN - A parallel implementation of the Empirical Bayesian Elastic Net in R
[![Build Status](https://travis-ci.org/colbyford/parEBEN.svg?branch=master)](https://travis-ci.org/colbyford/parEBEN)

## Abstract

The Empirical Bayesian Elastic Net (EBEN) algorithm was developed by [Huang et al.](https://www.nature.com/articles/hdy201479) for handling multicollinearity in generalized linear regression models. Historically, this has been used in the analysis of quantitative trait loci (QTLs) and gene-gene interactions (epistasis). In addition to the algorithm, the group also created the [EBEN package for R](https://cran.r-project.org/package=EBEN). This package includes functions to generate the elastic nets for both binomial and gaussian models. These functions are efficient and do not require large amounts of computational time. However, the package also includes functions for the cross-validation of those models. While essential, this step is a considerably more complex task. The cross-validation functions perform a sweep to determine hyperparameters and minimize prediction error. More specifically, an n-fold cross-validation sweep is performed to minimize error by trying combinations of two parameters (λ<sub>1</sub> and λ<sub>2</sub>) in a stepped manner. Experimentally, it has been shown that this can an extended amount of time, especially on larger dataset such in genomics problems.

To combat this complexity issue, the parallelization of the cross-validation functions was performed by employing parallel packages in R. By parallelizing the iterations of the cross-validation over multiple CPU cores or multiple machines of a computing clusters, a drastic time reduction can seen with no negative effect on the resulting EBEN models. By reducing the computation time, regression models on larger, more complex data can be completed without such a delay. This also opens the door for larger genetic datasets to be analyzed as opposed to limiting the research due to time and computing resource constraints. Thus, parallelizing the cross-validation of the EBEN models will prove to be greatly beneficial in future research using cross-validated Bayesian elastic nets.

## Installation

You can install the latest stable version from CRAN using the following command:
```r
library(devtools)
install_github("colbyford/parEBEN")
library(parEBEN)
```

## To Do List

- [x] Binomial prior cross-validation script with doParallel.
- [x] Gaussian prior cross-validation script with doParallel.
- [ ] Binomial prior cross-validation script with Microsoft ML Server (RevoScaleR).
- [ ] Gaussian prior cross-validation script with Microsoft ML Server (RevoScaleR).
- [ ] Binomial prior cross-validation script with SparkR.
- [ ] Gaussian prior cross-validation script with SparkR.
- [ ] Binomial prior cross-validation script with CUDA.
- [ ] Gaussian prior cross-validation script with CUDA.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
