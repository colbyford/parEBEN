# *parEBEN* - Parallel Implementations of the Empirical Bayesian Elastic Net Cross-Validation in R
<h3 align = "right">Colby T. Ford, Ph.D.</h3>

<img align="right" src="https://raw.githubusercontent.com/colbyford/parEBEN/master/img/parEBEN_icon.png" alt="parEBEN icon" width="200">

## Abstract

The Empirical Bayesian Elastic Net (EBEN) algorithm was developed by [Huang et al.](https://www.nature.com/articles/hdy201479) for handling multicollinearity in generalized linear regression models. Historically, this has been used in the analysis of quantitative trait loci (QTLs) and gene-gene interactions (epistasis). In addition to the algorithm, the group also created the [EBEN package for R](https://cran.r-project.org/package=EBEN). This package includes functions to generate the elastic nets for both binomial and gaussian priors. These functions are efficient and do not require large amounts of computational time. However, the package also includes functions for the cross-validation of those models. While essential, this step is a considerably more complex task. The cross-validation functions perform a sweep to determine hyperparameters and minimize prediction error. More specifically, an n-fold cross-validation sweep is performed to minimize error by trying combinations of two parameters (α and λ) in a stepped manner. Experimentally, it has been shown that this can take a rather extended amount of time, especially on larger datasets (as seen in genomics problems).

<img align="center" src="https://raw.githubusercontent.com/colbyford/parEBEN/master/img/timebottleneck_nfold.png" alt="CV Bottleneck">

To combat this complexity issue, the parallelization of the cross-validation functions was performed by employing parallel packages in R. By parallelizing the iterations of the cross-validation over multiple CPU cores or multiple machines of a computing clusters, a drastic time reduction can seen with no negative effect on the resulting EBEN models. By reducing the computation time, regression models on larger, more complex data can be completed without such a delay. This also opens the door for larger datasets to be analyzed as opposed to limiting the research due to time and computing resource constraints. Thus, parallelizing the cross-validation of the EBEN models will prove to be greatly beneficial in future research using cross-validated Bayesian elastic nets.

## Time Reduction Benchmark

To interactively view cross-validation time benchmarks between parEBEN and the original EBEN package, click [here](https://public.tableau.com/profile/cford38#!/vizhome/parEBEN-Benchmarks/BinomialCross-Validation)

## Installation

You can install the latest stable version from GitHub using the following command:
```r
library(devtools)
install_github("colbyford/parEBEN")
library(parEBEN)
```

## Usage
First, select the parallelization method you wish to use. Currently, all *foreach*-related methods are supported such as *doParallel*, *doMPI*, and *doSNOW*.
### Initialize The Cluster
Note: Refer to the manual for your desired *foreach* parallelization package as the initialization may differ between methods.
##### Local Parallel
```r
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
#clusterExport(cl, c("CrossValidate"))
registerDoParallel(cl)
```
##### Cluster Distribution
```r
library(doMPI)
# create and register a doMPI cluster if necessary
if (!identical(getDoParName(), 'doMPI')) {
  # set count to (cores_requested-1)
  cl <- startMPIcluster(count=255,verbose=TRUE)
  registerDoMPI(cl)
}
```

##### Microsoft Machine Learning Server Distribution
```
## Set your compute contaxt as Spark, local parallel, MapReduce, etc.
### See: https://docs.microsoft.com/en-us/machine-learning-server/r-reference/revoscaler/rxspark
### Sample Code: https://gist.github.com/premalxyz/e97ae7823052b7a426cb816830c0188c#file-spark_compute_context-r

mySparkCluster <- RxSpark(ClusterInfo)
rxSetComputeContext(mySparkCluster)

## Register the context using doRSR
library(doRSR)
registerDoRSR()
```

### Begin the Cross-Validation
```r
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

parEBENcv <- CrossValidate(BASIS,
                           y,
                           nFolds = 3,
                           Epis = "no",
                           prior = "gaussian",
                           search = "global"
                           )

## Use the optimal values in the EBEN model
EBENoutput <- EBelasticNet.Gaussian(BASIS,
                                    y,
                                    lambda = parEBENcv$lambda.optimal,
                                    alpha = cv$alpha.optimal,
                                    Epis = "no",
                                    verbose = 1)
```

## To Do List

- [x] Binomial prior cross-validation script with doParallel.
- [x] Gaussian prior cross-validation script with doParallel.
- [x] Binomial prior cross-validation script with doMPI.
- [x] Gaussian prior cross-validation script with doMPI.
- [x] Binomial prior cross-validation script with Microsoft ML Server (RevoScaleR/doRSR).
- [x] Gaussian prior cross-validation script with Microsoft ML Server (RevoScaleR/doRSR).
- [ ] Binomial prior cross-validation script with SparkR.
- [ ] Gaussian prior cross-validation script with SparkR.
- [ ] Binomial prior cross-validation script with CUDA.
- [ ] Gaussian prior cross-validation script with CUDA.
- [x] Manual File/Usage Instructions.

## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](LICENSE) file for details
