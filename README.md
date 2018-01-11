# parEBEN - A parallel implementation of the Empirical Bayesian Elastic Net in R

## Abstract

The Empirical Bayesian Elastic Net (EBEN) algorithm was developed by [Huang et al.](https://www.nature.com/articles/hdy201479) for the analysis of quantitative trait loci (QTLs) and gene-gene interactions (epistasis). In addition to the algorithm, the group also created the [EBEN package for R](https://cran.r-project.org/package=EBEN). This package includes functions to generate the elastic nets for both binomial and gaussian models. These functions are efficient and do not require large amounts of computational time. However, the package also includes functions for the cross validation of those models. While essential, this step is a considerably more complex task. The cross-validation functions perform a sweep to determine hyperparameters and minimize prediction error. More specifically, a five-fold cross-validation sweep is performed to minimize error by trying combinations of two parameters (λ<sub>1</sub> and λ<sub>2</sub>) in a stepped manner. Experimentally, it has been shown that this can take hours or days to perform, even on normal-sized genetics datasets.

To combat this complexity issue, the parallelization of the cross-validation functions was performed by employing parallel packages in R. By parallelizing over multiple cores of a CPU then over multiple nodes of a Hadoop cluster, a drastic time reduction was seen with no negative effect on the resultant EBEN models. By reducing the computation time, more epistasis and QTL research can be completed without such a delay. This also opens the door for larger genetic datasets to be analyzed as opposed to limiting the research due to time and computing resource constraints. Thus, parallelizing the cross-validation of the EBEN models will prove to be greatly beneficial in future genomics research.

## To Do List
- [x] Binomial prior cross-validation script with doParallel.
- [ ] Gaussian prior cross-validation script with doParallel.
- [ ] Binomial prior cross-validation script with SparkR.
- [ ] Gaussian prior cross-validation script with SparkR.
- [ ] Binomial prior cross-validation script with Microsoft ML Server (RevoScaleR).
- [ ] Gaussian prior cross-validation script with Microsoft ML Server (RevoScaleR).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
