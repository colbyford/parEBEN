# A Parallelized Strategy for Epistasis Analysis Based on Empirical Bayesian Elastic Net Models

<h4 align = "right">Jia Wen*, Colby T. Ford*, Daniel Janies, Xinghua Shi</h4>
<h4 align = "right">Department of Bioinformatics and Genomics, University of North Carolina at Charlotte, 28223 USA</h4>

## Abstract

Epistasis reflects the distortion that the combination effect of two or more genes or variants from the sum effect of individual gene or variant on human diseases or quantitative traits. Epistasis is an important genetic foundation underlying human diseases and quantitative traits. There are two barriers in identifying epistasis using large genomic data. One is that epistasis analysis will induce the over-fitting problem of over-saturated model with high-dimensionality of genomic dataset. It therefore demands efficient statistical methods to solve this problem. 
Another one is the intensive computing time to identify epistasis among two or more variants or genetic components even though the appropriate model and data are specified. In this study, we combine statistical techniques and computing techniques that matrix multiplication and parallelization processing using empirical Bayesian Elastic Net (EBEN) method on the epistasis analysis. In our workflow, we first apply the matrix manipulation strategy that pre-computes the correlation matrix using matrix multiplication between the features and phenotype. This strategy can help narrow down the search space and thus accelerates the modeling process of epistasis analysis. Next, a parallelized version of the EBEN algorithm, parEBEN, is developed to accelerate the computing time. The parallelized package can afford multi-fold speed up in comparison with the classical EBEN method. Real yeast data based simulation demonstrates the reality of epistasis analysis of larger and more complex genomic datasets. We then analyze the real yeast data to identify some marginal and epistasis effects of genetic variants associated with quantitative traits relevant to yeast fitness.  


## Dataset

The yeast dataset is obtained from Bloom et al (2015) including 4,390 samples and 28,220 SNPs [1].


## Reference

[1] Bloom, Joshua S., et al. Genetic interactions contribute less than additive effects to quantitative trait variation in yeast. Nature communications 6 (2015).
