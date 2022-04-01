
# How to install

The release version on CRAN:
```r
install.packages("CondCopulas")
```

From GitHub, using the `devtools` package:
``` r
# install.packages("devtools")
devtools::install_github("AlexisDerumigny/CondCopulas")
```


Conditional copulas
=====================


# With pointwise conditioning

## Tests of the simplifying assumption

* `simpA.NP`: in a purely nonparametric framework

* `simpA.param`: assuming that the conditional copula
belongs to a parametric family of copulas for all values of the conditioning variable

* `simpA.kendallReg`: test of the simplifying assumption based on the constancy
of the conditional Kendall's tau assuming that it satisfies a regression-like equation


## Estimation of conditional copulas (using kernel smoothing)


* `estimateNPCondCopula`: nonparametric estimation of conditional copulas

* `estimateParCondCopula`: parametric estimation of conditional copulas

* `estimateParCondCopula_ZIJ`: parametric estimation of conditional copulas
using (already computed) conditional pseudo-observations


## Estimation of conditional Kendall's tau (CKT)

A general wrapper function:

* `CKT.estimate`: that can be used for any method of estimating conditional Kendall's tau.
Each of these methods is detailed below and has its own function.

### Kernel-based estimation of conditional Kendall's tau

* `CKT.kernel`: for any number of variable and with possible choice of the bandwidth

### Kendall's regression

* `CKT.kendallReg.fit`: fit Kendall's regression, a regression-like method for the estimation of conditional Kendall's tau

* `CKT.kendallReg.predict`: for prediction of the new conditional Kendall's tau (given new covariates)


### Classification-based estimation of conditional Kendall's tau

* using tree:
    * `CKT.fit.tree`: for fitting a tree-based model for the conditional Kendall's tau
    * `CKT.predict.tree`: for prediction of new conditional Kendall's taus    
    
* using random forests:
    * `CKT.fit.randomForest`: for fitting a random forest-based model for the conditional Kendall's tau
    * `CKT.predict.randomForest`: for prediction of new conditional Kendall's taus    

* using nearest neighbors:
    * `CKT.predict.kNN`: for several numbers of nearest neighbors
    
* using neural networks:
    * `CKT.fit.nNets`: for fitting a neural networks-based model for the conditional Kendall's tau
    * `CKT.predict.nNets`: for prediction of new conditional Kendall's taus
    
* using GLM:
    * `CKT.fit.GLM`: for fitting a GLM-like model for the conditional Kendall's tau
    * `CKT.predict.GLM`: for prediction of new conditional Kendall's taus


### Advanced functions for manual hyperparameter choices

* `CKT.hCV.Kfolds`: for K-fold cross-validation choice of the bandwidth for kernel smoothing

* `CKT.hCV.l1out`: for leave-one-out cross-validation choice of the bandwidth for kernel smoothing

* `CKT.KendallReg.LambdaCV` : cross-validated choice of the penalization parameter lambda

* `CKT.adaptkNN`: for a (local) aggregation of the number of nearest neighbors based on Lepski's method


# With discrete conditioning by Borel sets

## Test of the assumption that the conditioning Borel subset has no influence on the conditional copula

* `bCond.simpA.param` : assuming that the copula belongs to a parametric family


## Estimation

* `bCond.pobs` : computation of the conditional pseudo-observations
$F_{1|A(i)}(X_{i,1} | A(i))$ and $F_{2|A(i)}(X_{i,2} | A(i))$ for every $i=1, \dots, n$.

* `bCond.estParamCopula` : estimation of a conditional parametric copula,
i.e. for every set $A$, a conditional parameter $\theta(A)$ is estimated.


# References

Derumigny, A., & Fermanian, J. D. (2017).
About tests of the “simplifying” assumption for conditional copulas.
*Dependence Modeling*, 5(1), 154-197.

Derumigny, A., & Fermanian, J. D. (2019).
A classification point-of-view about conditional Kendall’s tau.
*Computational Statistics & Data Analysis*, 135, 70-94.

Derumigny, A., & Fermanian, J. D. (2019).
On kernel-based estimation of conditional Kendall’s tau:
finite-distance bounds and asymptotic behavior.
*Dependence Modeling*, 7(1), 292-321.

Derumigny, A., & Fermanian, J. D. (2020).
On Kendall’s regression.
*Journal of Multivariate Analysis*, 178, 104610.

