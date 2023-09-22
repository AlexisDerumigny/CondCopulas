
# How to install

The release version on CRAN:
```r
install.packages("CondCopulas")
```

The development version from GitHub, using the `devtools` package:
``` r
# install.packages("devtools")
devtools::install_github("AlexisDerumigny/CondCopulas")
```

If you have any questions or suggestions, feel free to
[open an issue](https://github.com/AlexisDerumigny/CondCopulas/issues/new).


Conditional copulas with pointwise conditioning
=================================================

In this first part, we are interesting in the inference of the conditional copula
of a random vector $X$ given the pointwise conditioning $Z = z$,
where $Z$ is another random vector and $z$ is a fixed value.

## Tests of the simplifying assumption

These functions perform a test of the "simplifying assumption"
that the conditional copula $C_{X | Z = z}$ does not depend on the value of $z$.

* `simpA.NP`: in a purely nonparametric framework

* `simpA.param`: assuming that the conditional copula
belongs to a parametric family of copulas for all values of the conditioning variable

* `simpA.kendallReg`: test of the simplifying assumption based on the constancy
of the conditional Kendall's tau assuming that it satisfies a regression-like equation


## Estimation of conditional copulas (using kernel smoothing)

These functions estimate the conditional copula $C_{X | Z = z}$
in different frameworks.

* `estimateNPCondCopula`: nonparametric estimation of conditional copulas.

* `estimateParCondCopula`: parametric estimation of conditional copulas.

* `estimateParCondCopula_ZIJ`: parametric estimation of conditional copulas
using (already computed) conditional pseudo-observations.


## Estimation of conditional Kendall's tau (CKT)

In this part, we assume that the dimension of $X$ is $2$, i.e. $X = (X_1, X_2)$.
Instead of estimating the conditional copula $C_{X | Z = z}$,
which is an infinite-dimensional object,
it is possible to estimate the conditional Kendall's tau (CKT) $\tau_{1,2|Z=z}$,
which is a real number for every given value of $z$.


To estimate the conditional Kendall's tau,
the package provides a general wrapper function:

* `CKT.estimate`: that can be used for any method of estimating conditional Kendall's tau.
Each of these methods is detailed below and has its own function.


### Kernel-based estimation of conditional Kendall's tau

* `CKT.kernel`: use kernel smoothing to estimate the conditional Kendall's tau.
The bandwidth can be given by the user or determined by cross-validation.


### Kendall's regression

* `CKT.kendallReg.fit`: fit Kendall's regression, a regression-like method for the estimation of conditional Kendall's tau.

* `CKT.kendallReg.predict`: predict the conditional Kendall's tau
given new values $z$ of the covariates.


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



Conditional copulas with discrete conditioning by Borel sets
==============================================================

In this second part, we are interesting in the inference of the conditional copula
of a random vector $X$ given the discrete conditioning $Z \in A$,
where $Z$ is another random vector and $A$ is a Borel subset of possible values of $A$.


## Test of the hypothesis that the conditioning Borel subset has no influence on the conditional copula

These functions perform a test of the hypothesis
that the conditional copula $C_{X | Z \in A}$ does not depend on the value of $A$
for different choices of the conditioning set $A$.

* `bCond.simpA.param` : test of this hypothesis, assuming that the copula belongs to a parametric family

* `bCond.simpA.CKT`: test of the hypothesis that conditional Kendall's tau are equal
over all the different conditioning subsets.


## Estimation

* `bCond.pobs` : computation of the conditional pseudo-observations
$F_{1|A(i)}(X_{i,1} | A(i))$ and $F_{2|A(i)}(X_{i,2} | A(i))$ for every $i=1, \dots, n$.

* `bCond.estParamCopula` : estimation of a conditional parametric copula,
i.e. for every set $A$, a conditional parameter $\theta(A)$ is estimated.


## Data-driven choice of conditioning subsets

* `bCond.treeCKT`: construction of binary tree whose leaves corresponds to the most relevant conditioning subsets
(in the sense of maximizing the difference between estimated conditional Kendall's taus).


# References

Derumigny, A., & Fermanian, J. D. (2017).
About tests of the “simplifying” assumption for conditional copulas.
*Dependence Modeling*, 5(1), 154-197.
[pdf](https://doi.org/10.1515/demo-2017-0011)

Derumigny, A., & Fermanian, J. D. (2019).
A classification point-of-view about conditional Kendall’s tau.
*Computational Statistics & Data Analysis*, 135, 70-94.
[pdf](https://doi.org/10.1016/j.csda.2019.01.013)

Derumigny, A., & Fermanian, J. D. (2019).
On kernel-based estimation of conditional Kendall’s tau:
finite-distance bounds and asymptotic behavior.
*Dependence Modeling*, 7(1), 292-321.
[pdf](https://doi.org/10.1515/demo-2019-0016)

Derumigny, A., & Fermanian, J. D. (2020).
On Kendall’s regression.
*Journal of Multivariate Analysis*, 178, 104610.
[pdf](https://doi.org/10.1016/j.jmva.2020.104610)

Derumigny, A., & Fermanian, J. D. (2022).
Conditional empirical copula processes and generalized dependence measures.
*Electronic Journal of Statistics*, 16(2), 5692-5719.
[pdf](https://doi.org/10.1214/22-EJS2075)

Derumigny, A., Fermanian, J. D., & Min, A. (2022).
Testing for equality between conditional copulas
given discretized conditioning events.
*Canadian Journal of Statistics*.
[pdf](https://doi.org/10.1002/cjs.11742)

van der Spek, R., & Derumigny, A. (2022). Fast estimation of Kendall’s
Tau and conditional Kendall’s Tau matrices under structural assumptions.
[arXiv:2204.03285](https://arxiv.org/pdf/2204.03285.pdf).

