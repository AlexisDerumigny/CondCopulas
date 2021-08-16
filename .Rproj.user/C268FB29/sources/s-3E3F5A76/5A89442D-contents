
Conditional copulas
=====================


# With pointwise conditioning


## Estimation of conditional copulas

* estimationCondCopulas.R
    * estimateNPCondCopula: nonparametric estimation of conditional copulas
    * estimateParCondCopula: parametric estimation of conditional copulas
    * estimateParCondCopula_ZIJ: parametric estimation of conditional copulas
    using (already computed) conditional pseudo-observations


## Estimation of conditional Kendall's tau (CKT)

A general wrapper function:

* CKT.estimate: that can be used for any method of estimating conditional Kendall's tau

### Kernel-based estimation of conditional Kendall's tau

* CKT.kernel: for any number of variable and with possible choice of the bandwidth

### Kendall's regression

* CKT.kendallReg.fit: fit Kendall's regression, a regression-like method for the estimation of conditional Kendall's tau
* CKT.kendallReg.predict: for prediction of the new conditional Kendall's tau (given new covariates)


### Classification-based estimation of conditional Kendall's tau

* using tree:
    * CKT.fit.tree: for fitting the tree
    * CKT.predict.tree: for prediction of the tree
    
* using random forests:
    * CKT.fit.randomForest: for fitting the random forest
    * CKT.predict.randomForest: for prediction of new conditional Kendall's taus    

* using nearest neighbors:
    * CKT.predict.kNN: for several numbers of nearest neighbors
    
* using neural networks:
    * CKT.fit.nNets: for fitting the neural networks
    * CKT.predict.nNets: for prediction of new conditional Kendall's taus
    
* using GLM:
    * CKT.fit.GLM: for fitting the GLM
    * CKT.predict.GLM: for prediction of the GLM

### Advanced functions for manual hyperparameter choices

* CKT.hCV.Kfolds: for K-fold cross-validation choice of the bandwidth for kernel smoothing

* CKT.hCV.l1out: for leave-one-out cross-validation choice of the bandwidth for kernel smoothing

* CKT.KendallReg.LambdaCV : cross-validated choice of the penalization parameter lambda

* CKT.adaptkNN: for a (local) aggregation of the number of nearest neighbors based on Lepski's method


