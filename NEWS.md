
* Function `simpA.kendallReg` can handle the case where only one regressor is given.
It also uses `stats::lm.fit` for the unpenalized regression.
A typo in the Wald test statistic has been fixed.
Its output is an `S3` object of class `simpA_kendallReg_test` with `print` and `plot` methods.

* Function `CKT.kernel` has now more options to control the possible display of
the progress bar to show the progress of the computation.

* New dependency: `testthat` has been added to `Suggests`.


# CondCopulas 0.1.3

* Adding a warning to `CKT.kernel()` when some estimated conditional Kendall's
taus are `NA` because of a too small bandwidth.

* Fixing a bug for `CKT.kernel()` when the conditioning variable is multivariate.

* Adding and updating references for conditional copulas with discretized conditioning events.

* Fix an error when running `bCond.simpA.CKT()`.

* Fix default value of the argument `minSize` in `bCond.treeCKT()` to be
`minSize = minProb * nrow(XI)` as intended.

* Functions `CKT.kernel()` and `CKT.estimate()` now warn and return `numeric(0)`
when the argument `newZ` is `numeric(0)`.


# CondCopulas 0.1.2

* New dependence `wdm` instead of `pcaPP` for fast computation of Kendall's tau.

* Functions `CKT.kernel`, `CKT.hCV.l1out`, `CKT.hCV.Kfolds` and `CKTmatrix.kernel`
gain a new choice `"wdm"` for the argument `typeEstCKT`. This new choice allows
for faster computation using the package `wdm`. For observations without ties,
it is equivalent to the previous choice `typeEstCKT = 4`.


# CondCopulas 0.1.1

* Adding DOIs in the DESCRIPTION and small documentation fixes.


# CondCopulas 0.1.0

* Initial release
