

# Weighting ------------------------------------------------------------------


#' Computes kernel weights (univariate)
#'
#' @param vectorZ vector of observed data
#' @param pointZ point at which the weights should be computed
#' @param h bandwidth
#' @param kernel.name name of the kernel. Possible choices are
#' "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#' @param normalization if TRUE, normalize by the sum of the weights
#'
#' @return a vector of the same size as vectorZ containing the weights
#' for each point.
#'
#' @examples
#' vectorZ = seq(0,1, by = 0.1)
#' pointZ = 0.3
#' h = 0.2
#' my_weights = computeWeights.univariate(
#'   vectorZ = vectorZ, h = h, pointZ = pointZ,
#'   kernel.name = "Gaussian")
#'
#' @noRd
#'
computeWeights.univariate <- function(vectorZ, h, pointZ,
                                      kernel.name, normalization = TRUE)
{
  u = (vectorZ - pointZ) / h

  switch (
    kernel.name,

    "Gaussian" = {listWeights = exp(- u^2)},

    "Epa" = {listWeights = 0.75 * (1 - u^2) * as.numeric(abs(u) <= 1)},

    {stop(paste0("kernel.name: ", kernel.name,
                 " does not belong to the list of possible names for the kernels." ))}
  )

  if (normalization) {
    listWeights = listWeights/sum(listWeights)
  }

  return (listWeights)
}


#' Computes kernel weights (multivariate)
#'
#' @param matrixZ matrix of observed data of dimension n*p
#' @param pointZ point of dimension p at which the weights should be computed
#' @param h bandwidth
#' @param kernel.name name of the kernel. Possible choices are
#' "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel)
#' @param normalization if TRUE, normalize by the sum of the weights
#'
#' @return a vector of length n containing the weights
#' for each point.
#'
#' @examples
#' n = 20
#' matrixZ = cbind(rnorm(20), rnorm(20))
#' h = 0.2
#' pointZ = c(2.1, 3.2)
#' my_weights = computeWeights.multivariate(
#'   matrixZ = matrixZ, h = 0.2, pointZ = pointZ,
#'   kernel.name = "Gaussian")
#'
#' n = 20
#' matrixZ = matrix(rnorm(20), ncol = 1)
#' h = 0.2
#' pointZ = c(2.1)
#' my_weights = computeWeights.multivariate(
#'   matrixZ = matrixZ, h = 0.2, pointZ = pointZ,
#'   kernel.name = "Gaussian")
#' my_weights_ = computeWeights.univariate(
#'   vectorZ = as.numeric(matrixZ), h = 0.2, pointZ = pointZ,
#'   kernel.name = "Gaussian")
#' my_weights == my_weights_
#'
#' @noRd
#'
computeWeights.multivariate <- function(matrixZ, h, pointZ,
                                        kernel.name, normalization = TRUE)
{
  u = sweep(matrixZ, MARGIN = 2, STATS = pointZ) / h

  switch (
    kernel.name,

    "Gaussian" = {listWeights = apply(X = exp(- u^2),
                                      MARGIN = 1, FUN = prod)},

    "Epa" = {listWeights = apply(X = 0.75 * (1 - u^2) * (abs(u) <= 1),
                                 MARGIN = 1, FUN = prod)},

    {stop(paste0("kernel.name: ", kernel.name,
                 " does not belong to the list of possible names for the kernels." ))}
  )

  if (normalization) {
    listWeights = listWeights/sum(listWeights)
  }
  return (listWeights)
}



# Pointwise estimation of CKT ----------------------------------------------------


#' Estimate the conditional Kendall's tau of X1 and X2
#' at a fixed univariate point Z = pointZ
#'
#' @param X1,X2 vectors of observations of the conditioned variables
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#' @param vectorZ vector of observed points of Z.
#' It shall have the same length as the number of rows of matrixSignsPairs.
#' @param h bandwidth
#' @param pointZ point at which the conditional Kendall's tau is computed.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#'
#' @return an estimator of the conditional Kendall's tau
#' of X1 and X2 given Z = z.
#'
#' @noRd
#'
CKT.kernelPointwise.univariate <- function(X1, X2, matrixSignsPairs, vectorZ,
                                           h, pointZ, kernel.name, typeEstCKT)
{
  listWeights = computeWeights.univariate(vectorZ, h, pointZ, kernel.name)

  if(typeEstCKT == "wdm"){

    estimate = wdm::wdm(x = X1, y = X2, weights = listWeights, method = "kendall")

  } else {
    matrixWeights = outer(listWeights, listWeights)

    switch (
      typeEstCKT,
      # 1
      { estimate = 4 * sum(matrixWeights * matrixSignsPairs) - 1 } ,
      # 2
      { estimate = sum(matrixWeights * matrixSignsPairs) },
      # 3
      { estimate = 1 - 4 * sum(matrixWeights * matrixSignsPairs) },
      # 4
      { estimate = sum(matrixWeights * matrixSignsPairs) / (1 - sum(listWeights^2)) },
      {stop(paste0("typeEstCKT: ", typeEstCKT, " is not in {1,2,3,4}" ) ) }
    )
  }

  return(estimate)
}



#' Estimate the conditional Kendall's tau of X1 and X2
#' at a fixed multivariate point Z = pointZ
#'
#' @param X1,X2 vectors of observations of the conditioned variables
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by \code{computeMatrixSignPairs}.
#' @param vectorZ vector of observed points of Z.
#' It shall have the same length as the number of rows of \code{matrixSignsPairs}.
#' @param h bandwidth
#' @param pointZ point at which the conditional Kendall's tau is computed.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#'
#' @return an estimator of the conditional Kendall's tau
#' of \eqn{X_1} and \eqn{X_2} given \eqn{Z = z}.
#'
#' @noRd
#'
CKT.kernelPointwise.multivariate <- function(X1, X2, matrixSignsPairs, matrixZ,
                                             h, pointZ, kernel.name, typeEstCKT)
{
  if (kernel.name == "Epa"){
    # For faster computation, only uses points with non-zero values for the kernel
    u = sweep(matrixZ, MARGIN = 2, STATS = pointZ)

    isSmaller_h = apply(X = u, MARGIN = 1, FUN = function(x){
      return (all(abs(x) <= h))
    })

    whichNonZero = which( isSmaller_h )

    # We need at least 2 points, i.e. at least 1 pair,
    # to estimate (conditional) Kendall's tau
    if (length(whichNonZero) <= 1){
      return (NA)
    } # now `length(whichNonZero)` is at least 2.

    matrixZ = matrixZ[whichNonZero, ]
    matrixSignsPairs = matrixSignsPairs[whichNonZero, whichNonZero]
    X1 = X1[whichNonZero]
    X2 = X2[whichNonZero]
  }

  listWeights = computeWeights.multivariate(matrixZ, h, pointZ, kernel.name)

  if(typeEstCKT == "wdm"){

    estimate = wdm::wdm(x = X1, y = X2, weights = listWeights, method = "kendall")

  } else {
    matrixWeights = outer(listWeights, listWeights)

    switch (
      typeEstCKT,
      # 1
      { estimate = 4 * sum(matrixWeights * matrixSignsPairs) - 1 } ,
      # 2
      { estimate = sum(matrixWeights * matrixSignsPairs) },
      # 3
      { estimate = 1 - 4 * sum(matrixWeights * matrixSignsPairs) },
      # 4
      { estimate = sum(matrixWeights * matrixSignsPairs) / (1 - sum(listWeights^2)) },

      {stop(paste0("typeEstCKT: ", typeEstCKT, " is not in {1,2,3,4}" ) ) }
    )
  }
  # if(! is.finite(estimate) ) {print(whichNonZero) ; print(listWeights) ; print(pointZ)}
  return(estimate)
}



# Estimation of CKT at multiple points ----------------------------------------------


#' Estimate the conditional Kendall's tau of X1 and X2
#' at different points
#'
#'
#' @param X1,X2 vectors of observations of the conditioned variables
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#'
#' @param Z vector of observed points of Z.
#' It shall have the same length as the number of rows of \code{matrixSignsPairs}.
#'
#' @param h bandwidth. It can be a real, in this case the same \code{h}
#' will be used for every element of \code{vectorZToEstimate}.
#' If \code{h}is a vector then its elements are recycled to match the length of
#' \code{vectorZToEstimate}.
#'
#' @param ZToEstimate points at which the conditional Kendall's tau is computed.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#' @param progressBar if TRUE, a progressbar is displayed to show the computation.
#'
#' @return a vector of the same length as vectorZToEstimate whose elements
#' are the estimated conditional Kendall's taus of X1 and X2 given Z = z.
#'
#' @noRd
#'
CKT.kernel.univariate <- function(X1, X2, matrixSignsPairs, Z,
                                  h, ZToEstimate,
                                  kernel.name = "Epa", typeEstCKT = 4,
                                  progressBar = TRUE)
{
  if (typeEstCKT != "wdm"){
    .check_MatrixSignPairs(matrixSignsPairs)

    if (nrow(matrixSignsPairs) != length(Z)){
      stop(paste0("Z must have the same length ",
                  "as the number of rows of matrixSignsPairs."))
    }
  }

  n_prime = length(ZToEstimate)
  if (length(h) == 1) {
    h_vect = rep(h, n_prime)
  } else {
    h_vect = h
  }

  if (progressBar) {
    estimates = pbapply::pbapply(
      X = array(1:n_prime), MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.univariate(
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i],
        matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], vectorZ = Z,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  } else {
    estimates = apply(
      X = array(1:n_prime), MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.univariate(
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i],
        matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], vectorZ = Z,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  }

  return(estimates)
}


#' Estimate the conditional Kendall's tau of X1 and X2
#' at different points
#'
#'
#' @param X1,X2 vectors of observations of the conditioned variables
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#'
#' @param Z matrix of observed points of Z.
#' It shall have the number of rows as \code{matrixSignsPairs}.
#'
#' @param h bandwidth. It can be a real, in this case the same \code{h}
#' will be used for every element of \code{vectorZToEstimate}.
#' If \code{h}is a vector then its elements are recycled to match the length of
#' \code{vectorZToEstimate}.
#'
#' @param ZToEstimate points at which the conditional Kendall's tau is computed.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#' @param progressBar if TRUE, a progressbar is displayed to
#' show the progress of the computation.
#'
#' @return a vector of the same length as vectorZToEstimate whose elements
#' are the estimated conditional Kendall's taus of X1 and X2 given Z = z.
#'
#' @noRd
#'
CKT.kernel.multivariate <- function(X1, X2, matrixSignsPairs, Z,
                                    h, ZToEstimate,
                                    kernel.name = "Epa", typeEstCKT = 4,
                                    progressBar = TRUE)
{
  if (typeEstCKT != "wdm"){
    .check_MatrixSignPairs(matrixSignsPairs)

    if (NROW(matrixSignsPairs) != NROW(Z)){
      stop("'Z' and 'matrixSignsPairs' must have the same number of rows.")
    }
  }

  .checkSame_ncols_Z_newZ(Z, ZToEstimate, name_Z = "Z", name_newZ = "ZToEstimate")

  dim_Z = ncol(Z)
  n_prime = nrow(ZToEstimate)

  if (length(h) == 1) {
    h_vect = rep(h, n_prime)
  } else {
    h_vect = h
  }

  if (progressBar) {
    estimates = pbapply::pbapply(
      X = array(1:n_prime), MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.multivariate(
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i,],
        matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], matrixZ = Z,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  } else {
    estimates = apply(
      X = 1:n_prime, MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.multivariate(
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i,],
        matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], matrixZ = Z,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  }

  return(estimates)
}



#' Estimation of conditional Kendall's tau using kernel smoothing
#'
#'
#' Let \eqn{X_1} and \eqn{X_2} be two random variables.
#' The goal of this function is to estimate the conditional Kendall's tau
#' (a dependence measure) between \eqn{X_1} and \eqn{X_2} given \eqn{Z=z}
#' for a conditioning variable \eqn{Z}.
#' Conditional Kendall's tau between \eqn{X_1} and \eqn{X_2} given \eqn{Z=z}
#' is defined as:
#' \deqn{P( (X_{1,1} - X_{2,1})(X_{1,2} - X_{2,2}) > 0 | Z_1 = Z_2 = z)}
#' \deqn{- P( (X_{1,1} - X_{2,1})(X_{1,2} - X_{2,2}) < 0 | Z_1 = Z_2 = z),}
#' where \eqn{(X_{1,1}, X_{1,2}, Z_1)} and \eqn{(X_{2,1}, X_{2,2}, Z_2)}
#' are two independent and identically distributed copies of \eqn{(X_1, X_2, Z)}.
#' For this, a kernel-based estimator is used, as described in
#' (Derumigny, & Fermanian (2019)).
#'
#'
#' \strong{Choice of the bandwidth \code{h}}.
#' The choice of the bandwidth must be done carefully.
#' In the univariate case, the default kernel (Epanechnikov kernel) has a support
#' on \eqn{[-1,1]}, so for a bandwidth \code{h}, estimation of conditional Kendall's
#' tau at \eqn{Z=z} will only use points for which \eqn{Z_i \in [z \pm h]}.
#' As usual in nonparametric estimation, \code{h} should not be too small
#' (to avoid having a too large variance) and should not be large
#' (to avoid having a too large bias).
#'
#' We recommend that for each \eqn{z} for which the conditional Kendall's tau
#' \eqn{\tau_{X_1, X_2 | Z=z}} is estimated, the set
#' \eqn{\{i: Z_i \in [z \pm h] \}}
#' should contain at least 20 points and not more than 30\% of the points of
#' the whole dataset.
#' Note that for a consistent estimation, as the sample size \eqn{n} tends
#' to the infinity, \code{h} should tend to \eqn{0} while the size of the set
#' \eqn{\{i: Z_i \in [z \pm h]\}} should also tend to the infinity.
#' Indeed the conditioning points should be closer and closer to the point of
#' interest \eqn{z} (small \code{h}) and more and more numerous
#' (\code{h} tending to 0 slowly enough).
#'
#' In the multivariate case, similar recommendations can be made.
#' Because of the curse of dimensionality, a larger sample will be necessary to
#' reach the same level of precision as in the univariate case.
#'
#'
#' @param X1 a vector of n observations of the first variable
#' (or a 1-column matrix)
#'
#' @param X2 a vector of n observations of the second variable
#' (or a 1-column matrix)
#'
#' @param Z a vector of n observations of the conditioning variable,
#' or a matrix with n rows of observations of the conditioning vector
#' (in the case that several conditioning variables are given; in this case,
#' each column corresponds to 1 conditioning variable). It can also be a
#' \code{data.frame}, provided that all the entries are \code{numeric}.
#'
#' @param newZ the new data of observations of Z at which the conditional
#' Kendall's tau should be estimated. It must have the same number of column
#' as \code{Z}. It can be a vector (if \code{NCOL(Z) == 1}), a \code{matrix} or
#' a \code{data.frame}, provided that all the entries are \code{numeric}.
#'
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' Possible choices are \itemize{
#'   \item \code{1} and \code{3} produced biased estimators.
#'   \code{2} does not attain the full range \eqn{[-1,1]}.
#'   Therefore these 3 choices are not recommended for applications on real data.
#'
#'   \item \code{4} is an improved version of \code{1,2,3} that has less bias
#'   and attains the full range \eqn{[-1,1]}.
#'
#'   \item \code{"wdm"} is the default version and produces the same results
#'   as \code{4} when they are no ties in the data.
#' }
#'
#' @param methodCV method used for the cross-validation.
#' Possible choices are \code{"leave-one-out"} and \code{"Kfolds"}.
#'
#' @param nPairs number of pairs used in the cross-validation criteria,
#' if \code{methodCV = "leave-one-out"}. Use \code{nPairs = "all"} to choose
#' all pairs. The default is \code{nPairs = 10 * n}, where \code{n} is the
#' sample size.
#'
#' @param Kfolds number of subsamples used,
#' if \code{methodCV = "Kfolds"}.
#'
#' @param h the bandwidth used for kernel smoothing.
#' If this is a vector, then cross-validation is used following the method
#' given by argument \code{methodCV} to choose the best bandwidth
#' before doing the estimation.
#'
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are \code{"Gaussian"} (Gaussian kernel)
#' and \code{"Epa"} (Epanechnikov kernel).
#'
#' @param progressBar control the display of progress bars.
#' Possible choices are: \itemize{
#'   \item \code{0} no progress bar is displayed
#'
#'   \item \code{1} a general progress bar is displayed
#'
#'   \item \code{2} and larger values:
#'   a general progress bar is displayed, and additionally,
#'   a progressbar for each value of \code{h} is displayed
#'   to show the progress of the computation.
#'   This only applies when the bandwidth is chosen by cross-validation
#'   (i.e. when \code{h} is a vector).
#' }
#'
#' @param warnNA a Boolean to indicate whether warnings should be raised if
#' \code{NA}s are produced. By default it is \code{TRUE}. If \code{warnNA=FALSE},
#' then no warning is raised even if \code{NA}s are produced. This is the case
#' usually if either the bandwidth \code{h} is too small, or if there are already
#' \code{NA}s in (some of) the inputs.
#'
#' @param observedX1,observedX2,observedZ old parameter names for \code{X1},
#' \code{X2}, \code{Z}. Support for this will be removed at a later version.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' On kernel-based estimation of conditional Kendall’s tau:
#' finite-distance bounds and asymptotic behavior.
#' Dependence Modeling, 7(1), 292-321.
#' \doi{10.1515/demo-2019-0016}
#'
#' @return a list with components:
#' \itemize{
#'    \item \code{estimatedCKT} the vector of size \code{NROW(newZ)}
#'    containing the values of the estimated conditional Kendall's tau.
#'
#'    \item \code{finalh} the bandwidth \code{h} that was finally used
#'    for kernel smoothing (either the one specified by the user
#'    or the one chosen by cross-validation if multiple bandwidths were given.)
#'
#'    \item \code{resultCV} (only in case of cross-validation). This gives the
#'    output of the cross-validation function that is used, i.e. the output of
#'    either \code{\link{CKT.hCV.l1out}} or \code{\link{CKT.hCV.Kfolds}}.
#' }
#'
#' @seealso \code{\link{CKT.estimate}} for other estimators
#' of conditional Kendall's tau.
#' \code{\link{CKTmatrix.kernel}} for a generalization of this function
#' when the conditioned vector is of dimension \code{d}
#' instead of dimension \code{2} here.
#'
#' See \code{\link{CKT.hCV.l1out}} for manual selection of the bandwidth \code{h}
#' by leave-one-out or K-folds cross-validation.
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 100
#' # This is a small example for performance reason.
#' # For a better example, use:
#' # N = 800
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' newZ = seq(2,10,by = 0.1)
#' estimatedCKT_kernel <- CKT.kernel(
#'    X1 = X1, X2 = X2, Z = Z,
#'    newZ = newZ, h = 0.1, kernel.name = "Epa")$estimatedCKT
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in red)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col = "black",
#'      type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_kernel, col = "red")
#'
#' # Multivariate example
#' N = 100
#' # This is a small example for performance reason.
#' # For a better example, use:
#' # N = 1000
#' Z1 = rnorm(n = N, mean = 5, sd = 2)
#' Z2 = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z1 - Z2, mean = 2, sd = 2)
#' simCopula = VineCopula::BiCopSim(N = N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' Z = cbind(Z1, Z2)
#'
#' newZ = expand.grid(Z1 = seq(2,8,by = 0.5),
#'                    Z2 = seq(2,8,by = 1))
#' estimatedCKT_kernel <- CKT.kernel(
#'    X1 = X1, X2 = X2, Z = Z,
#'    newZ = newZ, h = 1, kernel.name = "Epa")$estimatedCKT
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   df = rbind(
#'     data.frame(newZ, CKT = estimatedCKT_kernel,
#'                type = "estimated CKT") ,
#'     data.frame(newZ, CKT = -0.9 + 1.8 * pnorm(newZ$Z1 - newZ$Z2,
#'                                               mean = 2, sd = 2),
#'                type = "true CKT")
#'   )
#'
#'   ggplot2::ggplot(df) +
#'   ggplot2::geom_tile(ggplot2::aes(x = Z1, y = Z2, fill = CKT)) +
#'   ggplot2::facet_grid(as.formula("~type"))
#' }
#'
#' @export
#'
CKT.kernel <- function(X1 = NULL, X2 = NULL, Z = NULL, newZ,
                       h, kernel.name = "Epa",
                       se = FALSE, confint = FALSE, level = 0.95,
                       methodCV = "Kfolds",
                       Kfolds = 5, nPairs = NULL,
                       typeEstCKT = "wdm", progressBar = 1,
                       warnNA = TRUE,
                       observedX1 = NULL, observedX2 = NULL, observedZ = NULL)
{
  if (length(newZ) == 0){
    warning("'newZ' is of length 0, therefore, no estimation is done.")
    return (numeric(0))
  }

  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  env = environment()
  .observedX1X2_to_X1X2(env)
  .observedZ_to_Z(env)

  # Checking for same number of observations
  .checkSame_nobs_X1X2Z(X1, X2, Z)

  # Checking that X1 and X2 are univariate
  .checkUnivX1X2(X1, X2)
  X1 = as.numeric(X1)
  X2 = as.numeric(X2)

  # Checking that the number of columns of Z and of newZ are the same
  .checkSame_ncols_Z_newZ(Z, newZ, name_Z = "Z", name_newZ = "newZ")

  # Checking the class of Z, newZ, and converting them to the right class
  Z = .ensure_Z_numeric_vector_or_matrix(Z = Z, nameZ = "Z")
  newZ = .ensure_Z_numeric_vector_or_matrix(Z = newZ, nameZ = "newZ")


  if (typeEstCKT == "wdm") {
    matrixSignsPairs = NULL
  } else {
    matrixSignsPairs = computeMatrixSignPairs(
      vectorX1 = X1, vectorX2 = X2, typeEstCKT = typeEstCKT)
  }

  if (length(h) == 1){
    finalh = h
    resultCV = NULL

  } else {

    # Do the cross-validation

    switch (
      methodCV,

      "Kfolds" = {
        resultCV = CKT.hCV.Kfolds(
          X1 = X1, X2 = X2, Z = Z,
          range_h = h, matrixSignsPairs = matrixSignsPairs,
          ZToEstimate = newZ,
          typeEstCKT = typeEstCKT, kernel.name = kernel.name, Kfolds = Kfolds,
          progressBar = progressBar > 1)
      },

      "leave-one-out" = {
        resultCV = CKT.hCV.l1out(
          X1 = X1, X2 = X2, Z = Z,
          range_h = h, matrixSignsPairs = matrixSignsPairs,
          typeEstCKT = typeEstCKT, kernel.name = kernel.name,
          nPairs = nPairs, progressBar = progressBar > 1)
      }

    )
    finalh = resultCV$hCV
  }

  # We finally computed the estimated conditional Kendall's tau
  # using the selected value for h and returns it

  if (is.vector(Z)){
    estCKT = CKT.kernel.univariate(
      X1 = X1, X2 = X2, matrixSignsPairs = matrixSignsPairs, Z = Z,
      h = finalh, ZToEstimate = newZ,
      kernel.name = kernel.name, typeEstCKT = typeEstCKT,
      progressBar = progressBar > 0)


  } else {

    estCKT = CKT.kernel.multivariate(
      X1 = X1, X2 = X2, matrixSignsPairs = matrixSignsPairs, Z = Z,
      h = finalh, ZToEstimate = newZ,
      kernel.name = kernel.name, typeEstCKT = typeEstCKT,
      progressBar = progressBar > 0)
  }

  result = list(estimatedCKT = estCKT, h = finalh,
                resultCV = resultCV,
                matrixSignsPairs = matrixSignsPairs,
                X1 = X1, X2 = X2, Z = Z, newZ = newZ,
                kernel.name = kernel.name,
                typeEstCKT = typeEstCKT)

  class(result) <- "estimated_CKT_kernel"

  # Adding additional components to result as requested
  if (confint){
    # Confidence intervals requires the knowledge of the standard error
    se <- TRUE
  }
  if (se){
    result$se <- se(result, progressBar = progressBar)
  }
  if (confint){
    result$confint = confint(object = result, level = level,
                             progressBar = progressBar)
  }

  # Raise warnings for NA values, if any  ======================================

  if (warnNA){
    if (anyNA(estCKT)){
      warnNA_CKT.kernel(
        X1 = X1, X2 = X2, Z = Z, newZ = newZ,
        estimator = estCKT,
        nameEstimator = "estimated conditional Kendall's tau")

    } else if (anyNA(result$se)){
      warnNA_CKT.kernel(
        X1 = X1, X2 = X2, Z = Z, newZ = newZ,
        estimator = result$se,
        nameEstimator = "standard error of conditional Kendall's tau")
    }
  }

  return (result)
}


warnNA_CKT.kernel <- function(X1, X2, Z, newZ, estimator, nameEstimator){
  n_NA = length(which(is.na(estimator)))

  if (!anyNA(X1) && !anyNA(X2) && !anyNA(Z) && !anyNA(newZ)){

    message = paste0(
      "NA in ", nameEstimator ," (", n_NA, " out of ", NROW(newZ), ").\n",
      "This often happens when the bandwidth h is too small, ",
      "consider using a bigger bandwidth ",
      "(see the documentation for advice on the choice of h).\n",
      "You can disable this warning using the input `warnNA = FALSE`.")

  } else {
    n_NA_X1 = length(which(is.na(X1)))
    n_NA_X2 = length(which(is.na(X2)))
    n_NA_Z = length(which(is.na(Z)))
    n_NA_newZ = length(which(is.na(newZ)))

    message = paste0(
      "NA in ", nameEstimator ," (", n_NA, " out of ", NROW(newZ), "). \n",
      "Here there are also missing values in the following inputs: \n",
      "* X1: "  , n_NA_X1  , " missing out of ", length(X1)  , "\n",
      "* X2: "  , n_NA_X2  , " missing out of ", length(X2)  , "\n",
      "* Z: "   , n_NA_Z   , " missing out of ", length(Z)   , "\n",
      "* newZ: ", n_NA_newZ, " missing out of ", length(newZ), "\n",
      "This can also happens if the bandwidth is too small ",
      "(see the documentation for advice on the choice of h).\n",
      "You can disable this warning using the input `warnNA = FALSE`.")
  }

  warning(CondCopulas_warning_condition_base(
    message = message,
    subclass = "NA_ProducedWarning")
  )
}



#' @export
#'
#' @rdname plot.estimated_CKT_kernel
se.estimated_CKT_kernel <- function(object, progressBar = TRUE, ...)
{
  if ( !is.null(object$se) ){
    return(object$se)
  }

  if ( is.null(object$matrixSignPairs)){
    typeEstCKT = if(object$typeEstCKT == "wdm") 4 else object$typeEstCKT

    matrixSignsPairs = computeMatrixSignPairs(
      vectorX1 = object$X1, vectorX2 = object$X2, typeEstCKT = typeEstCKT)
  } else {
    matrixSignsPairs = object$matrixSignsPairs
  }

  temp <- compute_all_Gn_H_ii(vectorZ = object$Z,
                              vectorZToEstimate = object$newZ,
                              vector_hat_CKT_NP = object$estimatedCKT,
                              matrixSignsPairs = matrixSignsPairs,
                              h = object$h, kernel.name = "Epa", intK2 = 3/5,
                              progressBar = progressBar)

  n = length(object$X1)

  asympt_se_np = sqrt(temp$vect_H_ii) / sqrt(n * object$h)

  return (asympt_se_np)
}


#' @export
#'
#' @rdname plot.estimated_CKT_kernel
confint.estimated_CKT_kernel <- function(object, level = 0.95,
                                         progressBar = TRUE, ...)
{
  if (is.null(object$se)){
    object$se = se(object, progressBar = progressBar)
  }
  # Now se is available

  # Quantile of the standard normal distribution at level 1 - alpha / 2
  alpha = 1 - level
  q_1_alpha_2 <- stats::qnorm(1 - alpha / 2)

  estCKT = object$estimatedCKT
  nprime = length(estCKT)

  result <- matrix(nrow = nprime, ncol = 2)
  result[, 1] = pmax(-1, estCKT - q_1_alpha_2 * object$se)
  result[, 2] = pmin(1, estCKT + q_1_alpha_2 * object$se)

  colnames(result) <- paste(c(alpha / 2, 1 - alpha / 2) , "%")
  rownames(result) <- object$newZ

  return (result)
}


#' Methods for class `estimated_CKT_kernel`
#'
#' @param object,x an \code{S3} object of class \code{estimated_CKT_kernel}.
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 100
#' # This is a small example for performance reasons.
#' # For a better example, use:
#' # N = 800
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' newZ = seq(2, 10, by = 1)
#' estimatedCKT_kernel <- CKT.kernel(
#'    X1 = X1, X2 = X2, Z = Z,
#'    newZ = newZ, h = 0.2, kernel.name = "Epa", se = TRUE)
#'
#' se(estimatedCKT_kernel)
#' confint(estimatedCKT_kernel, level = 0.9)
#'
#' plot(estimatedCKT_kernel, confint = TRUE)
#'
#'
#' @export
plot.estimated_CKT_kernel <- function(x, confint = NULL, level = NULL,
                                      xlim = NULL, ylim = c(-1.2, 1.2),
                                      progressBar = TRUE,
                                      color_CKT = "black", color_confint = "red",
                                      ...)
{
  plot(x$newZ, x$estimatedCKT, type = "l", ylim = ylim, xlim = xlim,
       xlab = "z",
       ylab = "Conditional Kendall's tau given Z = z", col = color_CKT, ...)

  if (!isFALSE(confint)){

    # Easy case: a level is not specified but confint is available in the object.
    # then we just plot it.
    if (is.null(level) && !is.null(x$confint)){
      graphics::lines(x$newZ, x$confint[, 1], type = "l", col = color_confint)
      graphics::lines(x$newZ, x$confint[, 2], type = "l", col = color_confint)
    } else {
      # We will need to do computations. But first, let's see whether the user
      # intended to plot a confidence interval
      if (is.null(confint)){
        # If the user specifies a level then it shows intent to have a
        # confidence interval
        confint = !is.null(level)
      } else {
        # If the user specifies confint = TRUE but without giving an explicit
        # level, we use by default 95%
        if (is.null(level)){
          level = 0.95
        }
      }
      if (confint){
        x$confint = confint(x, level = level, progressBar = progressBar)
        graphics::lines(x$newZ, x$confint[, 1], type = "l", col = color_confint)
        graphics::lines(x$newZ, x$confint[, 2], type = "l", col = color_confint)
      }
    }
  }
}


compute_all_Gn_H_ii <- function(vectorZ, vectorZToEstimate, matrixSignsPairs,
                                vector_hat_CKT_NP,
                                h, kernel.name, intK2, progressBar){
  nprime = length(vectorZToEstimate)

  # 1 - Computation of G_n(z'_i) for all i
  vectorZToEstimate_arr = array(vectorZToEstimate)
  matrixSignsPairsSymmetrized = (matrixSignsPairs + t(matrixSignsPairs)) / 2
  if (progressBar){
    Gn_zipr = pbapply::pbapply(
      X = vectorZToEstimate_arr, MARGIN = 1,
      FUN = function(pointZ) {compute_vect_Gn_zipr(
        pointZ, vectorZ, h,
        kernel.name, matrixSignsPairsSymmetrized) } )
  } else {
    Gn_zipr = apply(
      X = vectorZToEstimate_arr, MARGIN = 1,
      FUN = function(pointZ) {compute_vect_Gn_zipr(
        pointZ, vectorZ, h,
        kernel.name, matrixSignsPairsSymmetrized) } )
  }

  # 2 - Computation of H_(i,i) under the hypothesis that all z'_i are distinct
  vect_H_ii = rep(NA, nprime)
  for (iprime in 1:nprime) {
    pointZ = vectorZToEstimate[iprime]
    listKh = computeWeights.univariate(
      vectorZ = vectorZ, h = h, pointZ = pointZ,
      kernel.name = kernel.name, normalization = FALSE)
    estimator_fZ = mean(listKh)
    vect_H_ii[iprime] = 4 * (intK2 / estimator_fZ) *
      abs(Gn_zipr[iprime] - (vector_hat_CKT_NP[iprime])^2)
  }

  return (list(Gn_zipr = Gn_zipr,
               vect_H_ii = vect_H_ii))
}

