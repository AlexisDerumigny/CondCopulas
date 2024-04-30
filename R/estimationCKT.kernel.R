

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
    # For faster computation, only uses points with non-zero
    # values for the kernel
    u = sweep(matrixZ, MARGIN = 2, STATS = pointZ)
    isSmaller_h = apply(X = u, MARGIN = 1, FUN = function(x){
      return (all(abs(x) <= h))
    })
    whichNonZero = which( isSmaller_h )
    listWeights = computeWeights.multivariate(
      matrixZ[whichNonZero, ], h, pointZ, kernel.name)

  } else {
    listWeights = computeWeights.multivariate(
      matrixZ[ , ], h, pointZ, kernel.name)
  }


  if(typeEstCKT == "wdm"){

    estimate = wdm::wdm(x = X1, y = X2, weights = listWeights, method = "kendall")

  } else {
    matrixWeights = outer(listWeights, listWeights)

    switch (
      typeEstCKT,
      # 1
      { estimate =
        4 * sum(matrixWeights * matrixSignsPairs[whichNonZero, whichNonZero]) - 1 } ,
      # 2
      { estimate =sum(matrixWeights * matrixSignsPairs[whichNonZero, whichNonZero]) },
      # 3
      { estimate =
        1 - 4 * sum(matrixWeights * matrixSignsPairs[whichNonZero, whichNonZero]) },
      # 4
      { estimate =
        sum(matrixWeights * matrixSignsPairs[whichNonZero, whichNonZero]) /
        (1 - sum(listWeights^2)) },

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
    if (nrow(matrixSignsPairs) != ncol(matrixSignsPairs)){
      stop("matrixSignsPairs must be a square matrix.")
    } else if (nrow(matrixSignsPairs) != length(Z)){
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
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i], matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], vectorZ = Z,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  } else {
    estimates = apply(
      X = array(1:n_prime), MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.univariate(
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i], matrixSignsPairs = matrixSignsPairs,
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
    if (nrow(matrixSignsPairs) != ncol(matrixSignsPairs)){
      stop("matrixSignsPairs must be a square matrix.")
    } else if (nrow(matrixSignsPairs) != nrow(Z)){
      stop(paste0("Z and matrixSignsPairs must have",
                  "the same number of rows."))
    } else if (ncol(Z) != ncol(ZToEstimate)){
      stop(paste0("Z and ZToEstimate must have",
                  "the same number of columns."))
    }
  }

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
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i,], matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], matrixZ = Z,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  } else {
    estimates = apply(
      X = 1:n_prime, MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.multivariate(
        X1 = X1, X2 = X2, pointZ = ZToEstimate[i,], matrixSignsPairs = matrixSignsPairs,
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
#' Indeed the conditioning points should be closer and closer to the point of interest \eqn{z}
#' (small \code{h}) and more and more numerous (\code{h} tending to 0 slowly enough).
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
#'
#' @param newZ the new data of observations of Z at which
#' the conditional Kendall's tau should be estimated.
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
#' if \code{methodCV = "leave-one-out"}.
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
#' @return a list with two components
#' \itemize{
#'    \item \code{estimatedCKT} the vector of size \code{NROW(newZ)}
#'    containing the values of the estimated conditional Kendall's tau.
#'
#'    \item \code{finalh} the bandwidth \code{h} that was finally used
#'    for kernel smoothing (either the one specified by the user
#'    or the one chosen by cross-validation if multiple bandwidths were given.)
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
#' N = 800
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
#' @export
#'
CKT.kernel <- function(X1 = NULL, X2 = NULL, Z = NULL, newZ,
                       h, kernel.name = "Epa",
                       methodCV = "Kfolds",
                       Kfolds = 5, nPairs = 10*length(observedX1),
                       typeEstCKT = "wdm", progressBar = TRUE,
                       observedX1 = NULL, observedX2 = NULL, observedZ = NULL)
{
  if (length(newZ) == 0){
    warning("'newZ' is of length 0, therefore, no estimation is done.")
    return (numeric(0))
  }

  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  if (is.null(X1)){
    if (is.null(observedX1)){
      stop("X1 must be of non-zero length.")
    } else {
      X1 = observedX1
    }
  }
  if (is.null(X2)){
    if (is.null(observedX2)){
      stop("X2 must be of non-zero length.")
    } else {
      X2 = observedX2
    }
  }
  if (is.null(Z)){
    if (is.null(observedX2)){
      stop("Z must be of non-zero length.")
    } else {
      Z = observedZ
    }
  }

  # Checking for same number of observations
  .checkSame_nobs_X1X2Z(X1, X2, Z)

  # Checking that X1 and X2 are univariate
  .checkUnivX1X2(X1, X2)
  X1 = as.numeric(X1)
  X2 = as.numeric(X2)

  # Putting Z as a column vector if it has only one column
  if (NCOL(Z) == 1 && is.matrix(Z)){
    Z = as.numeric(Z)
  }

  if (typeEstCKT == "wdm") {
    matrixSignsPairs = NULL
  } else {
    matrixSignsPairs = computeMatrixSignPairs(
      vectorX1 = X1, vectorX2 = X2, typeEstCKT = typeEstCKT)
  }

  if (length(h) == 1){
    finalh = h
  } else {

    # Do the cross-validation

    switch (
      methodCV,

      "Kfolds" = {
        resultCV = CKT.hCV.Kfolds(
          observedX1 = X1, observedX2 = X2, observedZ = Z,
          range_h = h, matrixSignsPairs = matrixSignsPairs,
          ZToEstimate = newZ,
          typeEstCKT = typeEstCKT, kernel.name = kernel.name, Kfolds = Kfolds,
          progressBar = progressBar > 1)
      },

      "leave-one-out" = {
        resultCV = CKT.hCV.l1out(
          observedX1 = X1, observedX2 = X2, observedZ = Z,
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

  if (anyNA(estCKT)){
    if (!anyNA(X1) && !anyNA(X2) && !anyNA(Z)){
      warning("NA in estimated CKT. ",
              "This often happens when the bandwidth h is too small, ",
              "consider using a bigger bandwidth ",
              "(see the documentation for advice on the choice of h).")
    }
  }

  return (list(estimatedCKT = estCKT, h = finalh))
}


