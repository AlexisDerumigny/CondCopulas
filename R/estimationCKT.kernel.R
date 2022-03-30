

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
#' @noRd
#'
computeWeights.multivariate <- function(matrixZ, h, pointZ,
                                        kernel.name, normalization = TRUE)
{
  u = sweep(matrixZ, MARGIN = 2, STATS = pointZ)

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
CKT.kernelPointwise.univariate <- function(matrixSignsPairs, vectorZ,
                                           h, pointZ, kernel.name, typeEstCKT)
{
  listWeights = computeWeights.univariate(vectorZ, h, pointZ, kernel.name)
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

  return(estimate)
}



#' Estimate the conditional Kendall's tau of X1 and X2
#' at a fixed multivariate point Z = pointZ
#'
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
CKT.kernelPointwise.multivariate <- function(matrixSignsPairs, matrixZ,
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
  # if(! is.finite(estimate) ) {print(whichNonZero) ; print(listWeights) ; print(pointZ)}
  return(estimate)
}



# Estimation of CKT at multiple points ----------------------------------------------


#' Estimate the conditional Kendall's tau of X1 and X2
#' at different points
#'
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#'
#' @param observedZ vector of observed points of Z.
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
CKT.kernel.univariate <- function(matrixSignsPairs, observedZ,
                                  h, ZToEstimate,
                                  kernel.name = "Epa", typeEstCKT = 4,
                                  progressBar = TRUE)
{
  if (nrow(matrixSignsPairs) != ncol(matrixSignsPairs)){
    stop("matrixSignsPairs must be a square matrix.")
  } else if (nrow(matrixSignsPairs) != length(observedZ)){
    stop(paste0("observedZ must have the same length ",
                "as the number of rows of matrixSignsPairs."))
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
        pointZ = ZToEstimate[i], matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], vectorZ = observedZ,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  } else {
    estimates = apply(
      X = array(1:n_prime), MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.univariate(
        pointZ = ZToEstimate[i], matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], vectorZ = observedZ,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  }

  return(estimates)
}


#' Estimate the conditional Kendall's tau of X1 and X2
#' at different points
#'
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#'
#' @param observedZ matrix of observed points of Z.
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
CKT.kernel.multivariate <- function(matrixSignsPairs, observedZ,
                                    h, ZToEstimate,
                                    kernel.name = "Epa", typeEstCKT = 4,
                                    progressBar = TRUE)
{
  if (nrow(matrixSignsPairs) != ncol(matrixSignsPairs)){
    stop("matrixSignsPairs must be a square matrix.")
  } else if (nrow(matrixSignsPairs) != nrow(observedZ)){
    stop(paste0("observedZ and matrixSignsPairs must have",
                "the same number of rows."))
  } else if (ncol(observedZ) != ncol(ZToEstimate)){
    stop(paste0("observedZ and ZToEstimate must have",
                "the same number of columns."))
  }

  dim_Z = ncol(observedZ)
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
        pointZ = ZToEstimate[i,], matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], matrixZ = observedZ,
        kernel.name = kernel.name, typeEstCKT = typeEstCKT) } )
  } else {
    estimates = apply(
      X = 1:n_prime, MARGIN = 1,
      FUN = function(i) {CKT.kernelPointwise.multivariate(
        pointZ = ZToEstimate[i,], matrixSignsPairs = matrixSignsPairs,
        h = h_vect[i], matrixZ = observedZ,
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
#' (Derumigny,& Fermanian (2019)).
#'
#' @param observedX1 a vector of n observations of the first variable
#'
#' @param observedX2 a vector of n observations of the second variable
#'
#' @param observedZ a vector of n observations of the conditioning variable,
#' or a matrix with n rows of observations of the conditioning vector
#'
#' @param newZ the new data of observations of Z at which
#' the conditional Kendall's tau should be estimated.
#'
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' Default is 4. 1 and 3 produced biased estimator
#' while 2 does not attain the full range \eqn{[-1,1]}.
#'
#' @param methodCV method used for the cross-validation.
#' Possible choices are \code{leave-one-out} and \code{Kfolds}.
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
#' @param progressBar if \code{TRUE}, a progressbar for each h is displayed
#' to show the progress of the computation.
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
#' See \code{\link{CKT.hCV.l1out}} for manual selection of the bandwidth
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
#'    observedX1 = X1, observedX2 = X2, observedZ = Z,
#'    newZ = newZ, h = 0.07, kernel.name = "Epa")$estCKT
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in red)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col="black",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_kernel, col = "red")
#'
#' @export
#'
CKT.kernel <- function(observedX1, observedX2, observedZ, newZ,
                       h, kernel.name = "Epa",
                       methodCV = "Kfolds",
                       Kfolds = 5, nPairs = 10*length(observedX1),
                       typeEstCKT = 4, progressBar = TRUE)
{
  matrixSignsPairs = computeMatrixSignPairs(
    vectorX1 = observedX1, vectorX2 = observedX2, typeEstCKT = typeEstCKT)

  if (length(h) == 1){
    finalh = h
  } else {

    # Do the cross-validation

    switch (
      methodCV,

      "Kfolds" = {
        resultCV = CKT.hCV.Kfolds(
          range_h = h, matrixSignsPairs = matrixSignsPairs,
          observedZ = observedZ, ZToEstimate = newZ,
          typeEstCKT = typeEstCKT, kernel.name = kernel.name, Kfolds = Kfolds,
          progressBar = progressBar)
      },

      "leave-one-out" = {
        resultCV = CKT.hCV.l1out(
          observedX1 = observedX1, observedX2 = observedX2,
          observedZ = observedZ,
          range_h = h, matrixSignsPairs = matrixSignsPairs,
          typeEstCKT = typeEstCKT, kernel.name = kernel.name,
          nPairs = nPairs, progressBar = progressBar)
      }

    )
    finalh = resultCV$hCV
  }

  # We finally computed the estimated conditional Kendall's tau
  # using the selected value for h and returns it

  if (is.vector(observedZ)){
    estCKT = CKT.kernel.univariate(
      matrixSignsPairs = matrixSignsPairs, observedZ = observedZ,
      h = finalh, ZToEstimate = newZ,
      kernel.name = kernel.name, typeEstCKT = typeEstCKT,
      progressBar = TRUE)

    return (list(estimatedCKT = estCKT, h = finalh))
  } else {

    estCKT = CKT.kernel.multivariate(
      matrixSignsPairs = matrixSignsPairs, observedZ = observedZ,
      h = finalh, ZToEstimate = newZ,
      kernel.name = kernel.name, typeEstCKT = typeEstCKT,
      progressBar = TRUE)

    return (list(estimatedCKT = estCKT, h = finalh))
  }
}


