
#' Compute a kernel-based estimator of the conditional copula
#'
#' Assuming that we observe a sample \eqn{(X_{i,1}, X_{i,2}, X_{i,3}), i=1, \dots, n},
#' this function returns a array
#' \eqn{\hat C_{1,2|3}(u_1, u_2 | X_3 = x_3)}
#' for each choice of (u_1, u_2, x_3).
#'
#' @param X1,X2,X3 vectors of observations of size \code{n}
#' @param U1_  a vector of numbers in [0, 1]
#' @param U2_  a vector of numbers in [0, 1]
#' @param newX3  a vector of new values for the conditioning variable \code{X3}
#' @param kernel a character string describing the kernel to be used.
#' Possible choices are \code{Gaussian}, \code{Triangular} and \code{Epanechnikov}.
#' @param h the bandwidth to use in the estimation.
#'
#' @param observedX1,observedX2,observedX3 old parameter names for \code{X1},
#' \code{X2}, \code{X3}. Support for this will be removed at a later version.
#'
#' @return An array of dimension \code{(length(U1_, U2_, newX3))}
#' whose element in position (i, j, k) is
#' \eqn{\hat C_{1,2|3}(u_1, u_2 | X_3 = x_3)}
#' where \eqn{u_1} = U1_[i], \eqn{u_2} = U2_[j] and \eqn{x_3} = newX3[k]
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#' \doi{10.1515/demo-2017-0011}
#'
#' @seealso \code{\link{estimateParCondCopula}} for estimating a conditional
#' copula in a parametric setting ( = where the conditional copula is assumed to
#' belong to a parametric class).
#' \code{\link{simpA.NP}} for a test that this conditional copula is
#' constant with respect to the value \eqn{x_3} of the conditioning variable.
#'
#' @examples
#' # We simulate from a conditional copula
#' N = 500
#' X3 = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 0.9 * pnorm(X3, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 3,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' # We do the estimation
#' grid = c(0.2, 0.4, 0.6, 0.8)
#' arrayEst = estimateNPCondCopula(
#'   X1 = X1, X2 = X2, X3 = X3,
#'   U1_ = grid, U2_ = grid, newX3 = c(2, 5, 7),
#'   kernel = "Gaussian", h = 0.8)
#' arrayEst
#'
#' @export
#'
estimateNPCondCopula <- function(X1 = NULL, X2 = NULL, X3 = NULL,
                                 U1_, U2_, newX3, kernel, h,
                                 observedX1 = NULL, observedX2 = NULL, observedX3 = NULL)
{
  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  env = environment()
  .observedX1X2_to_X1X2(env)
  .observedX3_to_X3(env)

  .checkSame_nobs_X1X2X3(X1, X2, X3)
  .checkUnivX1X2X3(X1, X2, X3)

  if (any(U1_ < 0 | U1_ > 1 | U2_ < 0 | U2_ > 1)){
    stop(errorCondition(
      message = "All the elements of U1_ and U2_ should all be in the interval [0,1]",
      class = "InvalidInputError"
    ))
  }

  nNewPoints1 = length(U1_)
  nNewPoints2 = length(U2_)
  nNewPoints3 = length(newX3)

  matrixK3 = computeKernelMatrix(
    observedX = X3, newX = newX3, kernel = kernel, h = h)

  # Sometimes, h is too small, and then there are some points newX3, that are
  # such that all all the kernel weights that correspond to them are 0.
  # These points should be found out, and a warning be raised.
  # They would cause anyway a problem further downstream with estimateCondQuantiles(),
  # so it is better to catch them now than later.
  which_zeroWeights = which(colSums(matrixK3) == 0)
  which_nonZeroWeights = which(colSums(matrixK3) > 0)

  if (length(which_zeroWeights) > 0){
    warning(CondCopulas_warning_condition_base(
      message = paste0(
        length(which_zeroWeights), " ",
        "points (out of ", length(newX3), ") have kernel weights that are all 0. ",
        "This is usually caused by a too small bandwidth and/or ",
        "by points that are too far away from the distribution of the conditioning variable."
      ),
      problematicPoints = newX3[which_zeroWeights],
      subclass = "ZeroWeights_KernelWarning"
    ))

    if (length(which_nonZeroWeights) == 0){
      # We now that all points have 0 weight. Estimation is therefore impossible.
      # The best we can do is return an array with NAs everywhere.

      array_C_IJ = array(dim = c(nNewPoints1, nNewPoints2, nNewPoints3))
      return (array_C_IJ)
    }
  }

  # Computation of conditional quantiles
  matrix_condQuantile13 = estimateCondQuantiles(
    observedX1 = X1, probsX1 = U1_,
    matrixK3 = matrixK3[ , which_nonZeroWeights, drop = FALSE])

  matrix_condQuantile23 = estimateCondQuantiles(
    observedX1 = X2, probsX1 = U2_,
    matrixK3 = matrixK3[ , which_nonZeroWeights, drop = FALSE])


  # Computation of C_(1,2|3) ( U_(i,1) , U_(j,2) | X_(k,3) )
  array_C_IJ = array(dim = c(nNewPoints1, nNewPoints2, nNewPoints3))

  for (k in 1:length(which_nonZeroWeights))
  {
    k_original = which_nonZeroWeights[k]
    for (i in 1:nNewPoints1)
    {
      for (j in 1:nNewPoints2)
      {
        array_C_IJ[i, j, k_original] = sum(
          matrixK3[ , k_original] * as.numeric(X1 <= matrix_condQuantile13[i , k] &
                                                 X2 <= matrix_condQuantile23[j , k] ) )
      }
    }
    array_C_IJ[, , k_original] = array_C_IJ[, , k_original] / sum(matrixK3[, k_original])
  }

  return (array_C_IJ)
}


#' Estimation of parametric conditional copulas
#'
#' The function \code{estimateParCondCopula} computes an estimate
#' of the conditional parameters in a conditional parametric copula model, i.e.
#' \deqn{C_{X_1, X_2 | X_3 = x_3} = C_{\theta(x_3)},}
#' for some parametric family \eqn{(C_\theta)}, some conditional
#' parameter \eqn{\theta(x_3)}, and a three-dimensional random
#' vector \eqn{(X_1, X_2, X_3)}. Remember that \eqn{C_{X_1,X_2 | X_3 = x_3}}
#' denotes the conditional copula of \eqn{X_1} and \eqn{X_2} given \eqn{X_3 = x_3}.
#'
#' @param X1 a vector of \code{n} observations
#' of the first conditioned variable
#'
#' @param X2 a vector of \code{n} observations
#' of the second conditioned variable
#'
#' @param X3 a vector of \code{n} observations
#' of the conditioning variable
#'
#' @param newX3 a vector of new observations of \eqn{X3}
#'
#' @param family an integer indicating the parametric family of copulas to be used,
#' following the conventions of the \link[VineCopula:VineCopula]{VineCopula}
#' package, see e.g. \link[VineCopula:BiCop]{VineCopula::BiCop}.
#'
#' @param method the method of estimation of the conditional parameters.
#' Can be \code{"mle"} for maximum likelihood estimation
#' or \code{"itau"} for estimation by inversion of Kendall's tau.
#'
#' @param h bandwidth to be chosen
#'
#' @param observedX1,observedX2,observedX3 old parameter names for \code{X1},
#' \code{X2}, \code{X3}. Support for this will be removed at a later version.
#'
#'
#' @return a vector of size \code{length(newX3)} containing
#' the estimated conditional copula parameters for each value of \code{newX3}.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#' \doi{10.1515/demo-2017-0011}
#'
#' @seealso \code{\link{estimateNPCondCopula}} for estimating a conditional
#' copula in a nonparametric setting ( = without parametric assumption on the
#' conditional copula).
#' \code{\link{simpA.param}} for a test that this conditional copula is
#' constant with respect to the value \eqn{x_3} of the conditioning variable.
#'
#'
#' @examples
#'
#' # We simulate from a conditional copula
#' N = 500
#'
#' X3 = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 0.9 * pnorm(X3, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(
#'     N=N , family = 1, par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' gridnewX3 = seq(2, 8, by = 1)
#' conditionalTauNewX3 = 0.9 * pnorm(gridnewX3, mean = 5, sd = 2)
#'
#' vecEstimatedThetas = estimateParCondCopula(
#'   X1 = X1, X2 = X2, X3 = X3,
#'   newX3 = gridnewX3, family = 1, h = 0.1)
#'
#' # Estimated conditional parameters
#' vecEstimatedThetas
#' # True conditional parameters
#' VineCopula::BiCopTau2Par(1 , conditionalTauNewX3 )
#'
#' # Estimated conditional Kendall's tau
#' VineCopula::BiCopPar2Tau(1 , vecEstimatedThetas )
#' # True conditional Kendall's tau
#' conditionalTauNewX3
#'
#'
#' @export
#'
estimateParCondCopula <- function(X1 = NULL, X2 = NULL, X3 = NULL,
                                  newX3, family, method = "mle", h,
                                  observedX1 = NULL, observedX2 = NULL, observedX3 = NULL
)
{
  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  env = environment()
  .observedX1X2_to_X1X2(env)
  .observedX3_to_X3(env)

  .checkSame_nobs_X1X2X3(X1, X2, X3)
  .checkUnivX1X2X3(X1, X2, X3)

  # Computation of pseudo-observations
  U1 = stats::ecdf(X1)(X1)
  U2 = stats::ecdf(X2)(X2)
  U3 = stats::ecdf(X3)(X3)
  newU3 = stats::ecdf(X3)(newX3)

  # Computation of the kernel
  matrixK3 = computeKernelMatrix(observedX = U3, newX = U3, kernel = "Epanechnikov", h = h)

  # Computation of conditional cumulative distribution functions
  Z1_J = estimateCondCDF_vec(observedX1 = U1, newX1 = U1, matrixK3 = matrixK3)

  Z2_J = estimateCondCDF_vec(observedX1 = U2, newX1 = U2, matrixK3 = matrixK3)

  theta_xJ = estimateParCondCopula_ZIJ(
    Z1_J = Z1_J, Z2_J = Z2_J, observedX3 = U3, newX3 = newU3,
    family = family, method = method, h = h)

  return (theta_xJ)
}


#' Estimation of a parametric conditional copula using conditional pseudos-observations
#'
#' The function \code{estimateParCondCopula_ZIJ} is an auxiliary function
#' that is called when conditional pseudos-observations are already
#' available when one wants to estimate a parametric conditional copula.
#'
#' @param Z1_J the conditional pseudos-observations of the first variable,
#' i.e. \eqn{\hat F_{1|J}( x_{i,1} | x_J = x_{i,J})} for \eqn{i=1,\dots, n}.
#'
#' @param Z2_J the conditional pseudos-observations of the second variable,
#' i.e. \eqn{\hat F_{2|J}( x_{i,2} | x_J = x_{i,J})} for \eqn{i=1,\dots, n}.
#'
#'
#' @rdname estimateParCondCopula
#' @export
#'
estimateParCondCopula_ZIJ <- function (Z1_J, Z2_J, observedX3,
                                       newX3, family, method, h)
{
  # Computation of conditional parameters
  nNewPoints = length(newX3)
  theta_xJ = rep(NA, nNewPoints)

  for (i in 1:nNewPoints)
  {
    diff = abs(observedX3 - newX3[i])/h
    # To do later: implement other kernels
    weights = 3/4 * (1-diff^2)/h * as.numeric(diff < 1)
    theta_xJ[i] =
      try(VineCopula::BiCopEst(
        u1 = Z1_J[which(weights > 0)], u2 = Z2_J[which(weights > 0)],
        family = family , method = method,
        weights = weights[which(weights > 0)]
      )$par , silent = FALSE)

    if (!is.finite(theta_xJ[i]) | !is.numeric(theta_xJ[i]))
    {
      warnings(paste0("Problem of estimation in estimationParCondCopulas",
                      ". For newX3 =", newX3[i],
                      " , theta(xJ) is estimated by", theta_xJ[i]))
    }
  }
  theta_xJ = as.numeric(theta_xJ)

  return (theta_xJ)
}

