
#' Compute a kernel-based estimator of the conditional copula
#'
#' Assuming that we observe a sample \eqn{(X_{i,1}, X_{i,2}, X_{i,3}), i=1, \dots, n},
#' this function returns a array
#' \eqn{\hat C_{1,2|3}(u_1, u_2 | X_3 = x_3)}
#' for each choice of (u_1, u_2, x_3).
#'
#' @param observedX1 a vector of observations of size n
#' @param observedX2 a vector of observations of size n
#' @param observedX3 a vector of observations of size n
#' @param U1_  a vector of numbers in [0, 1]
#' @param U2_  a vector of numbers in [0, 1]
#' @param newX3  a vector of new values for the conditioning variable X3
#' @param kernel a character string describing the kernel to be used.
#' Possible choices are \code{Gaussian}, \code{Triangular} and \code{Epanechnikov}.
#' @param h the bandwidth to use in the estimation.
#'
#' @return An array of dimension \code{(length(U1_, U2_, newX3))}
#' whose element in position (i, j, k) is
#' \eqn{\hat C_{1,2|3}(u_1, u_2 | X_3 = x_3)}
#' where \eqn{u_1} = U1_[i], \eqn{u_2} = U2_[j] and \eqn{x_3} = newX3[k]
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
#' arrayEst = estimateNPCondCopula(observedX1 = X1,
#'   observedX2 = X2, observedX3 = X3,
#'   U1_ = grid, U2_ = grid, newX3 = c(2, 5, 7),
#'   kernel = "Gaussian", h = 0.8)
#' arrayEst
#'
#' @export
#'
estimateNPCondCopula <- function(observedX1, observedX2, observedX3,
                                 U1_, U2_, newX3, kernel, h)
{
  if (length(observedX1) != length(observedX2)) {
    stop("Error: observedX1 and observedX2 have different lengths.")
  }
  if (length(observedX1) != length(observedX3)) {
    stop("Error: observedX1 and observedX3 have different lengths.")
  }

  nNewPoints1 = length(U1_)
  nNewPoints2 = length(U2_)
  nNewPoints3 = length(newX3)

  matrixK3 = computeKernelMatrix(
    observedX = observedX3, newX = newX3, kernel = kernel, h = h)

  # Computation of conditional quantiles
  matrix_condQuantile13 = estimateCondQuantiles(
    observedX1 = observedX1, probsX1 = U1_, matrixK3 = matrixK3)

  matrix_condQuantile23 = estimateCondQuantiles(
    observedX1 = observedX2, probsX1 = U2_, matrixK3 = matrixK3)


  # Computation of C_(1,2|3) ( U_(i,1) , U_(j,2) | X_(k,3) )
  array_C_IJ = array(dim = c(nNewPoints1, nNewPoints2, nNewPoints3))

  for (k in 1:nNewPoints3)
  {
    for (i in 1:nNewPoints1)
    {
      for (j in 1:nNewPoints2)
      {
        array_C_IJ[i, j, k] =
          sum(matrixK3[ , k] * as.numeric(observedX1 <= matrix_condQuantile13[i , k] &
                                            observedX2 <= matrix_condQuantile23[j , k] ) )
      }
    }
    array_C_IJ[, , k] = array_C_IJ[, , k] / sum(matrixK3[, k])
  }

  return (array_C_IJ)
}


#' Estimation of parametric conditional copulas
#'
#' This function computes an estimate of the conditional parameters
#' of a conditional parametric copula model.
#'
#' We assume that we have a three-dimensional vector (X1, X2, X3).
#' We are interested in the estimation of the conditional copula of (X1, X2) given X3.
#' This function takes in parameter three vectors of observations of X1, X2, X3,
#' and a vector of conditioning points \code{newX3}
#' and compute estimators of the corresponding parametric conditional copulas.
#'
#' @param observedX1 a vector of n observations of the first variable
#' @param observedX2 a vector of n observations of the second variable
#' @param observedX3 a vector of n observations of the conditioning variable
#'
#' @param newX3 a vector of new observations of X3
#' @param family an integer indicating the parametric family of copulas to be used,
#' following the conventions of the package \code{VineCopula}.
#' @param method the method of estimation of the conditional parameters.
#' Can be \code{"mle"} for maximum likelihood estimation
#' or \code{"itau"} for estimation by inversion of Kendall's tau.
#' @param h bandwidth to be chosen
#'
#' @return a vector of size \code{length(newX3)} containing
#' the estimated conditional copula parameters for each value of newU3.
#'
#' @examples
#'
#' # We simulate from a conditional copula
#' N = 5000
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
#'   observedX1 = X1, observedX2 = X2, observedX3 = X3,
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
estimateParCondCopula <- function (observedX1, observedX2, observedX3,
                                   newX3, family, method = "mle", h)
{
  # Computation of pseudo-observations
  U1 = stats::ecdf(observedX1)(observedX1)
  U2 = stats::ecdf(observedX2)(observedX2)
  U3 = stats::ecdf(observedX3)(observedX3)
  newU3 = stats::ecdf(observedX3)(newX3)
  # Computation of the kernel
  matrixK3 = computeKernelMatrix(observedX = U3, newX = U3, kernel = "Epanechnikov", h = h)

  # Computation of conditional cumulative distribution functions
  Z1_J = estimateCondCDF_vec(observedX1 = U1, newX1 = U1, matrixK3 = matrixK3)

  Z2_J = estimateCondCDF_vec(observedX1 = U2, newX1 = U2, matrixK3 = matrixK3)

  theta_xJ = estimateParCondCopula_ZIJ(
    Z1_J = Z1_J, Z2_J = Z2_J, U3 = U3, newU3 = newU3,
    family = family, method = method, h = h)

  return (theta_xJ)
}


#' Estimation of a parametric conditional copula using conditional pseudos-observations
#'
#' This is an auxiliary function that is called when conditional pseudos-observations
#' are already available to estimate a parametric conditional copula.
#'
#' @param Z1_J the conditional pseudos-observations of the first variable,
#' i.e. \eqn{\hat F_{1|J}( x_{i,1} | x_J = x_{i,J})} for \eqn{i=1,\dots, n}.
#' @param Z2_J the conditional pseudos-observations of the second variable,
#' i.e. \eqn{\hat F_{2|J}( x_{i,2} | x_J = x_{i,J})} for \eqn{i=1,\dots, n}.
#' @param U3 the pseudos-observations of the conditional variable.
#' @param newU3 new observations of the conditionning variable on which we would like
#' to estimate the conditional parameter
#' @param family the family of conditional copulas used.
#' @param method method of estimation for the parametric conditional copula.
#' @param h the bandwidth.
#'
#' @return a vector of size \code{length(newX3)} containing
#' the estimated conditional copula parameters for each value of newU3.
#'
#' @export
#'
estimateParCondCopula_ZIJ <- function (Z1_J, Z2_J, U3,
                                       newU3, family, method, h)
{

  # Computation of conditional parameters
  nNewPoints = length(newU3)
  theta_xJ = rep(NA, nNewPoints)

  for (i in 1:nNewPoints)
  {
    diff = abs(U3 - newU3[i])/h
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
                      ". For newU3 =", newU3[i],
                      " , theta(xJ) is estimated by", theta_xJ[i]))
    }
  }
  theta_xJ = as.numeric(theta_xJ)

  return (theta_xJ)
}

