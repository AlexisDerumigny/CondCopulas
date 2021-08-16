
#' Computing the kernel matrix
#'
#' This function computes a matrix of dimensions \code{(length(observedX3), length(newX3))},
#' whose element at coordinate \code{(i,j)} is
#' \eqn{ K_{h}(}\code{observedX3}\eqn{[i] - }\code{newX3}\eqn{[j] )},
#' where \eqn{K_h(x) := K(x/h) / h} and \eqn{K} is the \code{kernel}.
#'
#' @param observedX a numeric vector of observations of X3.
#' on the interval \eqn{[0,1]}.
#' @param newX a numeric vector of points of X3.
#' @param kernel a character string describing the kernel to be used.
#' Possible choices are \code{Gaussian}, \code{Triangular} and \code{Epanechnikov}.
#' @param h the bandwidth
#'
#' @return a numeric matrix of dimensions \code{(length(observedX), length(newX))}
#'
#' @seealso \code{\link{estimateCondCDF_matrix}}, \code{\link{estimateCondCDF_vec}},
#'
#'
#' @examples
#' Y = MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = cbind(c(1, 0.9), c(0.9, 1)))
#' matrixK = computeKernelMatrix(observedX = Y[,2], newX = c(0, 1, 2.5),
#' kernel = "Gaussian", h = 0.8)
#'
#' # To have an estimator of the conditional expectation of Y1 given Y2 = 0, 1, 2.5
#' Y[,1] * matrixK[,1] / sum(matrixK[,1])
#' Y[,1] * matrixK[,2] / sum(matrixK[,2])
#' Y[,1] * matrixK[,3] / sum(matrixK[,3])
#'
#' @export
#'
computeKernelMatrix <- function(observedX, newX, kernel, h)
{
  matrixK = outer(observedX , newX , function(x,y) {return(x-y)} )
  if (kernel == "Gaussian"){
    matrixK = exp( - 0.5 * ( (matrixK / h)^2) ) / (h * sqrt(2*pi) )
  }
  else if (kernel == "Triangular") {
    matrixK_norm = abs(matrixK) / h
    matrixK = (1 - matrixK_norm) / h * as.numeric(matrixK_norm < 1)
  }
  else if (kernel == "Epanechnikov") {
    matrixK_norm = abs(matrixK)/h
    matrixK = 3 / 4 * (1 - matrixK_norm^2) / h *
      as.numeric(matrixK_norm < 1)
  } else {
    stop(paste0("Error, unknown kernel:", kernel,
                " . Acceptable kernels are Gaussian, Triangular and Epanechnikov"))
  }
  return (matrixK)
}


#' Compute kernel-based conditional marginal (univariate) cdfs
#'
#' This function computes an estimate of the conditional (marginal) cdf
#' of X1 given a conditioning variable X3.
#'
#' This function is supposed to be used with \code{\link{computeKernelMatrix}}.
#' Assume that we observe a sample \eqn{(X_{i,1}, X_{i,3}), i=1, \dots, n}.
#' We want to estimate the conditional cdf of \eqn{X_1} given \eqn{X_3 = x_3}
#' at point \eqn{x_1} using the following kernel-based estimator
#' \deqn{\hat P(X_1 \le x_1 | X_3 = x_3)
#' := \frac{\sum_{l=1}^n 1 \{X_{l,1} \leq x_1 \} K_h(X_{l,3} - x_3)}
#' {\sum_{l=1}^n K_h(X_{l,3} - x_3)},}
#' for every \eqn{x_1} in \code{newX1} and every \eqn{x_3} in \code{newX3}.
#' The \code{matrixK3} should be a matrix of the values \eqn{K_h(X_{l,3} - x_3)}
#' such as the one produced by
#' \code{\link{computeKernelMatrix}(observedX3, newX3, kernel, h)}.
#'
#' @param observedX1 a sample of observations of X1 of size \code{n}
#' @param newX1 a sample of new points for the variable X1, of size \code{p1}
#' @param matrixK3 a matrix of kernel values of dimension \code{(p3, n)}
#' \eqn{\big(K_h(X3[i] - U3[j])\big)_{i,j}}
#' such as given by \code{\link{computeKernelMatrix}}.
#'
#' @return A matrix of dimensions \code{(p1 = length(newX), p3 = length(matrixK3[,1]))} of estimators
#' \eqn{\hat P(X_1 \leq x_1 | X_3 = x_3)} for every possible choices
#' of \eqn{(x_1, x_3)}.
#'
#' @examples
#' Y = MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = cbind(c(1, 0.9), c(0.9, 1)))
#' newY1 = seq(-1, 1, by = 0.5)
#' newY2 = c(0, 1, 2)
#' matrixK = computeKernelMatrix(observedX = Y[,2], newX = newY2,
#'   kernel = "Gaussian", h = 0.8)
#' # In this matrix, there are the estimated conditionl cdf at points given by newY1
#' # conditionally to the points given by newY2.
#' matrixCondCDF = estimateCondCDF_matrix(observedX1 = Y[,1],
#'   newX1 = newY1, matrixK)
#' matrixCondCDF
#'
#' @export
#'
estimateCondCDF_matrix <- function(observedX1, newX1, matrixK3)
{
  nNewPoints1 = length(newX1)
  nNewPoints2 = dim(matrixK3)[2]
  stopifnot(length(observedX1) == dim(matrixK3)[1])

  matrix_mn_1 = matrix(nrow = nNewPoints1, ncol = nNewPoints2)
  for (i in 1:nNewPoints1)
  {
    for (j in 1:nNewPoints2)
    {
      matrix_mn_1[i,j] =
        matrixK3[,j] %*% (as.numeric(observedX1 <= newX1[i])) / sum(matrixK3[,j])
    }
  }
  return (matrix_mn_1)
}

#
# To do later
# #' Compute kernel-based conditional marginal (multivariate) cdfs
# #'
# #' This function computes an estimate of the conditional (marginal) cdf
# #' of X = \eqn{{X_1, ..., X_p}} given a conditioning vector of variables, X = \eqn{{X_{p+1},
# #' ..., X_d}}.
# #'
# #' @param ObservedX A numeric matrix of p observations of \eqn{ X_i \in {X_1, X_2, ..., X_p}} of
# #' sample size \code{n}.
# #' @param newX A numeric matrix of p new point of \eqn{ X_i \in {X_1, X_2, ..., X_p}} of
# #' sample size \code{p1}
# #' @param HDmatrixK A matrix of kernel values, such as given by \code{\link{computeKernelMatrixHD}}.
# #'
# #' @return
# #'
# #' @keywords internal
# #'
# estimateCondCDF_matrixHD <- function(observedX, newX, HDmatrixK)
# {
#
# }

#' Compute kernel-based conditional marginal (univariate) cdfs
#'
#' This function computes an estimate of the conditional (marginal) cdf
#' of X1 given a conditioning variable X3.
#'
#' This function is supposed to be used with \code{\link{computeKernelMatrix}}.
#' Assume that we observe a sample \eqn{(X_{i,1}, X_{i,3}), i=1, \dots, n}.
#' We want to estimate the conditional cdf of \eqn{X_1} given \eqn{X_3 = x_3}
#' at point \eqn{x_1} using the following kernel-based estimator
#' \deqn{\hat P(X_1 \leq x_1 | X_3 = x_3)
#' := \frac{\sum_{l=1}^n 1 \{X_{l,1} \leq x_1 \} K_h(X_{l,3} - x_3)}
#' {\sum_{l=1}^n K_h(X_{l,3} - x_3)},}
#' for every couple \eqn{(x_{j,1}, x_{j,3})} where
#' \eqn{x_{j,1}} in \code{newX1} and \eqn{x_{j,3}} in \code{newX3}.
#' The \code{matrixK3} should be a matrix of the values \eqn{K_h(X_{l,3} - x_3)}
#' such as the one produced by
#' \code{\link{computeKernelMatrix}(observedX3, newX3, kernel, h)}.
#'
#' @param observedX1 a sample of observations of X1 of size n
#' @param newX1 a sample of new points for the variable X1, of size p1
#' @param matrixK3 a matrix of kernel values of dimension (p2 , n)
#' \eqn{\big(K_h(X3[i] - U3[j])\big)_{i,j}}
#' such as given by \code{\link{computeKernelMatrix}}.
#'
#' @return It returns a vector of length \code{newX1} of estimators
#' \eqn{\hat P(X_1 \leq x_1 | X_3 = x_3)}
#' for every couple \eqn{(x_{j,1}, x_{j,3})}.
#'
#' @examples
#' Y = MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = cbind(c(1, 0.9), c(0.9, 1)))
#' newY1 = seq(-1, 1, by = 0.5)
#' newY2 = newY1
#' matrixK = computeKernelMatrix(observedX = Y[,2], newX = newY2,
#'   kernel = "Gaussian", h = 0.8)
#' vecCondCDF = estimateCondCDF_vec(observedX1 = Y[,1],
#'   newX1 = newY1, matrixK)
#' vecCondCDF
#'
#' @export
#'
estimateCondCDF_vec <- function(observedX1, newX1, matrixK3)
{
  nNewPoints = length(newX1)
  stopifnot(length(observedX1) == dim(matrixK3)[1] ,
            nNewPoints == dim(matrixK3)[2])

  vect_mn_1 = rep(NA , nNewPoints)
  for (i in 1:nNewPoints)
  {
    vect_mn_1[i] = matrixK3[,i] %*% (as.numeric(observedX1 <= newX1[i])) / sum(matrixK3[,i])
  }
  vect_mn_1 = pmax(0.000001,pmin(0.999999,vect_mn_1))
  return (vect_mn_1)
}


#' Compute kernel-based conditional quantiles
#'
#' This function is supposed to be used with \code{\link{computeKernelMatrix}}.
#' Assume that we observe a sample \eqn{(X_{i,1}, X_{i,3}), i=1, \dots, n}.
#' We want to estimate the conditional quantiles of \eqn{X_1} given \eqn{X_3 = x_3}
#' at point \eqn{u_1} using the following kernel-based estimator
#' \deqn{\hat Q(u_1 | X_3 = x_3) := \hat P^{(-1)}(u_1 \leq x_1 | X_3 = x_3),}
#' where
#' \deqn{\hat P(X_1 \leq x_1 | X_3 = x_3)
#' := \frac{\sum_{l=1}^n 1 \{X_(l,1) \leq x_1 \} K_h(X_(l,3) - x_3)}
#' {\sum_{l=1}^n K_h(X_(l,3) - x_3)},}
#' for every \eqn{u_1} in \code{probsX1} and every \eqn{x_3} in \code{newX3}.
#' The \code{matrixK3} should be a matrix of the values \eqn{K_h(X_(l,3) - x_3)}
#' such as the one produced by
#' \code{\link{computeKernelMatrix}(observedX3, newX3, kernel, h)}.
#'
#' @param observedX1 a sample of observations of X1 of size n
#' @param probsX1 a sample of probabilities at which we want to compute
#' the quantiles for the variable X1, of size p1
#' @param matrixK3 a matrix of kernel values of dimension (p2 , n)
#' \eqn{\big(K_h(X3[i] - U3[j])\big)_{i,j}}
#' such as given by \code{\link{computeKernelMatrix}}.
#'
#' @return A matrix of dimensions \code{(p1,p2)} whose (i,j) entry is \eqn{\hat Q(u_1 | X_3 = x_3)}
#' with \eqn{u_1} = \code{probsX1[i]} and \eqn{x_3} = \code{newX3[j]},
#' where \code{newX3[j]} is the vector that was used to construct \code{matrixK3}.
#'
#' @examples
#' Y = MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = cbind(c(1, 0.9), c(0.9, 1)))
#' matrixK = computeKernelMatrix(observedX = Y[,2] , newX = c(0, 1, 2.5),
#'   kernel = "Gaussian", h = 0.8)
#' matrixnp = estimateCondQuantiles(observedX1 = Y[,2],
#'   probsX1 = c(0.3, 0.5) , matrixK3 = matrixK)
#' matrixnp
#'
#' @export
#'
estimateCondQuantiles <- function(observedX1, probsX1, matrixK3)
{
  nNewPoints1 = length(probsX1)
  nNewPoints2 = dim(matrixK3)[2]
  nObserved = length(observedX1)
  stopifnot(nObserved == dim(matrixK3)[1])

  observedX1_order = order(observedX1 , decreasing = FALSE)
  observedX1_sorted = observedX1[observedX1_order]

  matrix_condCDF = matrix(nrow = length(observedX1), ncol = nNewPoints2)
  # For each j, we compute the cumulated kernel conditional to X_3 = X_(j,3)
  for (j in 1:nNewPoints2)
  {
    matrix_condCDF[, j] = cumsum( matrixK3[observedX1_order, j] )
    matrix_condCDF[, j] = matrix_condCDF[,j] / matrix_condCDF[nObserved, j]
  }

  # Computation of the conditional quantiles for the probability probsX1[i]
  # by inversion of the conditional cdf (conditional to X_3 = X_(j,3))
  matrix_condQuantiles = matrix(nrow = nNewPoints1, ncol = nNewPoints2)
  for (i in 1:nNewPoints1)
  {
    for (j in 1:nNewPoints2)
    {
      matrix_condQuantiles[i,j] =
        observedX1_sorted[ min( which( matrix_condCDF[,j] >= probsX1[i] ) ) ]
    }
  }
  return (matrix_condQuantiles)
}

