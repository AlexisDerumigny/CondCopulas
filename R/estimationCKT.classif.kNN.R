


#' Prediction of conditional Kendall's tau using nearest neighbors
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
#' In other words, conditional Kendall's tau is the difference
#' between the probabilities of observing concordant and discordant pairs
#' from the conditional law of \deqn{(X_1, X_2) | Z=z.}
#' This function estimates conditional Kendall's tau using a
#' \strong{nearest neighbors}. This is possible by the relationship between
#' estimation of conditional Kendall's tau and classification problems
#' (see Derumigny and Fermanian (2019)): estimation of conditional Kendall's tau
#' is equivalent to the prediction of concordance in the space of pairs
#' of observations.
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#'
#' @param designMatrix the matrix of predictors.
#' They must have the same number of variables as \code{newZ} and
#' the same number of observations as \code{inputMatrix},
#' i.e. there should be one "multivariate observation" of the predictor for each pair.
#'
#' @param newZ the matrix of predictors for which we want to
#' estimate the conditional Kendall's taus at these values.
#'
#' @param number_nn vector of numbers of nearest neighbors to use.
#' If several number of neighbors are given (local) aggregation is performed
#' using Lepski's method on the subset determined by the \code{partition}.
#'
#' @param weightsVariables optional argument to give
#' different weights \eqn{w_j} to each variable.
#'
#' @param normLp the p in the weighted p-norm \eqn{|| x ||_p = \sum_j w_j * x_j^p}
#' used to determine the distance in the computation of the nearest neighbors.
#'
#' @param constantA a tuning parameter that controls the adaptation.
#' The higher, the smoother it is; while the smaller, the least smooth it is.
#'
#' @param partition used only if \code{length(number_nn) > 1}.
#' It is the number of subsets to consider for the local choice of the
#' number of nearest neighbors ; or a vector giving the id of each observations
#' among the subsets.
#' If \code{NULL}, only one set is used.
#'
#' @param verbose if TRUE, this print information each \code{lengthVerbose} iterations
#' @param lengthVerbose number of iterations at each time for which progress is printed.
#'
#' @param methodSort is the sorting method used to find the nearest neighbors.
#' Possible choices are
#' \code{ecdf} (uses the ecdf to order the points to find the neighbors)
#' and \code{partial.sort} uses a partial sorting algorithm.
#' This parameter should not matter except for the computation time.
#'
#' @return  a list with two components
#' \itemize{
#'     \item \code{estimatedCKT} the estimated conditional Kendall's tau, a vector of
#'     the same size as the number of rows in \code{newZ};
#'     \item \code{vect_k_chosen} the locally selected number of nearest neighbors,
#'     a vector of the same size as the number of rows in \code{newZ}.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 5)
#'
#' @seealso See also other estimators of conditional Kendall's tau:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.nNets}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.kernel}},
#' \code{\link{CKT.kendallReg.fit}},
#' and the more general wrapper \code{\link{CKT.estimate}}.
#'
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
#' datasetP = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
#' estimatedCKT_knn <- CKT.predict.kNN(
#'   datasetPairs = datasetP,
#'   newZ = matrix(newZ,ncol = 1),
#'   number_nn = c(50,80, 100, 120,200),
#'   partition = 8)
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in red)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col="black",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_knn$estimatedCKT, col = "red")
#'
#' @export
#'
CKT.predict.kNN <- function(
  datasetPairs,
  designMatrix = datasetPairs[,2:(ncol(datasetPairs)-3),drop=FALSE],
  newZ,
  number_nn, weightsVariables = 1, normLp = 2,
  constantA = 1, partition = NULL,
  verbose = 1, lengthVerbose = 100,
  methodSort = "partial.sort")
{
  if (length(number_nn) == 1){

    estimatedCKT = CKT.predict.kNN.1(
      datasetPairs = datasetPairs, designMatrix = designMatrix,
      newZ = newZ, number_nn = number_nn,
      weightsVariables = weightsVariables, normLp = normLp,
      verbose = verbose-1, lengthVerbose = lengthVerbose, methodSort = methodSort)

    return (list(estimatedCKT = estimatedCKT))
  }

  matrixKNN = matrix(nrow = NROW(newZ), ncol = length(number_nn))
  for (i_nn in 1:length(number_nn)){
    matrixKNN[, i_nn] = CKT.predict.kNN.1(
      datasetPairs = datasetPairs, designMatrix = designMatrix,
      newZ = newZ, number_nn = number_nn[i_nn],
      weightsVariables = weightsVariables, normLp = normLp,
      verbose = verbose, lengthVerbose = lengthVerbose, methodSort = methodSort)
  }

  result = CKT.adaptkNN(
    matrixKNN = matrixKNN, vect_k = number_nn,
    constantA = constantA, partition = partition, verbose = verbose-1)

  return (result)
}


#' Prediction of conditional Kendall's tau using nearest neighbors
#' (for one number of nearest neighbors)
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#'
#' @param designMatrix the matrix of predictors that will be used for the nearest neighbors
#' They must have the same number of variables as \code{newZ} and
#' the same number of observations as \code{inputMatrix},
#' i.e. there should be one "multivariate observation" of the predictor for each pair.
#'
#' @param newZ the matrix of predictors for which we want to
#' estimate the conditional Kendall's taus at these values.
#'
#' @param number_nn number of nearest neighbors to use.
#'
#' @param weightsVariables optional argument to give
#' different weights \eqn{w_j} to each variable.
#'
#' @param normLp the p in the weighted p-norm \eqn{|| x ||_p = \sum_j w_j * x_j^p}
#' used to determine the distance in the computation of the nearest neighbors.
#'
#' @param verbose if 0, this print information each \code{lengthVerbose} iterations
#' @param lengthVerbose number of iterations at each time for which progress is printed.
#'
#' @param methodSort is the sorting method used to find the nearest neighbors.
#' Possible choices are
#' \code{ecdf} (uses the ecdf to order the points to find the neighbors)
#' and \code{partial.sort} uses a partial sorting algorithm.
#' This parameter should not matter except for the computation time.
#'
#' @return a vector of size \code{nrow(newZ)} of estimated conditional Kendall's taus
#' corresponding to each (multivariate) predictor.
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 5)
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
#' datasetP = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
#' estimatedCKT_knn <- CKT.predict.kNN.1(
#'   datasetPairs = datasetP,
#'   newZ = matrix(newZ,ncol = 1),
#'   number_nn = 100)
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in red)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col="black",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_knn, col = "red")
#'
#'
#'
#' @noRd
#'
CKT.predict.kNN.1 <- function(datasetPairs,
                              designMatrix,
                              newZ,
                              number_nn, weightsVariables = 1, normLp = 2,
                              verbose = 1, lengthVerbose = 100,
                              methodSort = "partial.sort")
{
  n_data = nrow(designMatrix)
  n_newZ = nrow(newZ)
  pPrime = ncol(designMatrix)

  if (n_data != nrow(datasetPairs)){
    stop("designMatrix and datasetPairs should have the same number of rows.")
  } else if (pPrime != ncol(designMatrix)){
    stop("designMatrix and newZ should have the same number of columns.")
  }

  weightsVar = rep(weightsVariables, length.out = pPrime)
  prediction = rep(NA, n_newZ)

  prt = proc.time()
  for (i_data in 1:n_newZ)
  {
    pointToEstimate = newZ[i_data, ]
    distance = rep(0, n_data)
    for (i_var in 1:pPrime)
    {
      distance = distance + weightsVar[i_var] *
        (abs(designMatrix[,i_var] - pointToEstimate[i_var]))^normLp
    }
    # order_dist = order(distance)[1 : number_nn]
    if (methodSort == "partial.sort")
    {
      vect_ok_dist = which(distance <= sort(distance, partial=number_nn)[number_nn])
    } else if (methodSort == "ecdf"){
      vect_ok_dist = which(stats::ecdf(distance)(distance) <= number_nn / n_data)
    } else {
      stop(paste0("Method ", methodSort, " is not implemented yet."))
    }

    prediction[i_data] = stats::weighted.mean(
      x = datasetPairs[vect_ok_dist , "concordant"]
      , w = datasetPairs[vect_ok_dist , "kernel.value"] )

    if ((verbose>0) & i_data%%lengthVerbose == 0)
    {
      prt2 = proc.time()
      cat(i_data) ; cat(" out of ") ; cat(n_newZ)
      cat("  ") ; cat((prt2-prt)[3]) ; cat("s \n")
      prt = proc.time()
    }
  }

  predict_CKT = 2 * prediction - 1

  return (predict_CKT)
}


#' Local aggregation of nearest neighbor estimators
#' of the conditional Kendall's tau by Lepski's method.
#'
#' @param matrixKNN the matrix of nearest neighbor-estimated conditional Kendall's taus.
#' Each column should correspond to a choice of number of nearest neighbors
#' and each row should correspond to the conditional Kendall's taus
#' estimated at the new point.
#'
#' @param vect_k vector of the numbers of nearest neighbors
#' (corresponding to each column of the matrixKNN).
#'
#' @param constantA a tuning parameter that controls the adaptation.
#' The higher, the smoother it is; while the smaller, the least smooth it is.
#'
#' @param partition the number of subsets to consider;
#' or a vector giving the id of each observations.
#' If \code{NULL}, only one set is used.
#'
#' @param verbose if larger than 0, prints the evolution of the distances
#' among the path of numbers of nearest neighbors.
#'
#' @return a list with two components
#' \itemize{
#'     \item \code{estimatedCKT} the estimated conditional Kendall's tau, a vector of
#'     the same size as the number of rows in \code{matrixKNN};
#'     \item \code{vect_k_chosen} the locally selected number of nearest neighbors,
#'     a vector of the same size as the number of rows in \code{matrixKNN}.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 6)
#'
#' @noRd
#'
CKT.adaptkNN <- function(matrixKNN, vect_k,
                         constantA = 1, partition = NULL, verbose = 0)
{
  CKT_esti_min = min(matrixKNN)
  CKT_esti_max = max(matrixKNN)

  n_grid_z = nrow(matrixKNN)
  number_k = length(vect_k)

  if (ncol(matrixKNN) != number_k) {
    stop("Lengths are different for matrixKNN[1,] and vect_k")
  }

  if (is.null(partition)){
    partition_id = rep(1, n_grid_z)
  } else if (length(partition) == 1){
    partition_id = ceiling((1:n_grid_z)*(partition/n_grid_z))
  } else {
    partition_id = partition
  }
  sizePartition = max(partition)

  order_k = order(vect_k)

  vect_k_chosen = rep(NA, n_grid_z)

  for (iPartition in 1:sizePartition)
  {
    indexes = which(partition_id == iPartition)
    iChosen = 1
    for (i_k in number_k:2)
    {
      accepted = TRUE
      for (j_k in (i_k-1):1)
      {
        vectDiff = abs(matrixKNN[indexes, order_k[i_k]] - matrixKNN[indexes, order_k[j_k]])
        vectDiffRenorm = vectDiff / (CKT_esti_max - CKT_esti_min)
        distance_ij = sqrt(mean(vectDiff^2))
        limitAccept = constantA / ( sqrt(vect_k[order_k[j_k]]) ) *
          sqrt(log( vect_k[order_k[number_k]] / vect_k[order_k[j_k]] ))
        if (distance_ij > limitAccept)
        {
          if (verbose>0) {
            cat(i_k) ; cat("  ") ; cat(order_k[i_k]) ; cat("  ")
            cat(prettyNum(distance_ij)); cat("  ")
            cat(prettyNum(limitAccept)); cat("\n")
          }
          accepted = FALSE ; break
        }
      }
      if (accepted){
        iChosen = i_k ; break
      }
    }
    vect_k_chosen[indexes] = order_k[iChosen]
  }
  estimator = diag(matrixKNN[, vect_k_chosen])
  return (list(estimatedCKT = estimator, vect_k_chosen = vect_k_chosen))
}


