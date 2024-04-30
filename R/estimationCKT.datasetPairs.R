

#' Construct a dataset of pairs of observations for the estimation
#' of conditional Kendall's tau
#'
#' In (Derumigny, & Fermanian (2019)), it is described how the problem
#' of estimating conditional Kendall's tau can be rewritten as a
#' classification task for a dataset of pairs (of observations).
#' This function computes such a dataset, that can be then used to
#' estimate conditional Kendall's tau using one of the following
#' functions:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.fit.nNets}},
#' \code{\link{CKT.predict.kNN}}.
#'
#' @param X1 vector of observations of the first conditioned variable.
#'
#' @param X2 vector of observations of the second conditioned variable.
#'
#' @param Z vector or matrix of observations of the conditioning variable(s),
#' of dimension \code{dimZ}.
#'
#' @param h the bandwidth. Can be a vector; in this case,
#' the components of \code{h} will be reused to match the dimension of \code{Z}.
#'
#' @param cut the cutting level to keep a given pair or not.
#' Used only if no \code{nPairs} is provided.
#'
#' @param onlyConsecutivePairs if \code{TRUE}, only consecutive pairs are used.
#  else every pair is used.

#' @param nPairs number of most relevant pairs to keep in the final datasets.
#' If this is different than the default \code{NULL},
#' the cutting level \code{cut} is not used.
#'
#' @return A matrix with \code{(4+dimZ)} columns and \code{n*(n-1)/2} rows
#' if \code{onlyConsecutivePairs=FALSE} and else \code{(n/2)} rows.
#' It is structured in the following way:
#' \itemize{
#'     \item column \code{1} contains the information
#'     about the concordance of the pair (i,j) ;
#'
#'     \item columns \code{2} to \code{1+dimZ} contain
#'     the mean value of Z (the conditioning variables) ;
#'
#'     \item column \code{2+dimZ} contains the value of the kernel K_h(Z_j - Z_i) ;
#'
#'     \item column \code{3+dimZ} and \code{4+dimZ}
#'     contain the corresponding values of i and j.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 1 for all pairs and Algorithm 8
#' for the case of only consecutive pairs)
#' \doi{10.1016/j.csda.2019.01.013}
#'
#' @seealso the functions that require such a dataset of pairs
#' to do the estimation of conditional Kendall's tau:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.fit.nNets}},
#' \code{\link{CKT.predict.kNN}}, and \code{\link{CKT.fit.randomForest}}.
#'
#' @examples
#' # We simulate from a conditional copula
#' N = 500
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 0.9 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N = N , family = 3,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau) )
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' datasetP = datasetPairs(
#' X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
#'
#' @export
#'
datasetPairs <- function(X1, X2, Z, h, cut = 0.9,
                         onlyConsecutivePairs = FALSE, nPairs = NULL)
{
  .checkSame_nobs_X1X2Z(X1, X2, Z)

  .checkUnivX1X2(X1, X2)

  n = length(X1)
  dimZ = NCOL(Z)
  if (dimZ == 1 && is.matrix(Z)){
    Z = as.numeric(Z)
  }

  if (onlyConsecutivePairs){

    datasetPairs = matrix(nrow = n/2, ncol = 1 + dimZ + 3)
    hVect = rep(h, length.out = dimZ)
    indexes = 1:(n/2)
    vect_i = 2 * indexes - 1
    vect_j = 2 * indexes
    datasetPairs[indexes, 1] =
      as.numeric( (X1[vect_j]-X1[vect_i]) * (X2[vect_j]-X2[vect_i]) > 0 )
    - as.numeric( (X1[vect_j]-X1[vect_i]) * (X2[vect_j]-X2[vect_i]) < 0 )
    if (dimZ == 1){
      datasetPairs[indexes,2] = (Z[vect_i]+Z[vect_j])/2
      normDiffZ = ( (Z[vect_i]-Z[vect_j]) / h )^2
    } else {
      normDiffZ = 0
      for (k in 1:dimZ)
      {
        datasetPairs[indexes, 1+k] = (Z[vect_i, k]+Z[vect_j, k])/2
        normDiffZ = normDiffZ + ( (Z[vect_i, k]-Z[vect_j, k]) / hVect[k] )^2
      }
    }
    datasetPairs[indexes, 1 + dimZ + 1] = exp( - normDiffZ)
    datasetPairs[indexes, 1 + dimZ + 2] = vect_i
    datasetPairs[indexes, 1 + dimZ + 3] = vect_j

  } else {

    datasetPairs = matrix(nrow = n*(n-1)/2, ncol = 1 + dimZ + 3)
    hVect = rep(h, length.out = dimZ)
    lastPosition = 0
    for (i in 1:(n-1))
    {
      j = (i+1):n
      indexes = lastPosition + 1:(n-i)
      datasetPairs[indexes,1] =
        as.numeric( (X1[j]-X1[i])*(X2[j]-X2[i]) > 0 )
      - as.numeric( (X1[j]-X1[i])*(X2[j]-X2[i]) < 0 )
      if (dimZ == 1){
        datasetPairs[indexes,2] = (Z[i]+Z[j])/2
        normDiffZ = ( (Z[i]-Z[j]) / h )^2
      } else {
        normDiffZ = 0
        for (k in 1:dimZ)
        {
          datasetPairs[indexes,1+k] = (Z[i, k]+Z[j, k])/2
          normDiffZ = normDiffZ + ( (Z[i, k]-Z[j, k]) / hVect[k] )^2
        }
      }
      datasetPairs[indexes, 1 + dimZ + 1] = exp( - normDiffZ)
      datasetPairs[indexes, 1 + dimZ + 2] = i
      datasetPairs[indexes, 1 + dimZ + 3] = j

      lastPosition = lastPosition + (n-i)
    }
  }
  if (!is.null(nPairs)) {
    if (nPairs < nrow(datasetPairs)) {
      cut = stats::quantile(x = datasetPairs[, dimZ+2],
                            probs = nPairs / nrow(datasetPairs))
    } else {
      warning("nPairs larger or equal to the created number of pairs: ",
              nPairs, " , ", nrow(datasetPairs),
              " No cut will be realized.")
    }
  }

  datasetPairs = datasetPairs[ datasetPairs[, dimZ+2]>=cut ,]
  if (is.null(colnames(Z))){
    colnames_Z = paste(rep("Z",dimZ), 1:dimZ, sep="_")
  } else {
    colnames_Z = colnames(Z)
  }
  colnames(datasetPairs) <- c("concordant", colnames_Z,
                              "kernel.value", "iUsed", "jUsed")

  return (datasetPairs)
}


#' Construct a dataset of pairs to be used for the cross-validation
#' choice of the bandwidth for kernel estimation of conditional Kendall's tau
#'
#' @param observedX1 a vector of n observations of the first variable.
#' @param observedX2 a vector of n observations of the second variable.
#' @param nPairs number of most relevant pairs to keep in the final datasets.
#' If \code{NULL}, no cut is realized.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#'
#' @return the corresponding dataset of pairs
#'
#' @noRd
#'
datasetPairs_hCV <- function(X1, X2, Z, nPairs = NULL, typeEstCKT = 2)
{
  n = length(X1)
  if (length(X1) != length(X2)) {
    stop(paste0("X1 and X2 have different lengths: ", length(X1), " and ", length(X2)))
  }
  if (is.vector(Z)) {
    dimZ = 1
    if (length(X1) != length(Z)) {
      stop(paste0("X1 and Z have different lengths: ", length(X1), " and ", length(Z)))
    }
  } else {
    if (length(X1) != length(Z[,1])) {
      stop(paste0("X1 and Z have different lengths: ", length(X1), " and ", length(Z)))
    }
    dimZ = length(Z[1,])
  }

  # matrix_diff = outer(
  #   X = vectorZ, Y = vectorZ, FUN = function (u1, u2) {
  #     return(abs(u1 - u2) ) } )
  # diag(matrix_diff) <- NA
  #
  # if (2*nPairs < length(matrix_diff)-n){
  #   cut = stats::quantile(x = as.numeric(matrix_diff),
  #                  probs = (2*nPairs / (length(matrix_diff)-n)) , na.rm = TRUE)
  # } else { cut = +Inf }


  dataMatrix = matrix(nrow = n*(n-1)/2, ncol = 1 + dimZ + 3)
  lastPosition = 0

  if (typeEstCKT == 4 || typeEstCKT == "wdm"){typeEstCKT <- 2} # Same function for types 2 and 4

  switch(
    typeEstCKT,

    # 1
    {for (i in 1:(n-1)){
      j = (i+1):n
      indexes = lastPosition + 1:(n-i)

      dataMatrix[indexes,1] = 2 * as.numeric( (X1[j]-X1[i])*(X2[j]-X2[i]) > 0 ) - 1

      if (dimZ == 1){
        dataMatrix[indexes,2] = (Z[i]+Z[j])/2
        normDiffZ = ( (Z[i]-Z[j]) )^2
      } else {
        normDiffZ = 0
        for (k in 1:dimZ)
        {
          dataMatrix[indexes,1+k] = (Z[i, k]+Z[j, k])/2
          normDiffZ = normDiffZ + ( (Z[i, k]-Z[j, k])  )^2
        }
      }
      dataMatrix[indexes, 1 + dimZ + 1] = sqrt(normDiffZ)
      dataMatrix[indexes, 1 + dimZ + 2] = i
      dataMatrix[indexes, 1 + dimZ + 3] = j

      lastPosition = lastPosition + (n-i)
    } },

    # 2
    {for (i in 1:(n-1)){
      j = (i+1):n
      indexes = lastPosition + 1:(n-i)

      dataMatrix[indexes,1] = as.numeric( (X1[j]-X1[i])*(X2[j]-X2[i]) > 0 )
      - as.numeric( (X1[j]-X1[i])*(X2[j]-X2[i]) < 0 )

      if (dimZ == 1){
        dataMatrix[indexes,2] = (Z[i]+Z[j])/2
        normDiffZ = ( (Z[i]-Z[j]) )^2
      } else {
        normDiffZ = 0
        for (k in 1:dimZ)
        {
          dataMatrix[indexes,1+k] = (Z[i, k]+Z[j, k])/2
          normDiffZ = normDiffZ + ( (Z[i, k]-Z[j, k])  )^2
        }
      }
      dataMatrix[indexes, 1 + dimZ + 1] = sqrt(normDiffZ)
      dataMatrix[indexes, 1 + dimZ + 2] = i
      dataMatrix[indexes, 1 + dimZ + 3] = j

      lastPosition = lastPosition + (n-i)
    } },

    # 3
    {for (i in 1:(n-1)){
      j = (i+1):n
      indexes = lastPosition + 1:(n-i)

      dataMatrix[indexes,1] = 1 - 2 * as.numeric( (X1[j]-X1[i])*(X2[j]-X2[i]) < 0 )

      if (dimZ == 1){
        dataMatrix[indexes,2] = (Z[i]+Z[j])/2
        normDiffZ = ( (Z[i]-Z[j]) )^2
      } else {
        normDiffZ = 0
        for (k in 1:dimZ)
        {
          dataMatrix[indexes,1+k] = (Z[i, k]+Z[j, k])/2
          normDiffZ = normDiffZ + ( (Z[i, k]-Z[j, k])  )^2
        }
      }
      dataMatrix[indexes, 1 + dimZ + 1] = sqrt(normDiffZ)
      dataMatrix[indexes, 1 + dimZ + 2] = i
      dataMatrix[indexes, 1 + dimZ + 3] = j

      lastPosition = lastPosition + (n-i)
    } }
  )

  # if (nPairs < nrow(dataMatrix)){
  #   cut = stats::quantile(
  #     x = as.numeric(dataMatrix[indexes, 1 + dimZ + 1]) ,
  #     probs = nPairs / dataMatrix[indexes, 1 + dimZ + 1] , na.rm = TRUE)
  # } else { cut = +Inf }
  if (!is.null(nPairs) && nPairs < nrow(dataMatrix)){
    cut = stats::quantile(x = dataMatrix[, dimZ+2], probs = nPairs / nrow(dataMatrix))
  }

  dataMatrix = dataMatrix[ dataMatrix[,dimZ+2] <= cut ,]
  colnames(dataMatrix) <- c("concordant", paste(rep("Z",dimZ), 1:dimZ, sep="_"),
                            "normDiff", "iUsed", "jUsed")
  return (dataMatrix)
}


#' Compute the matrix of signs of pairs
#'
#' Compute a matrix giving the concordance or discordance
#' of each pair of observations.
#'
#' @param vectorX1 vector of observed data (first coordinate)
#'
#' @param vectorX2 vector of observed data (second coordinate)
#'
#' @param typeEstCKT if typeEstCKT = 2 or 4, compute the matrix whose term (i,j) is :
#'    \deqn{1 \{ (X_{i,1} - X_{j,1}) * (X_{i,2} - X_{j,2}) > 0 \}
#'  - 1 \{ (X_{i,1} - X_{j,1}) * (X_{i,2} - X_{j,2}) < 0 \},}
#'  where \eqn{1} is the indicator function.
#'
#' For \code{typeEstCKT = 1} (respectively \code{typeEstCKT = 3})
#' a negatively biased (respectively positively) matrix is given.
#'
#' @return an \code{n * n} matrix with the signs of each pair
#' of observations.
#'
#' @examples
#' # We simulate from a conditional copula
#' N = 500
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 0.9 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N = N , family = 3,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau) )
#' matrixPairs = computeMatrixSignPairs(vectorX1 = simCopula[,1],
#'                                      vectorX2 = simCopula[,2])
#'
#' @export
#'
computeMatrixSignPairs <- function(vectorX1, vectorX2, typeEstCKT = 4)
{
  if (length(vectorX1) != length(vectorX2)) {
    stop(errorCondition(
      message = paste0("vectorX1 and vectorX2 have different lengths: ",
                       length(vectorX1), " and ", length(vectorX2)),
      class = "DifferentLengthsError") )
  }

  n = length(vectorX1)
  matrixSignsPairs = matrix(nrow = n, ncol = n)
  if (typeEstCKT == 1)
  {
    for (i in (1:(n-1)))
    {
      matrixSignsPairs[i,i] = 0
      matrixSignsPairs[i,(i+1):n] = as.numeric(
        vectorX1[i] < vectorX1[(i+1):n] & vectorX2[i] < vectorX2[(i+1):n] )
      matrixSignsPairs[(i+1):n,i] = as.numeric(
        vectorX1[(i+1):n] < vectorX1[i] & vectorX2[(i+1):n] < vectorX2[i] )
    }
  } else if (typeEstCKT %in% c(2,4) ) {
    for (i in (1:(n-1)))
    {
      matrixSignsPairs[i,i] = 0
      signVector = (vectorX1[i] - vectorX1[(i+1):n]) *
        (vectorX2[i] - vectorX2[(i+1):n])
      matrixSignsPairs[i,(i+1):n] =
        as.numeric( signVector > 0) - as.numeric( signVector < 0)
      matrixSignsPairs[(i+1):n,i] = matrixSignsPairs[i,(i+1):n]

    }
  } else if (typeEstCKT == 3) {
    for (i in (1:(n-1)))
    {
      matrixSignsPairs[i,i] = 0
      matrixSignsPairs[i,(i+1):n] = as.numeric(
        vectorX1[i] < vectorX1[(i+1):n] & vectorX2[i] > vectorX2[(i+1):n] )
      matrixSignsPairs[(i+1):n,i] = as.numeric(
        vectorX1[(i+1):n] < vectorX1[i] & vectorX2[(i+1):n] > vectorX2[i] )
    }
  } else {stop(paste0("typeEstCKT: ", typeEstCKT, " is not in {1,2,3,4}" ) ) }

  matrixSignsPairs[n,n] = 0
  return (matrixSignsPairs)
}
