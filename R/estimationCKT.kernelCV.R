
#' Choose the bandwidth for kernel estimation of
#' conditional Kendall's tau using cross-validation
#'
#'
#' @description
#' Let \eqn{X_1} and \eqn{X_2} be two random variables.
#' The goal here is to estimate the conditional Kendall's tau
#' (a dependence measure) between \eqn{X_1} and \eqn{X_2} given \eqn{Z=z}
#' for a conditioning variable \eqn{Z}.
#' Conditional Kendall's tau between \eqn{X_1} and \eqn{X_2} given \eqn{Z=z}
#' is defined as:
#' \deqn{P( (X_{1,1} - X_{2,1})(X_{1,2} - X_{2,2}) > 0 | Z_1 = Z_2 = z)}
#' \deqn{- P( (X_{1,1} - X_{2,1})(X_{1,2} - X_{2,2}) < 0 | Z_1 = Z_2 = z),}
#' where \eqn{(X_{1,1}, X_{1,2}, Z_1)} and \eqn{(X_{2,1}, X_{2,2}, Z_2)}
#' are two independent and identically distributed copies of \eqn{(X_1, X_2, Z)}.
#' For this, a kernel-based estimator is used, as described in
#' (Derumigny & Fermanian (2019)).
#' These functions aims at finding the best bandwidth \code{h} among a given
#' \code{range_h} by cross-validation. They use either:
#' \itemize{
#'    \item \strong{leave-one-out} cross-validation:
#'    function \code{CKT.hCV.l1out}
#'
#'    \item or \strong{K-folds} cross-validation:
#'    function \code{CKT.hCV.Kfolds}
#' }
#'
#'
#' @param observedX1 a vector of \code{n} observations of the first variable
#'
#' @param observedX2 a vector of \code{n} observations of the second variable
#'
#' @param observedZ observedZ vector of observed values of Z.
#' If Z is multivariate, then this is a matrix whose rows correspond
#' to the observations of Z
#'
#' @param range_h vector containing possible values for the bandwidth.
#'
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by \code{\link{computeMatrixSignPairs}(observedX1, observedX2)}.
#' Only needed if \code{typeEstCKT} is not the default 'wdm'.
#'
#' @param nPairs number of pairs used in the cross-validation criteria.
#'
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#'
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are \code{"Gaussian"} (Gaussian kernel)
#' and \code{"Epa"} (Epanechnikov kernel).
#'
#' @param progressBar if \code{TRUE}, a progressbar for each h is displayed
#' to show the progress of the computation.
#'
#' @param verbose if \code{TRUE}, print the score of each h during the procedure.
#'
#'
#' @return Both functions return a list with two components:
#' \itemize{
#'     \item \code{hCV}: the chosen bandwidth
#'     \item \code{scores}: vector of the same length as range_h giving the
#'     value of the CV criteria for each of the h tested.
#'     Lower score indicates a better fit.
#' }
#'
#' @seealso \code{\link{CKT.kernel}} for the corresponding
#' estimator of conditional Kendall's tau by kernel smoothing.
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' On kernel-based estimation of conditional Kendall’s tau:
#' finite-distance bounds and asymptotic behavior.
#' Dependence Modeling, 7(1), 292-321.
#' Page 296, Equation (4).
#' \doi{10.1515/demo-2019-0016}
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 200
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' newZ = seq(2,10,by = 0.1)
#' range_h = 3:10
#'
#' resultCV <- CKT.hCV.l1out(observedX1 = X1, observedX2 = X2,
#'   range_h = range_h, observedZ = Z, nPairs = 100)
#'
#' resultCV <- CKT.hCV.Kfolds(observedX1 = X1, observedX2 = X2,
#'   range_h = range_h, observedZ = Z, ZToEstimate = newZ)
#'
#' plot(range_h, resultCV$scores, type = "b")
#'
#' @export
#'
CKT.hCV.l1out <- function (observedX1, observedX2, observedZ,
                           range_h, matrixSignsPairs = NULL,
                           nPairs = 10*length(observedX1),
                           typeEstCKT = "wdm", kernel.name = "Epa",
                           progressBar = TRUE, verbose = FALSE)
{
  n_h = length(range_h)
  scores = rep(NA, n_h)

  .checkSame_nobs_X1X2Z(observedX1, observedX2, observedZ)
  .checkUnivX1X2(observedX1, observedX2)

  # Construction of the matrix of selected pairs
  dataMatrix = datasetPairs_hCV(
    X1 = observedX1, X2 = observedX2, Z = observedZ,
    nPairs = nPairs, typeEstCKT = typeEstCKT)

  nPairsReal = nrow(dataMatrix)
  if(nPairsReal != nPairs) {
    warning(paste0("nPairs =", nPairs, " but real number of pairs is ", nPairsReal))
  }

  if (is.vector(observedZ)){
    computeScore <- function(i, i_h, h){
      toBeRemoved = - c(dataMatrix[i,4], dataMatrix[i,5])
      matrixScore[i, i_h] <<- CKT.kernelPointwise.univariate(
        X1 = observedX1[toBeRemoved], X2 = observedX2[toBeRemoved],
        matrixSignsPairs = matrixSignsPairs[toBeRemoved , toBeRemoved],
        vectorZ = observedZ[toBeRemoved],
        h = h, pointZ = dataMatrix[i,2],
        kernel.name = kernel.name, typeEstCKT = typeEstCKT)
    }
  } else {
    dimZ = ncol(observedZ)
    computeScore <- function(i, i_h, h){
      toBeRemoved = - c(dataMatrix[i,4], dataMatrix[i,5])
      matrixScore[i, i_h] <<- CKT.kernelPointwise.multivariate(
        X1 = observedX1[toBeRemoved], X2 = observedX2[toBeRemoved],
        matrixSignsPairs = matrixSignsPairs[toBeRemoved , toBeRemoved],
        matrixZ = observedZ[toBeRemoved, ],
        h = h, pointZ = dataMatrix[i,1+1:dimZ],
        kernel.name = kernel.name, typeEstCKT = typeEstCKT)
    }
  }
  matrixScore = matrix(nrow = nPairsReal, ncol = n_h)
  for (i_h in 1:n_h)
  {
    h = range_h[i_h]
    if (progressBar){
      pbapply::pbapply(X = array(1:nPairsReal),
                       FUN = computeScore, MARGIN = 1, i_h = i_h, h = h)
    } else {
      apply(X = array(1:nPairsReal),
            FUN = computeScore, MARGIN = 1, i_h = i_h, h = h)
    }

    scores[i_h] = mean( (dataMatrix[,1] - matrixScore[, i_h])^2, na.rm = TRUE )
    if (verbose){print(paste0("h = ", h, " ; score = ", scores[i_h]))}
  }
  hCV = range_h[which.min(scores)]
  return(list(hCV = hCV, scores = scores))
}


#' Choose the bandwidth for estimation of conditional Kendall's tau
#' using K-fold cross-validation
#'
#'
#' @param ZToEstimate vector of fixed conditioning values at which
#' the difference between the two conditional Kendall's tau should be computed.
#' Can also be a matrix whose lines are the conditioning vectors at which
#' the difference between the two conditional Kendall's tau should be computed.
#'
#' @param Kfolds number of subsamples used.
#'
#' @rdname CKT.hCV.l1out
#' @export
#'
CKT.hCV.Kfolds <- function(observedX1, observedX2, observedZ, ZToEstimate,
                           range_h, matrixSignsPairs = NULL,
                           typeEstCKT = "wdm", kernel.name = "Epa",
                           Kfolds = 5, progressBar = TRUE, verbose = FALSE)
{
  n_h = length(range_h)
  scores = rep(NA, n_h)
  n = NROW(observedZ)

  .checkSame_nobs_X1X2Z(observedX1, observedX2, observedZ)
  .checkUnivX1X2(observedX1, observedX2)

  if (is.vector(observedZ)){
    estimCKTNP <- CKT.kernel.univariate
    observedZ = matrix(observedZ, ncol = 1)
  } else {
    estimCKTNP <- CKT.kernel.multivariate
  }

  computeScore <- function(i_h) {
    foldid = sample(rep(seq(Kfolds), length = n))
    list_vectorEstimate = as.list(seq(Kfolds))
    list_vectorEstimate_comp = as.list(seq(Kfolds))
    list_resultEstimation = as.list(seq(Kfolds))

    distanceTot = 0
    h = range_h[i_h]

    for (i in seq(Kfolds)) {
      which = foldid == i
      list_vectorEstimate[[i]] = estimCKTNP(
        X1 = observedX1[!which], X2 = observedX2[!which],
        matrixSignsPairs = matrixSignsPairs[!which, !which],
        Z = observedZ[!which,], ZToEstimate = ZToEstimate,
        typeEstCKT = typeEstCKT, h = h, kernel.name = kernel.name, progressBar = FALSE)

      list_vectorEstimate_comp[[i]] = estimCKTNP(
        X1 = observedX1[which], X2 = observedX2[which],
        matrixSignsPairs = matrixSignsPairs[which, which],
        Z = observedZ[which,], ZToEstimate = ZToEstimate,
        typeEstCKT = typeEstCKT, h = h, kernel.name = kernel.name, progressBar = FALSE)

      if (all(!is.finite(list_vectorEstimate[[i]])))
      {
        distanceTot = NA ; break
      }
      distance = sum((list_vectorEstimate[[i]] - list_vectorEstimate_comp[[i]])^2)
      if (verbose) {
        print(distance)
      }
      distanceTot = distanceTot + distance
    }
    scores[i_h] <<- distanceTot

    return()
  }

  if (progressBar){
    pbapply::pbapply(X = array(1:n_h), FUN = computeScore, MARGIN = 1 )
  } else {
    apply(X = array(1:n_h), FUN = computeScore, MARGIN = 1 )
  }

  hCV = range_h[which.min(scores)]
  return(list(hCV = hCV, scores = scores))
}





# Old code --------------------------------------------------------



# # This function takes in parameter
# # range_h the set of all h tested.
# # It returns a list where hCV is the best h
# # scores is the vector of the same length as range_h giving the
# # CV criteria for each of the h tested
# choose_hCV <- function (range_h, matrixSignsPairs, vectorZ, typeEstCKT, kernel.name)
# {
#   n_h = length(range_h)
#   scores = rep(NA, n_h)
#
#   if (typeEstCKT == 1){
#     matrix_g = 4 * matrixSignsPairs - 1
#   } else if (typeEstCKT == 2 || typeEstCKT == 4){
#     matrix_g = matrixSignsPairs
#   } else if (typeEstCKT == 3){
#     matrix_g = 1 - 4 * matrixSignsPairs
#   } else {stop(paste0("typeEstCKT: ", typeEstCKT, " is not in {1,2,3}" ) ) }
#
#   for (i_h in 1:n_h)
#   {
#     h = range_h[i_h]
#
#     matrix_weights = outer(
#       X = vectorZ, Y = vectorZ, FUN = function (u1, u2) {
#         delta_u = (u1 - u2)/h
#         return(0.75 * (1 - delta_u^2) * as.numeric(abs(delta_u) <= 1) / h ) } )
#
#     matrix_CKT_esti = matrix(nrow = n, ncol = n)
#
#     pbapply(X = array(n:1), FUN = function(i)
#     {
#       if (i == 1) {
#         matrix_CKT_esti[i,i] <<-
#           computeEstimateNP_tau_pointwise(
#             matrixSignsPairs = matrixSignsPairs[-c(i), -c(i)], vectorZ = vectorZ[-c(i)],
#             h = h, pointZ = vectorZ[i], kernel.name = kernel.name,
#             typeEstCKT = typeEstCKT)
#       } else {
#         for (j in 1:(i-1))
#         {
#           matrix_CKT_esti[i,j] <<-
#             computeEstimateNP_tau_pointwise(
#               matrixSignsPairs = matrixSignsPairs[-c(i,j), -c(i,j)], vectorZ = vectorZ[-c(i,j)],
#               h = h, pointZ = (vectorZ[i] + vectorZ[j])/2, kernel.name = kernel.name,
#               typeEstCKT = typeEstCKT)
#           matrix_CKT_esti[j,i] <<- matrix_CKT_esti[i,j]
#         }
#         matrix_CKT_esti[i,i] <<-
#           computeEstimateNP_tau_pointwise(
#             matrixSignsPairs = matrixSignsPairs[-c(i), -c(i)], vectorZ = vectorZ[-c(i)],
#             h = h, pointZ = vectorZ[i], kernel.name = kernel.name,
#             typeEstCKT = typeEstCKT)
#       }
#       return() }, MARGIN = 1 )
#
#     scores[i_h] = sum( (matrix_g - matrix_CKT_esti)^2 * matrix_weights ) / length(matrix_g)
#   }
#   hCV = range_h[which.min(scores)]
#   return(list(hCV = hCV, scores = scores))
# }

# # This function takes in parameter
# # range_h the set of all h tested.
# # It returns a list where hCV is the best h
# # scores is the vector of the same length as range_h giving the
# # CV criteria for each of the h tested
# choose_hCV_indK <- function (liste_h, matrixSignsPairs, vectorZ, typeEstCKT, kernel.name, nPairs)
# {
#   n_h = length(liste_h)
#   scores = rep(NA, n_h)
#
#   if (typeEstCKT == 1){
#     matrix_g = 4 * matrixSignsPairs - 1
#   } else if (typeEstCKT == 2 || typeEstCKT == 4){
#     matrix_g = matrixSignsPairs
#   } else if (typeEstCKT == 3){
#     matrix_g = 1 - 4 * matrixSignsPairs
#   } else {stop(paste0("typeEstCKT: ", typeEstCKT, " is not in {1,2,3}" ) ) }
#
#   matrix_diff = outer(
#     X = vectorZ, Y = vectorZ, FUN = function (u1, u2) {
#       return(abs(u1 - u2) ) } )
#
#   seuil_a = quantile(x = as.numeric(matrix_diff), probs = (2*nPairs / length(matrix_diff)) )
#   matrix_weights = as.numeric(matrix_diff <= seuil_a)
#
#   # matrix_weights = outer(
#   #   X = vectorZ, Y = vectorZ, FUN = function (u1, u2) {
#   #     delta_u = (u1 - u2)/h
#   #     return(0.75 * (1 - delta_u^2) * as.numeric(abs(delta_u) <= 1) / h ) } )
#
#   for (i_h in 1:n_h)
#   {
#     h = liste_h[i_h]
#
#     matrix_CKT_esti = matrix(nrow = n, ncol = n)
#
#     pbapply(X = array(n:1), FUN = function(i)
#     {
#       if (i == 1) {
#         matrix_CKT_esti[i,i] <<-
#           computeEstimateNP_tau_pointwise(
#             matrixSignsPairs = matrixSignsPairs[-c(i), -c(i)], vectorZ = vectorZ[-c(i)],
#             h = h, pointZ = vectorZ[i], kernel.name = kernel.name,
#             typeEstCKT = typeEstCKT)
#       } else {
#         for (j in 1:(i-1))
#         {
#           matrix_CKT_esti[i,j] <<-
#             computeEstimateNP_tau_pointwise(
#               matrixSignsPairs = matrixSignsPairs[-c(i,j), -c(i,j)], vectorZ = vectorZ[-c(i,j)],
#               h = h, pointZ = (vectorZ[i] + vectorZ[j])/2, kernel.name = kernel.name,
#               typeEstCKT = typeEstCKT)
#           matrix_CKT_esti[j,i] <<- matrix_CKT_esti[i,j]
#         }
#         matrix_CKT_esti[i,i] <<-
#           computeEstimateNP_tau_pointwise(
#             matrixSignsPairs = matrixSignsPairs[-c(i), -c(i)], vectorZ = vectorZ[-c(i)],
#             h = h, pointZ = vectorZ[i], kernel.name = kernel.name,
#             typeEstCKT = typeEstCKT)
#       }
#       return() }, MARGIN = 1 )
#
#     scores[i_h] = mean( (matrix_g - matrix_CKT_esti)^2 * matrix_weights, na.rm = TRUE )
#   }
#   hCV = liste_h[which.min(scores)]
#   return(list(hCV = hCV, scores = scores))
# }


# # This function takes in parameter
# # range_h the set of all h tested.
# # It returns a list where hCV is the best h
# # scores is the vector of the same length as range_h giving the
# # CV criteria for each of the h tested
# choose_hCV_indK_fast <- function (
#   observedX1, observedX2, observedZ,
#   liste_h, matrixSignsPairs, vectorZ, typeEstCKT,
#                                   kernel.name, nPairs, verbose = TRUE)
# {
#   n_h = length(liste_h)
#   scores = rep(NA, n_h)
#
#   # if (typeEstCKT == 1){
#   #   matrix_g = 4 * matrixSignsPairs - 1
#   # } else if (typeEstCKT == 2 || typeEstCKT == 4){
#   #   matrix_g = matrixSignsPairs
#   # } else if (typeEstCKT == 3){
#   #   matrix_g = 1 - 4 * matrixSignsPairs
#   # } else {stop(paste0("typeEstCKT: ", typeEstCKT, " is not in {1,2,3}" ) ) }
#
#   # Construction of the matrix of selected pairs
#   switch(
#     typeEstCKT,
#
#     1 = {dataMatrix <- datasetPairs_hCV1(
#       X1 = observedX1, X2 = observedX2, Z = observedZ,
#       h, cut = 0.5, onlyConsecutivePairs = FALSE, nPairs = NULL)},
#     datasetPairs_hCV1(X1, X2, Z, h, cut = 0.5, onlyConsecutivePairs = FALSE, nPairs = NULL)
#   )
#   dataMatrix = computeDataMatrix_hCV(X1, X2, Z = vectorZ, cut = cut_a, typeEstCKT = typeEstCKT)
#   nPairsReal = nrow(dataMatrix)
#   if(nPairsReal != nPairs) {
#     warning(paste0("nPairs =", nPairs, " but real number of pairs is ", nPairsReal))
#   }
#
#   matrixScore = matrix(nrow = nPairsReal, ncol = n_h)
#   for (i_h in 1:n_h)
#   {
#     h = liste_h[i_h]
#     pbapply::pbapply(X = array(1:nPairsReal), FUN = function(i)
#     {
#       matrixScore[i, i_h] <<- computeEstimateNP_tau_pointwise(
#         matrixSignsPairs = matrixSignsPairs[-c(dataMatrix[i,4], dataMatrix[i,5]) ,
#                                             -c(dataMatrix[i,4], dataMatrix[i,5])],
#         vectorZ = vectorZ[-c(dataMatrix[i,4], dataMatrix[i,5])],
#         h = h, pointZ = dataMatrix[i,2], kernel.name = kernel.name,
#         typeEstCKT = typeEstCKT)
#     }, MARGIN = 1 )
#     scores[i_h] = mean( (dataMatrix[,1] - matrixScore[, i_h])^2, na.rm = TRUE )
#     if (verbose){print(paste0("h = ", h, " ; score = ", scores[i_h]))}
#   }
#   hCV = liste_h[which.min(scores)]
#   return(list(hCV = hCV, scores = scores))
# }
