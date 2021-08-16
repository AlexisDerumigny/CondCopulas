


# Cross-validation procedures --------------------------------------------------------


#' Leave-one-out cross validation for choosing the bandwidth
#' for kernel-based range_h of conditional Kendall's tau
#'
#' @param observedX1 a vector of n observations of the first variable
#' @param observedX2 a vector of n observations of the second variable
#' @param observedZ a vector of n observations of the conditioning variable,
#' or a matrix with n rows of observations of the conditioning vector
#' @param range_h vector containing possible values for the bandwidth.
#'
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#' @param nPairs number of pairs used in the cross-validation criteria.
#'
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#'
#' @param progressBar if TRUE, a progressbar for each h is displayed
#' to show the progress of the computation.
#' @param verbose if TRUE, print the score of each h during the procedure.
#'
#'
#' @return A list with components:
#' \itemize{
#'     \item \code{hCV}: the chosen bandwidth
#'     \item \code{scores}: vector of the same length as range_h giving the
#'     value of the CV criteria for each of the h tested
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' On kernel-based estimation of conditional Kendall’s tau:
#' finite-distance bounds and asymptotic behavior.
#' Dependence Modeling, 7(1), 292-321.
#' Page 296, Equation (4).
#'
#' @export
#'
CKT.hCV.l1out <- function (observedX1, observedX2, observedZ,
                           range_h, matrixSignsPairs,
                           nPairs = 10*length(observedX1),
                           typeEstCKT = 4, kernel.name = "Epa",
                           progressBar = TRUE, verbose = TRUE)
{
  n_h = length(range_h)
  scores = rep(NA, n_h)

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
      matrixScore[i, i_h] <<- CKT.kernelPointwise.univariate(
        matrixSignsPairs = matrixSignsPairs[-c(dataMatrix[i,4], dataMatrix[i,5]) ,
                                            -c(dataMatrix[i,4], dataMatrix[i,5])],
        vectorZ = observedZ[-c(dataMatrix[i,4], dataMatrix[i,5])],
        h = h, pointZ = dataMatrix[i,2],
        kernel.name = kernel.name, typeEstCKT = typeEstCKT)
    }
  } else {
    dimZ = ncol(observedZ)
    computeScore <- function(i, i_h, h){
      matrixScore[i, i_h] <<- CKT.kernelPointwise.multivariate(
        matrixSignsPairs = matrixSignsPairs[-c(dataMatrix[i,4], dataMatrix[i,5]) ,
                                            -c(dataMatrix[i,4], dataMatrix[i,5])],
        matrixZ = observedZ[-c(dataMatrix[i,4], dataMatrix[i,5])],
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
#' using n-fold cross-validation
#'
#'
#' @param range_h vector containing possible values for the bandwidth.
#'
#' @param matrixSignsPairs square matrix of signs of all pairs,
#' produced by computeMatrixSignPairs.
#'
#' @param observedZ vector of observed values of Z.
#' If Z is multivariate, then this is a matrix whose rows correspond
#' to the observations of Z
#'
#' @param ZToEstimate vector of fixed conditioning values at which
#' the difference between the two conditional Kendall's tau should be computed.
#' Can also be a matrix whose lines are the conditioning vectors at which
#' the difference between the two conditional Kendall's tau should be computed.
#'
#' @param Kfolds number of subsamples used.
#'
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#'
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#' @param progressBar if TRUE, a progressbar for each h is displayed
#' to show the progress of the computation.
#'
#'
#' @return A list with components:
#' \itemize{
#'     \item \code{hCV}: the chosen bandwidth
#'     \item \code{scores}: vector of the same length as range_h giving the
#'     value of the CV criteria for each of the h tested
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' On kernel-based estimation of conditional Kendall’s tau:
#' finite-distance bounds and asymptotic behavior.
#' Dependence Modeling, 7(1), 292-321.
#'
#' @export
#'
CKT.hCV.Kfolds <- function(range_h, matrixSignsPairs, observedZ, ZToEstimate,
                           typeEstCKT = 4, kernel.name = "Epa",
                           Kfolds = 5, progressBar = TRUE)
{
  n_h = length(range_h)
  scores = rep(NA, n_h)
  n = NROW(observedZ)

  if (is.vector(observedZ)){
    estimCKTNP <- CKT.kernel.univariate
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
        matrixSignsPairs = matrixSignsPairs[!which, !which],
        observedZ = observedZ[!which,], ZToEstimate = ZToEstimate,
        typeEstCKT = typeEstCKT, h = h, kernel.name = kernel.name, progressBar = FALSE)

      list_vectorEstimate_comp[[i]] = estimCKTNP(
        matrixSignsPairs = matrixSignsPairs[which, which],
        observedZ = observedZ[which,], ZToEstimate = ZToEstimate,
        typeEstCKT = typeEstCKT, h = h, kernel.name = kernel.name, progressBar = FALSE)

      if (all(!is.finite(list_vectorEstimate[[i]])))
      {
        distanceTot = NA ; break
      }
      distance = sum((list_vectorEstimate[[i]] - list_vectorEstimate_comp[[i]])^2)
      print(distance)
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
