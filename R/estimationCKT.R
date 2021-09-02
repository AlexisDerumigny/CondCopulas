


#' Estimation of conditional Kendall's tau
#' between two variables X1 and X2 given Z = z
#'
#'
#' @param observedX1 a vector of n observations of the first variable
#' @param observedX2 a vector of n observations of the second variable
#' @param observedZ a vector of n observations of the conditioning variable,
#' or a matrix with n rows of observations of the conditioning vector
#' @param newData the new data of observations of Z at which
#' the conditional Kendall's tau should be estimated.
#' @param methodEstimation method for estimating the conditional Kendall's tau
#' @param h the bandwidth
#' @param listPhi the list of transformations to be applied
#' (in case of regression-like models)
#'
#' @param ... other parameters passed to the estimating functions
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.fit.nNets}},
#' \code{\link{CKT.predict.kNN}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.kernel}} and \code{\link{CKT.kendallReg.fit}}.
#'
#' @return the vector of estimated conditional Kendall's tau
#' at each of the observations of \code{newData}.
#'
#' @seealso the specialized functions for estimating
#' conditional Kendall's tau:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.fit.nNets}},
#' \code{\link{CKT.predict.kNN}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.kernel}} and \code{\link{CKT.kendallReg.fit}}.
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#'
#' Derumigny, A., & Fermanian, J. D. (2019).
#' On kernel-based estimation of conditional Kendall’s tau:
#' finite-distance bounds and asymptotic behavior.
#' Dependence Modeling, 7(1), 292-321.
#'
#' Derumigny, A., & Fermanian, J. D. (2020).
#' On Kendall’s regression.
#' Journal of Multivariate Analysis, 178, 104610.
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
#' estimatedCKT_tree <- CKT.estimate(
#'   observedX1 = X1, observedX2 = X2, observedZ = Z,
#'   newData = newZ,
#'   methodEstimation = "tree", h = 0.07)
#'
#' estimatedCKT_rf <- CKT.estimate(
#'   observedX1 = X1, observedX2 = X2, observedZ = Z,
#'   newData = newZ,
#'   methodEstimation = "randomForest", h = 0.07)
#'
#' estimatedCKT_GLM <- CKT.estimate(
#'   observedX1 = X1, observedX2 = X2, observedZ = Z,
#'   newData = newZ,
#'   methodEstimation = "logit", h = 0.07,
#'   listPhi = list(function(x){return(x)}, function(x){return(x^2)},
#'                  function(x){return(x^3)}) )
#'
#' estimatedCKT_kNN <- CKT.estimate(
#'   observedX1 = X1, observedX2 = X2, observedZ = Z,
#'   newData = newZ,
#'   methodEstimation = "nearestNeighbors", h = 0.07,
#'   number_nn = c(50,80, 100, 120,200),
#'   partition = 4
#'   )
#'
#' estimatedCKT_nNet <- CKT.estimate(
#'   observedX1 = X1, observedX2 = X2, observedZ = Z,
#'   newData = newZ,
#'   methodEstimation = "neuralNetwork", h = 0.07,
#'   )
#'
#' estimatedCKT_kernel <- CKT.estimate(
#'   observedX1 = X1, observedX2 = X2, observedZ = Z,
#'   newData = newZ,
#'   methodEstimation = "kernel", h = 0.07,
#'   )
#'
#' estimatedCKT_kendallReg <- CKT.estimate(
#'    observedX1 = X1, observedX2 = X2, observedZ = Z,
#'    newData = newZ,
#'    methodEstimation = "kendallReg", h = 0.07)
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in other colors)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col="black",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_tree, col = "red")
#' lines(newZ, estimatedCKT_rf, col = "blue")
#' lines(newZ, estimatedCKT_GLM, col = "green")
#' lines(newZ, estimatedCKT_kNN, col = "purple")
#' lines(newZ, estimatedCKT_nNet, col = "coral")
#' lines(newZ, estimatedCKT_kernel, col = "skyblue")
#' lines(newZ, estimatedCKT_kendallReg, col = "darkgreen")
#'
#' @export
#'
CKT.estimate <- function (
  observedX1, observedX2, observedZ,
  newData = observedZ, methodEstimation, h,
  listPhi = if(methodEstimation == "kendallReg")
            {list(function(x){return(x)}, function(x){return(x^2)},
                  function(x){return(x^3)})} else
              {list(identity)} , ...)
{
  if (methodEstimation %in% c("tree", "randomForest", "logit", "probit",
                              "nearestNeighbors", "neuralNetwork")){
    datasetPairs = datasetPairs(X1 = observedX1, X2 = observedX2,
                                Z = observedZ, h = h)
  }

  if (methodEstimation %in% c("logit", "probit",
                              "nearestNeighbors", "neuralNetwork")){
    designMatrix =
      sapply(listPhi,
             function(x) sapply(datasetPairs[,2:(ncol(datasetPairs)-3),drop=FALSE],
                                x))
    newDesignMatrix =
      sapply(listPhi,
             function(x) sapply(newData, x))
  }

  switch (
    methodEstimation,

    "tree" = {
      fit = CKT.fit.tree(datasetPairs = datasetPairs, ...)
      estCKT = CKT.predict.tree(fit = fit, newData = newData)
    },

    "randomForest" = {
      fit = CKT.fit.randomForest(datasetPairs = datasetPairs, n = length(observedX1), ...)
      estCKT = CKT.predict.randomForest(fit = fit, newData = newData)
    },

    "logit" = {
      fit = CKT.fit.GLM(datasetPairs = datasetPairs,
                        designMatrix = designMatrix, link = "logit", ...)

      estCKT = CKT.predict.GLM(fit = fit, newData = newDesignMatrix)
    },

    "probit" = {
      fit = CKT.fit.GLM(datasetPairs = datasetPairs,
                        designMatrix = designMatrix, link = "probit", ...)

      estCKT = CKT.predict.GLM(fit = fit, newData = newDesignMatrix)
    },

    "nearestNeighbors" = {
      result = CKT.predict.kNN(datasetPairs = datasetPairs,
                               designMatrix = designMatrix,
                               newData = newDesignMatrix, ...)
      estCKT = result$estimatedCKT
    },

    "neuralNetwork" = {
      fit = CKT.fit.nNets(datasetPairs = datasetPairs,
                          designMatrix = designMatrix, ...)

      estCKT = CKT.predict.nNets(fit = fit, newData = newDesignMatrix)
    },

    "kernel" = {
      result = CKT.kernel(observedX1 = observedX1, observedX2 = observedX2,
                          observedZ = observedZ, newZ = newData,
                          h = h, ...)
      estCKT = result$estimatedCKT
    },

    "kendallReg" = {

      ZToEstimate <- newData

      designMatrixZ =
          sapply(listPhi,
                 function(x) sapply(ZToEstimate, x))

      result = CKT.kendallReg.fit(observedX1 = observedX1, observedX2 = observedX2,
                                  observedZ = observedZ, newData = designMatrixZ,
                                  ZToEstimate = ZToEstimate,
                                  designMatrixZ = designMatrixZ,
                                  h_kernel = h, h_lambda = h, ...)
      estCKT = result$estimatedCKT
    },

    {
      stop("'methodEstimation' not found. Possible choices are:",
           "'tree', 'randomForest', 'logit', 'probit',
           'nearestNeighbors', 'neuralNetwork', 'kernel', 'kendallReg'.")
    }
  )

  return (estCKT)
}

