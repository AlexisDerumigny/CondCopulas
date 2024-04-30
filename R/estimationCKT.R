
#' Estimation of conditional Kendall's tau
#' between two variables X1 and X2 given Z = z
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
#' This function can use different estimators for conditional Kendall's tau,
#' see the description of the parameter \code{methodEstimation}
#' for a complete list of possibilities.
#'
#'
#' @param X1 a vector of \eqn{n} observations of the first variable
#'
#' @param X2 a vector of \eqn{n} observations of the second variable
#'
#' @param Z a vector of \eqn{n} observations of the conditioning variable,
#' or a matrix with \eqn{n} rows of observations of the conditioning vector
#' (if \eqn{Z} is multivariate).
#'
#' @param newZ the new values for the conditioning variable \eqn{Z}
#' at which the conditional Kendall's tau should be estimated.
#' \itemize{
#'    \item If \code{observedZ} is a vector,
#'    then \code{newZ} must be a vector as well.
#'    \item If \code{observedZ} is a matrix,
#'    then \code{newZ} must be a matrix as well, with the same number of columns
#'    ( = the dimension of \eqn{Z}).
#' }
#'
#' @param methodEstimation method for estimating the conditional Kendall's tau.
#' Possible estimation methods are:
#' \itemize{
#'    \item \code{"kernel"}: kernel smoothing,
#'    as described in (Derumigny, & Fermanian (2019a))
#'
#'    \item \code{"kendallReg"}: regression-type model,
#'    as described in (Derumigny, & Fermanian (2020))
#'
#'    \item \code{"tree"}, \code{"randomForest"},
#'    \code{"logit"}, and \code{"neuralNetwork"}:
#'    use the relationship between conditional Kendall's tau
#'    and classification problems to use the respective classification algorithms
#'    for the estimation of conditional Kendall's tau,
#'    as described in (Derumigny, & Fermanian (2019b))
#' }
#'
#' @param h the bandwidth
#'
#' @param listPhi the list of transformations to be applied
#' to the conditioning variable \eqn{Z}
#' (in case of regression-type models).
#'
#' @param ... other parameters passed to the estimating functions
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.fit.nNets}},
#' \code{\link{CKT.predict.kNN}}, \code{\link{CKT.kernel}}
#' and \code{\link{CKT.kendallReg.fit}}.
#'
#' @param observedX1,observedX2,observedZ old parameter names for \code{X1},
#' \code{X2}, \code{Z}. Support for this will be removed at a later version.
#'
#'
#' @return the vector of estimated conditional Kendall's tau
#' at each of the observations of \code{newZ}.
#'
#' @seealso the specialized functions for estimating
#' conditional Kendall's tau for each method:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.fit.nNets}},
#' \code{\link{CKT.predict.kNN}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.kernel}} and \code{\link{CKT.kendallReg.fit}}.
#'
#' See also the nonparametric estimator of conditional copula models
#' \code{\link{estimateNPCondCopula}},
#' and the parametric estimators of conditional copula models
#' \code{\link{estimateParCondCopula}}.
#'
#' In the case where \eqn{Z} is discrete
#' or in the case of discrete conditioning events, see
#' \code{\link{bCond.treeCKT}}.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019a).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' \doi{10.1016/j.csda.2019.01.013}
#'
#' Derumigny, A., & Fermanian, J. D. (2019b).
#' On kernel-based estimation of conditional Kendall’s tau:
#' finite-distance bounds and asymptotic behavior.
#' Dependence Modeling, 7(1), 292-321.
#' \doi{10.1515/demo-2019-0016}
#'
#' Derumigny, A., & Fermanian, J. D. (2020).
#' On Kendall’s regression.
#' Journal of Multivariate Analysis, 178, 104610.
#' \doi{10.1016/j.jmva.2020.104610}
#'
#' @usage
#' CKT.estimate(
#'   X1 = NULL, X2 = NULL, Z = NULL,
#'   newZ = Z, methodEstimation, h,
#'   listPhi = if(methodEstimation == "kendallReg")
#'                {list( function(x){return(x)}   ,
#'                       function(x){return(x^2)} ,
#'                       function(x){return(x^3)} )
#'                } else {list(identity)} ,
#'   ... ,
#'   observedX1 = NULL, observedX2 = NULL, observedZ = NULL )
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 300
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' newZ = seq(2,10,by = 0.1)
#' h = 0.1
#' estimatedCKT_tree <- CKT.estimate(
#'   X1 = X1, X2 = X2, Z = Z,
#'   newZ = newZ,
#'   methodEstimation = "tree", h = h)
#'
#' estimatedCKT_rf <- CKT.estimate(
#'   X1 = X1, X2 = X2, Z = Z,
#'   newZ = newZ,
#'   methodEstimation = "randomForest", h = h)
#'
#' estimatedCKT_GLM <- CKT.estimate(
#'   X1 = X1, X2 = X2, Z = Z,
#'   newZ = newZ,
#'   methodEstimation = "logit", h = h,
#'   listPhi = list(function(x){return(x)}, function(x){return(x^2)},
#'                  function(x){return(x^3)}) )
#'
#' estimatedCKT_kNN <- CKT.estimate(
#'   X1 = X1, X2 = X2, Z = Z,
#'   newZ = newZ,
#'   methodEstimation = "nearestNeighbors", h = h,
#'   number_nn = c(50,80, 100, 120,200),
#'   partition = 4
#'   )
#'
#' estimatedCKT_nNet <- CKT.estimate(
#'   X1 = X1, X2 = X2, Z = Z,
#'   newZ = newZ,
#'   methodEstimation = "neuralNetwork", h = h,
#'   )
#'
#' estimatedCKT_kernel <- CKT.estimate(
#'   X1 = X1, X2 = X2, Z = Z,
#'   newZ = newZ,
#'   methodEstimation = "kernel", h = h,
#'   )
#'
#' estimatedCKT_kendallReg <- CKT.estimate(
#'    X1 = X1, X2 = X2, Z = Z,
#'    newZ = newZ,
#'    methodEstimation = "kendallReg", h = h)
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
  X1 = NULL, X2 = NULL, Z = NULL,
  newZ = Z, methodEstimation, h,
  listPhi = if(methodEstimation == "kendallReg")
            {list(function(x){return(x)}, function(x){return(x^2)},
                  function(x){return(x^3)})} else
              {list(identity)} , ... ,
  observedX1 = NULL, observedX2 = NULL, observedZ = NULL
  )
{
  if (length(newZ) == 0){
    warning("'newZ' is of length 0, therefore, no estimation is done.")
    return (numeric(0))
  }

  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  env = environment()
  .observedX1X2_to_X1X2(env)
  .observedZ_to_Z(env)

  if (methodEstimation %in% c("tree", "randomForest", "logit", "probit",
                              "nearestNeighbors", "neuralNetwork")){
    datasetPairs = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = h)
  }

  if (methodEstimation %in% c("logit", "probit",
                              "nearestNeighbors", "neuralNetwork")){
    designMatrix =
      sapply(listPhi,
             function(x) sapply(datasetPairs[,2:(ncol(datasetPairs)-3),drop=FALSE],
                                x))
    newDesignMatrix =
      sapply(listPhi,
             function(x) sapply(newZ, x))
  }

  switch (
    methodEstimation,

    "tree" = {
      fit = CKT.fit.tree(datasetPairs = datasetPairs, ...)
      estCKT = CKT.predict.tree(fit = fit, newZ = newZ)
    },

    "randomForest" = {
      fit = CKT.fit.randomForest(datasetPairs = datasetPairs, n = length(X1), ...)
      estCKT = CKT.predict.randomForest(fit = fit, newZ = newZ)
    },

    "logit" = {
      fit = CKT.fit.GLM(datasetPairs = datasetPairs,
                        designMatrix = designMatrix, link = "logit", ...)

      estCKT = CKT.predict.GLM(fit = fit, newZ = newDesignMatrix)
    },

    "probit" = {
      fit = CKT.fit.GLM(datasetPairs = datasetPairs,
                        designMatrix = designMatrix, link = "probit", ...)

      estCKT = CKT.predict.GLM(fit = fit, newZ = newDesignMatrix)
    },

    "nearestNeighbors" = {
      result = CKT.predict.kNN(datasetPairs = datasetPairs,
                               designMatrix = designMatrix,
                               newZ = newDesignMatrix, ...)
      estCKT = result$estimatedCKT
    },

    "neuralNetwork" = {
      fit = CKT.fit.nNets(datasetPairs = datasetPairs,
                          designMatrix = designMatrix, ...)

      estCKT = CKT.predict.nNets(fit = fit, newZ = newDesignMatrix)
    },

    "kernel" = {
      result = CKT.kernel(X1 = X1, X2 = X2, Z = Z, newZ = newZ, h = h, ...)
      estCKT = result$estimatedCKT
    },

    "kendallReg" = {

      ZToEstimate <- newZ

      designMatrixZ =
          sapply(listPhi,
                 function(x) sapply(ZToEstimate, x))

      result = CKT.kendallReg.fit(X1 = X1, X2 = X2, Z = Z,
                                  newZ = designMatrixZ,
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

