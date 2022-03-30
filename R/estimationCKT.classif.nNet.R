

#' Estimation of conditional Kendall's taus by model averaging of neural networks
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
#' This function estimates conditional Kendall's tau using
#' \strong{neural networks}. This is possible by the relationship between
#' estimation of conditional Kendall's tau and classification problems
#' (see Derumigny and Fermanian (2019)): estimation of conditional Kendall's tau
#' is equivalent to the prediction of concordance in the space of pairs
#' of observations.
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#' @param designMatrix the matrix of predictor to be used for the fitting of the tree
#'
#' @param vecSize vector with the number of neurons for each network
#'
#' @param nObs_per_NN number of observations used for each neural network.
#'
#' @param verbose a number indicated what to print
#' \itemize{
#'     \item \code{0}: nothing printed at all.
#'     \item \code{1}: a message is printed at the convergence of each neural network.
#'     \item \code{2}: details are printed for each optimization of each network.
#' }
#'
#' @return \code{CKT.fit.nNets} returns a list of the fitted neural networks
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 7)
#' \doi{10.1016/j.csda.2019.01.013}
#'
#' @seealso See also other estimators of conditional Kendall's tau:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.GLM}}, \code{\link{CKT.predict.kNN}},
#' \code{\link{CKT.kernel}}, \code{\link{CKT.kendallReg.fit}},
#' and the more general wrapper \code{\link{CKT.estimate}}.
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
#'
#' fitCKT_nets <- CKT.fit.nNets(datasetPairs = datasetP)
#' estimatedCKT_nNets <- CKT.predict.nNets(
#'   fit = fitCKT_nets, newZ = matrix(newZ, ncol = 1))
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in red)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col="black",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_nNets, col = "red")
#'
#'
#' @export
#'
CKT.fit.nNets <- function(datasetPairs,
                          designMatrix = datasetPairs[,2:(ncol(datasetPairs)-3),drop=FALSE],
                          vecSize = rep(3, times = 10),
                          nObs_per_NN = 0.9*nrow(designMatrix),
                          verbose = 1)
{
  # Initialization
  length_list_nnet = length(vecSize)
  n = nrow(designMatrix)

  list_nnet = as.list(seq(length_list_nnet))

  for (i in 1:length_list_nnet)
  {
    # We choose at random the initial sample
    sampleId = sample.int(n = n, size = nObs_per_NN)

    # We compute the data matrix from this sample
    whichPairs = which( (datasetPairs[,"iUsed"] %in% sampleId)
                        & (datasetPairs[,"jUsed"] %in% sampleId) )
    datasetPairs_sample = datasetPairs[whichPairs,]

    designMatrix_sample = cbind(designMatrix[whichPairs,])

    fit_nnet = nnet::nnet(as.factor(datasetPairs_sample[ ,1]) ~ .,
                          data = designMatrix_sample, trace = (verbose == 2),
                          weights = datasetPairs_sample[, "kernel.value"],
                          size = vecSize[i], maxit = 1000)

    list_nnet[[i]] <- fit_nnet

    if (verbose >= 1)
    {
      cat(i) ; cat(" -- size = ") ; cat(vecSize[i]) ; cat("\n")
    }
  }

  return(list_nnet)
}


#' Predict the values of conditional Kendall's tau
#' using Model Averaging of Neural Networks
#'
#' @param fit result of a call to \code{CKT.fit.nNet}
#'
#' @param newZ new matrix of observations, with the same number of variables.
#' and same names as the \code{designMatrix} that was used to fit the neural networks.
#'
#' @param aggregationMethod the method to be used to aggregate all the predictions
#' together. Can be \code{"mean"} or \code{"median"}.
#'
#' @return \code{CKT.predict.nNets} returns
#' a vector of (predicted) conditional Kendall's taus of the same size
#' as the number of rows of the matrix \code{newZ}.
#'
#' @importFrom nnet nnet
#'
#' @rdname CKT.fit.nNet
#' @export
#'
CKT.predict.nNets <- function(fit, newZ, aggregationMethod = "mean")
{
  length_list_nn = length(fit)
  length_toPredict = length(newZ[,1])

  matrixPrediction = matrix(nrow = length_toPredict, ncol = length_list_nn)

  # Prediction
  if (ncol(newZ) > 1) {
    for (i_nn in 1:length_list_nn)
    {
      matrixPrediction[, i_nn] = stats::predict(fit[[i_nn]], newZ)
    }
  } else {
    colnames(newZ) <- "V1"
    for (i_nn in 1:length_list_nn)
    {
      matrixPrediction[, i_nn] = stats::predict(fit[[i_nn]], newZ)
    }
  }

  # Aggregation
  prediction = rep(NA, length_toPredict)
  if (aggregationMethod == "mean")
  {
    prediction = rep(0, length_toPredict)
    for (i_nn in 1:length_list_nn)
    {
      prediction = prediction + matrixPrediction[, i_nn]
    }
    prediction = prediction / length_list_nn
  } else if (aggregationMethod == "median") {
    for (iPredict in 1:length_toPredict)
    {
      prediction[iPredict] = stats::median(matrixPrediction[iPredict, ])
    }
  }

  predictCKt = 2 * prediction - 1
  return(predictCKt)
}






