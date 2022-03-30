
#' Estimation of conditional Kendall's taus using a classification tree
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
#' These functions estimate and predict conditional Kendall's tau using a
#' \strong{classification tree}. This is possible by the relationship between
#' estimation of conditional Kendall's tau and classification problems
#' (see Derumigny and Fermanian (2019)): estimation of conditional Kendall's tau
#' is equivalent to the prediction of concordance in the space of pairs
#' of observations.
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#'
#' @param mindev a factor giving the minimum deviation for a node to be splitted.
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#' @param mincut the minimum number of observations (of pairs) in a node
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#'
#' @return \code{CKT.fit.tree} returns the fitted tree.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Section 3.2)
#' \doi{10.1016/j.csda.2019.01.013}
#'
#' @seealso See also other estimators of conditional Kendall's tau:
#' \code{\link{CKT.fit.nNets}}, \code{\link{CKT.fit.randomForest}},
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
#' datasetP = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
#' est_Tree = CKT.fit.tree(datasetPairs = datasetP, mindev = 0.008)
#' print(est_Tree)
#'
#' newZ = seq(1,10,by = 0.1)
#' prediction = CKT.predict.tree(fit = est_Tree, newZ = data.frame(x=newZ))
#' # Comparison between true Kendall's tau (in red)
#' # and estimated Kendall's tau (in black)
#' plot(newZ, prediction, type = "l", ylim = c(-1,1))
#' lines(newZ, -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2), col="red")
#'
#' @export
#'
CKT.fit.tree <- function(datasetPairs, mindev = 0.008, mincut = 0)
{
  dim_Z = ncol(datasetPairs) - 4
  designMatrix_data = data.frame(x = datasetPairs[ , 2:(1+dim_Z)])

  fit_tree <- tree::tree(
    factor(datasetPairs[ , 1], levels = c(0,1)) ~ ., data = designMatrix_data,
    weights = datasetPairs[ , 2 + dim_Z], mindev = mindev, mincut = mincut)

  return(fit_tree)
}


#' Predict the values of conditional Kendall's tau using a classification tree
#'
#' @param fit result of a call to \code{CKT.fit.tree}
#'
#' @param newZ new matrix of observations, with the same number of variables.
#' and same names as the \code{designMatrix} that was used to fit the tree.
#'
#' @return \code{CKT.predict.tree} returns
#' a vector of (predicted) conditional Kendall's taus of the same size
#' as the number of rows of \code{newZ}.
#'
#'
#' @importFrom tree tree
#'
#' @rdname CKT.fit.tree
#' @export
#'
CKT.predict.tree <- function(fit, newZ){

  # if(length(newZ[1,])==1){colnames(newZ) <- "x"}
  pred = stats::predict(fit, newdata = data.frame(x=newZ))
  prediction = 2 * pred[,2] - 1

  return(prediction)
}


#' Fit a Random Forest that can be used for the estimation of conditional Kendall's tau.
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
#' These functions estimate and predict conditional Kendall's tau using a
#' \strong{random forest}. This is possible by the relationship between
#' estimation of conditional Kendall's tau and classification problems
#' (see Derumigny and Fermanian (2019)): estimation of conditional Kendall's tau
#' is equivalent to the prediction of concordance in the space of pairs
#' of observations.
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#' @param designMatrix the matrix of predictor to be used for the fitting of the tree
#' @param n the original sample size of the dataset
#'
#' @param nTree number of trees of the Random Forest.
#' @param nObs_per_Tree number of observations kept in each tree.
#' @param nVar_per_Tree number of variables kept in each tree.
#'
#' @param mindev a factor giving the minimum deviation for a node to be splitted.
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#'
#' @param mincut the minimum number of observations (of pairs) in a node
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#'
#' @param verbose if \code{TRUE}, a message is printed after fitting each tree.
#'
#' @param nMaxDepthAllowed the maximum number of errors of type
#' "the tree cannot be fitted" or "is too deep" before stopping the procedure.
#'
#' @return a list with two components
#' \itemize{
#'     \item \code{list_tree} a list of size \code{nTree}
#'     composed of all the fitted trees.
#'
#'     \item \code{list_variables} a list of size \code{nTree}
#'     composed of the (predictor) variables for each tree.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 4)
#' \doi{10.1016/j.csda.2019.01.013}
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
#' datasetP = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
#' est_RF = CKT.fit.randomForest(datasetPairs = datasetP, n = N,
#'   mindev = 0.008)
#'
#' newZ = seq(1,10,by = 0.1)
#' prediction = CKT.predict.randomForest(fit = est_RF,
#'    newZ = data.frame(x=newZ))
#' # Comparison between true Kendall's tau (in red)
#' # and estimated Kendall's tau (in black)
#' plot(newZ, prediction, type = "l", ylim = c(-1,1))
#' lines(newZ, -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2), col="red")
#'
#' @export
#'
CKT.fit.randomForest <- function(
  datasetPairs, designMatrix = data.frame(x = datasetPairs[ , 2:(ncol(datasetPairs) - 3)]),
  n, nTree = 10, mindev = 0.008, mincut = 0,
  nObs_per_Tree = ceiling(0.8*n),
  nVar_per_Tree = ceiling(0.8*(ncol(datasetPairs)-4)),
  verbose = FALSE, nMaxDepthAllowed = 10)
{
  mindevTree = rep(mindev, length.out = nTree)
  mincutTree = rep(mincut, length.out = nTree)
  dimZ = ncol(datasetPairs) - 4
  if (n < nObs_per_Tree){
    stop("Number of observations per tree greater ",
         "than the number of observations in the sample.")
  }

  list_tree = as.list(seq(nTree))
  list_variables = as.list(seq(nTree))
  iTree = 1
  nErrors = 0
  while (iTree <= nTree){
    tryCatch({
      # We choose at random the initial sample
      sampleId = sample.int(n = n, size = nObs_per_Tree)

      # We compute the data matrix from this sample
      whichPairs = which( (datasetPairs[, "iUsed"] %in% sampleId)
                          & (datasetPairs[, "jUsed"] %in% sampleId) )
      datasetPairs_sample = datasetPairs[whichPairs,]

      # We compute the inputMatrix
      inputMatrix_sample = data.frame(x = designMatrix[whichPairs,])

      # We select at random the variables used
      sampleVariable = sample.int(n = length(inputMatrix_sample[1,]), size = nVar_per_Tree)

      # We fit the tree using this variable
      list_tree[[iTree]] <- tree::tree(
        factor(datasetPairs_sample[ , 1], levels = c(0,1)) ~ . ,
        data = data.frame(x = inputMatrix_sample[,sampleVariable]) ,
        weights = datasetPairs_sample[ , 2 + dimZ] ,
        mindev = mindevTree[iTree],
        mincut = mincutTree[iTree])

      list_variables[[iTree]] <- sampleVariable
    },
    error = function(e){
      if(e$message %in% c("maximum depth reached\n","tree is too big")) {
        iTree <<- iTree - 1
        # warning("maximum depth reached for one tree. Trying again with another sample...")
        nErrors = nErrors + 1
        if (nErrors > nMaxDepthAllowed){
          stop("Maximum depth reach too often. Maybe mindev or nMaxDepthAllowed are too small.")
        }
      } else {
        stop(e)
      }
      e
    })
    if (verbose){ cat(iTree); cat(" ") }
    iTree = iTree + 1
  }
  if (nErrors > 0){
    warning("maximum depth reached for ", nErrors,
            " tree(s).")
  }
  return( list(list_tree = list_tree, list_variables = list_variables) )
}


#' Predict the values of conditional Kendall's tau using Random Forests
#'
#' @param fit result of a call to \code{CKT.fit.randomForest}.
#'
#' @param newZ new matrix of observations, with the same number of variables.
#' and same names as the \code{designMatrix} that was used to fit the Random Forest.
#'
#' @return \code{CKT.predict.randomForest} returns
#' a vector of (predicted) conditional Kendall's taus of the same size
#' as the number of rows of the newZ.
#'
#' @rdname CKT.fit.randomForest
#' @export
#'
CKT.predict.randomForest <- function(fit, newZ)
{
  prediction = 0
  for (iTree in 1:length(fit$list_tree)){
    prediction = prediction +
      stats::predict(
        fit$list_tree[[iTree]],
        newdata = data.frame(x = newZ)[,fit$list_variables[[iTree]]
                                          ,drop=FALSE]) [, 2]
  }
  prediction = prediction / length(fit$list_tree)
  return (2*prediction-1)
}



