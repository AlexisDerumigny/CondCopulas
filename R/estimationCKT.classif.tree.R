
#' Estimation of conditional Kendall's taus using a classification tree
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#'
#' @param mindev a factor giving the minimum deviation for a node to be splitted.
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#' @param mincut the minimum number of observations (of pairs) in a node
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#'
#' @return the fitted tree.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Section 3.2)
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


#' Predict the values of conditional Kendall's tau a classification tree
#'
#' @param fit result of a call to CKT.fit.tree
#' @param newData new matrix of observations, with the same number of variables.
#' and same names as the designMatrix that was used to fit the tree.
#'
#' @return a vector of (predicted) conditional Kendall's taus of the same size
#' as the number of rows of the newData.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Section 3.2)
#'
#' @importFrom tree tree
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
#'
#' newData = seq(1,10,by = 0.1)
#' prediction = CKT.predict.tree(fit = est_Tree, newData = data.frame(x=newData))
#' # Comparison between true Kendall's tau (in red)
#' # and estimated Kendall's tau (in black)
#' plot(newData, prediction, type = "l", ylim = c(-1,1))
#' lines(newData, -0.9 + 1.8 * pnorm(newData, mean = 5, sd = 2), col="red")
#'
#' @export
#'
CKT.predict.tree <- function(fit, newData){

  # if(length(newData[1,])==1){colnames(newData) <- "x"}
  pred = stats::predict(fit, newdata = data.frame(x=newData))
  prediction = 2 * pred[,2] - 1

  return(prediction)
}


#' Fit a Random Forest that can be used for the estimation of conditional Kendall's tau.
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#' @param designMatrix the matrix of predictor to be used for the fitting of the tree
#' @param n the original sample size of your dataset
#'
#' @param nTree number of trees of the Random Forest.
#' @param nObs_per_Tree number of observations kept in each tree.
#' @param nVar_per_Tree number of variables kept in each tree.
#'
#' @param mindev a factor giving the minimum deviation for a node to be splitted.
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#' @param mincut the minimum number of observations (of pairs) in a node
#' See \code{tree::\link[tree]{tree.control}()} for more details.
#'
#' @param verbose if \code{TRUE}, a message is printed after fitting each tree.
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
#' @param fit result of a call to CKT.fit.randomForest.
#' @param newData new matrix of observations, with the same number of variables.
#' and same names as the designMatrix that was used to fit the Random Forest.
#'
#' @return a vector of (predicted) conditional Kendall's taus of the same size
#' as the number of rows of the newData.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 4)
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
#'   mindev = 0.008, nVar_per_Tree = 1, nTree = 10, nObs_per_Tree = 600)
#'
#' newData = seq(1,10,by = 0.1)
#' prediction = CKT.predict.randomForest(fit = est_RF,
#'    newData = data.frame(x=newData))
#' # Comparison between true Kendall's tau (in red)
#' # and estimated Kendall's tau (in black)
#' plot(newData, prediction, type = "l", ylim = c(-1,1))
#' lines(newData, -0.9 + 1.8 * pnorm(newData, mean = 5, sd = 2), col="red")
#'
#' @export
#'
CKT.predict.randomForest <- function(fit, newData)
{
  prediction = 0
  for (iTree in 1:length(fit$list_tree)){
    prediction = prediction +
      stats::predict(
        fit$list_tree[[iTree]],
        newdata = data.frame(x = newData)[,fit$list_variables[[iTree]]
                                          ,drop=FALSE]) [, 2]
  }
  prediction = prediction / length(fit$list_tree)
  return (2*prediction-1)
}



