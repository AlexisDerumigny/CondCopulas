
#' Converting to matrix of indicators / matrix of conditional Kendall's tau
#'
#' The function \code{treeCKT2matrixInd}
#' takes as input a binary tree that has been returned
#' by the function \code{\link{bCond.treeCKT}}.
#' Since this tree describes a partition of the conditioning space,
#' it can be interesting to get, for a given dataset, the matrix
#' \deqn{1\{ X_{i,J} \in A_{j,J} \},}
#' where each \eqn{A_{j,J}} corresponds to a conditioning subset.
#' This is the so-called \code{matrixInd}.
#' Finally, it can be interesting to get the matrix of
#'
#' @param estimatedTree the tree that has been estimated before,
#' for example by \code{\link{bCond.treeCKT}}.
#'
#' @param newDataXI this is a matrix of size \code{N * |I|}
#' where \code{|I|} is the number of conditioned variables.
#' By default this is \code{NULL} meaning that
#' we return the matrix for the original data
#' used to compute the \code{estimatedTree}
#'
#' @param newDataXJ this is a matrix of size \code{N * |J|}
#' where \code{|J|} is the number of conditional variables used in the tree.
#' By default this is \code{NULL} meaning that
#' we return the matrix for the original data
#' (that was used to compute the \code{estimatedTree}).
#'
#' @param matrixInd a matrix of indexes of size (n, N.boxes) describing
#' for each observation i to which box ( = event) it belongs.
#'
#' @return \itemize{
#'   \item The function \code{treeCKT2matrixInd} returns
#'   a matrix of size \code{N * m} which component \code{[i,j]}
#'   is \deqn{1\{ X_{i,J} \in A_{j,J} \}}.
#'
#'   \item The function \code{matrixInd2matrixCKT} and \code{treeCKT2matrixCKT} return
#'   a matrix of size \code{|I| * (|I|-1) * m} where each component corresponds
#'   to a conditional Kendall's tau between a pair of conditional variables
#'   conditionally to the conditioned variables in one of the boxes
#' }
#'
#' @seealso \code{\link{bCond.treeCKT}} for the construction of such a binary tree.
#'
#' @examples
#' set.seed(1)
#' n = 200
#' XJ = MASS::mvrnorm(n = n, mu = c(3,3), Sigma = rbind(c(1, 0.2), c(0.2, 1)))
#' XI = matrix(nrow = n, ncol = 2)
#' high_XJ1 = which(XJ[,1] > 4)
#' XI[high_XJ1, ]  = MASS::mvrnorm(n = length(high_XJ1), mu = c(10,10),
#'                                 Sigma = rbind(c(1, 0.8), c(0.8, 1)))
#' XI[-high_XJ1, ] = MASS::mvrnorm(n = n - length(high_XJ1), mu = c(8,8),
#'                                 Sigma = rbind(c(1, -0.2), c(-0.2, 1)))
#'
#' result = bCond.treeCKT(XI = XI, XJ = XJ, minSize = 10, verbose = 2)
#'
#' treeCKT2matrixInd(result)
#'
#' matrixInd2matrixCKT(treeCKT2matrixInd(result), newDataXI = XI)
#'
#' treeCKT2matrixCKT(result)
#'
#' @name conv_treeCKT
NULL


#' Construction of the matrix of indicators given an estimated CKT tree
#'
#' @rdname conv_treeCKT
#' @export
treeCKT2matrixInd <- function(estimatedTree, newDataXJ = NULL)
{
  # We get all the leaves of the tree, which corresponds to the final boxes
  leaves <- data.tree::Traverse(estimatedTree, filterFun = data.tree::isLeaf)
  # We count the number of leaves
  m <- length(leaves)

  if (is.null(newDataXJ)) {
    # In this case, we just use the data stored in the tree
    n <- length(estimatedTree$condObs)
    matrixInd = matrix(data = 0, nrow = n, ncol = m)

    for (iBox in 1:m)
    {
      matrixInd[data.tree::GetAttribute(leaves[[iBox]], "condObs"), iBox] = 1
    }
  } else {
    # In this case, we use the new data to construct the matrix

    nprime <- nrow(newDataXJ)
    matrixInd = matrix(data = 0, nrow = nprime, ncol = m)

    # We fill each matrix column by column
    # for each box, which correspond to each leaf of the tree
    for (j in 1:m){
      leaf = leaves[[j]]
      vecIsInBoxj = rep(TRUE, nprime)
      # We go up in the parent chain to find all conditioning
      # i.e. all the edges of the box
      while( !leaf$isRoot )
      {
        if (leaf$sign == -1)
        {
          vecIsInBoxj = vecIsInBoxj & (newDataXJ[ ,leaf$jstar] <= leaf$xjstar)
        } else if (leaf$sign == 1){
          vecIsInBoxj = vecIsInBoxj & (newDataXJ[ ,leaf$jstar] > leaf$xjstar)
        }

        # We go up one step in the tree
        leaf = leaf$parent
      }
      # We fill the corresponding column
      # We don't need to fill the 0 because they are already present in the matrix
      matrixInd[which(vecIsInBoxj) , j] = 1
    }
  }

  return (matrixInd)
}


#' Construction of the matrix of estimated CKT given a matrix of indicators
#'
#' @rdname conv_treeCKT
#' @export
#'
matrixInd2matrixCKT <- function(matrixInd, newDataXI)
{
  if (nrow(newDataXI) != nrow(matrixInd)){
    stop("`newDataXI` and `matrixInd` should have the same number of rows, ",
         "as they come from the same sample.")
  }

  # We count the number of boxes in the partition
  m <- ncol(matrixInd)
  sizeI <- ncol(newDataXI)

  matrix_hat_CKT = matrix(data = NA, nrow = sizeI * (sizeI-1) / 2 , ncol = m)

  # Adding names
  if (is.null(colnames(newDataXI))){
    colnames(newDataXI) <- paste0("X", 1:sizeI)
  }
  index_pair = 1
  for (iVar in 1:(sizeI-1)){
    for (jVar in (iVar+1):sizeI){
      rownames(matrix_hat_CKT)[[index_pair]] <- paste0("CKT_", colnames(newDataXI)[[iVar]],
                                                       "_", colnames(newDataXI)[[jVar]])
      index_pair = index_pair + 1
    }
  }

  # We fill each matrix column by column
  # for each box, which correspond to each leaf of the tree
  for (j in 1:m){
    # Now we compute the conditional Kendall's tau for each couple of variables
    # given XJ in the j-th box

    index_pair = 1
    for (iVar in 1:(sizeI-1)){
      for (jVar in (iVar+1):sizeI){
        matrix_hat_CKT[index_pair, j] =
          wdm::wdm(newDataXI[which(matrixInd[, j] == 1) , iVar] ,
                   newDataXI[which(matrixInd[, j] == 1) , jVar] ,
                   method = "kendall" )
        index_pair = index_pair + 1
      }
    }
  }

  return (matrix_hat_CKT)
}


#' Construction of the matrix of estimated CKT given an estimated CKT tree
#'
#'
#' @rdname conv_treeCKT
#' @export
#'
treeCKT2matrixCKT <- function(estimatedTree, newDataXI = NULL, newDataXJ = NULL)
{
  if (is.null(newDataXI) & is.null(newDataXJ)) {
    # We get all the leaves of the tree, which corresponds to the final boxes
    leaves <- data.tree::Traverse(estimatedTree, filterFun = data.tree::isLeaf)

    matrix_hat_CKT = matrix(data.tree::Get(leaves, "CKT"),
                            ncol = length(leaves))
  } else if (!is.null(newDataXI)) {
    matrixInd = treeCKT2matrixInd(estimatedTree, newDataXJ = newDataXJ)
    matrix_hat_CKT = matrixInd2matrixCKT(matrixInd = matrixInd, newDataXI = newDataXI)
  } else {
    stop("`newDataXI` must be provided if `newDataXJ` is not null.")
  }

  rownames(matrix_hat_CKT) <- paste0("CKT_", names(estimatedTree$CKT))

  return (matrix_hat_CKT)
}
