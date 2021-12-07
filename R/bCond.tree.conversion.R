

#' Construction of the matrix of indicators given an estimated CKT tree
#'
#' @param estimatedTree the tree that has been estimated before,
#' for example by \code{\link{bCond.treeCKT}}.
#'
#' @param newDataXJ this is a matrix of size \eqn{n' \times |J|}
#' where \eqn{|J|} is the number of conditional variables used in the tree.
#' By default this is \code{NULL} meaning that
#' we return the matrix of indicators for the original data
#' (that was used to compute the \code{estimatedTree}).
#'
#' @return a matrix of size \eqn{n' \times m} which component \code{[i,j]}
#' is \eqn{1\{ X_{i,J} \in A_{j,J} \}}.
#'
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
#'
#' @param matrixInd a matrix of indexes of size (n, N.boxes) describing
#' for each observation i to which box ( = event) it belongs.
#'
#' @param newDataXI this is a matrix of size \eqn{n' \times |I|}
#' where \eqn{|I|} is the number of conditioned variables.
#' By default this is \code{NULL} meaning that
#' we return the matrix of conditional Kendall's tau for the original data
#' used to compute the \code{estimatedTree}
#'
#' @return a matrix of size \eqn{|I|(|I|-1) \times m} where each component corresponds
#' to a conditional Kendall's tau between a pair of conditional variables
#' conditionally to the conditioned variables in one of the boxes
#'
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

  # We fill each matrix column by column
  # for each box, which correspond to each leaf of the tree
  for (j in 1:m){
    # Now we compute the conditional Kendall's tau for each couple of variables
    # given XJ in the j-th box

    index_pair = 1
    for (iVar in 1:(sizeI-1)){
      for (jVar in (iVar+1):sizeI){
        matrix_hat_CKT[index_pair, j] =
          pcaPP::cor.fk( newDataXI[which(matrixInd[, j] == 1) , iVar] ,
                         newDataXI[which(matrixInd[, j] == 1) , jVar] )
        index_pair = index_pair + 1
      }
    }
  }

  return (matrix_hat_CKT)
}


#' Construction of the matrix of estimated CKT given an estimated CKT tree
#'
#' @param estimatedTree the tree that has been estimated before,
#' for example by \code{\link{bCond.treeCKT}}
#'
#' @param newDataXI this is a matrix of size \eqn{n' \times |I|}
#' where \eqn{|I|} is the number of conditioned variables.
#' @param newDataXJ this is a matrix of size \eqn{n' \times |J|}
#' where \eqn{|J|} is the number of conditional variables.
#'
#' By default both \code{newDataXI} and \code{newDataXJ} are \code{NULL}
#' meaning that we return the matrix of conditional Kendall's tau for the original data
#' used to compute the \code{estimatedTree}
#'
#' @return a matrix of size \eqn{|I|(|I|-1) \times m} where each component corresponds
#' to a conditional Kendall's tau between a pair of conditional variables
#' conditionally to the conditioned variables in one of the boxes
#'
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

  return (matrix_hat_CKT)
}
