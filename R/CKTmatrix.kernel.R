

#' Estimate the (averaging) conditional Kendall's tau matrix
#' at different points
#'
#' @param dataMatrix n x d matrix of d-dimensional multivariate random variable
#' of n timepoints.
#'
#' @param observedZ vector of observed points of Z.
#' It must have the same length as the number of rows of \code{dataMatrix}.
#'
#' @param h bandwidth. It can be a real, in this case the same \code{h}
#' will be used for every element of \code{gridZ}.
#' If \code{h}is a vector then its elements are recycled to match the length of
#' \code{gridZ}.
#'
#' @param gridZ points at which the conditional Kendall's tau is computed.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#'
#' @param typeMatrixCKT name of the matrix estimator used. Possible choices are
#' "all" (no averaging),
#' "aveDiag" (averaging along diagonal block elements) and
#' "aveAll" (averaging all CKT's within blocks)
#' @param blockStructure list of groups. A group is a vector with
#' variable numbers. \code{blockStructure} must be a partition of 1:d, where d is
#' the amount of columns in \code{dataMatrix}.
#'
#'
#' @return array with dimensions depending on \code{typeMatrixCKT}.
#' If \code{typeMatrixCKT} = "all", it returns an array of dimension
#' n x n x \code{length(gridZ)}, giving the conditional Kendall's tau matrix
#' given Z = z. Here, n is the number of rows in \code{dataMatrix}.
#' If \code{typeMatrixCKT} = "aveDiag" or "aveAll" the function returns an array
#' of dimension \code{length(blockStructure)} x
#' \code{length(blockStructure)} x  \code{length(gridZ)}, giving the block
#' estimates of the conditional Kendall's tau given Z = z with ones on the
#' diagonal.
#'
#' @export
#'
#' @author Rutger van der Spek, Alexis Derumigny
#'
CKTmatrix.kernel <- function(dataMatrix, observedZ, gridZ,
                             typeMatrixCKT = "all", blockStructure = NULL,
                             h, kernel.name = "Epa",
                             typeEstCKT = 4
)
{

  d = ncol( dataMatrix )
  n = nrow( dataMatrix )
  nz = length( gridZ )

  if(length(observedZ) != n)
  {
    stop("The length of observedZ and the number of rows in dataMatrix must be equal.")
  }

  arrayWeights = array(data = NA , dim = c(n ,n , nz) )
  matrixWeights = matrix(data = NA, nrow = n, ncol = nz)
  for(i in 1:length(gridZ))
  {
    matrixWeights[,i] = computeWeights.univariate(vectorZ = observedZ,
                                                  h = h,
                                                  pointZ = gridZ[i],
                                                  kernel.name = kernel.name,
                                                  normalization = TRUE)
    arrayWeights[,,i] = matrixWeights[,i] %*% t(matrixWeights[,i])
  }


  if(typeMatrixCKT == "all")
  {
    estimate = array(data = 1 , dim = c(d , d , nz))
    for(j1 in 2:d)
    {
      for(j2 in 1:(j1-1))
      {
        vectorX1 = dataMatrix[ , j1]
        vectorX2 = dataMatrix[ , j2]
        matrixSigns = computeMatrixSignPairs(vectorX1 = vectorX1,
                                             vectorX2 = vectorX2,
                                             typeEstCKT = typeEstCKT)
        estimate[j1,j2,] = apply(X = arrayWeights , MARGIN = 3 ,
                                 FUN = function(x){return(sum(x*matrixSigns))})
        estimate[j2,j1,] = estimate[j1,j2,]
      }
    }
  }


  else if(typeMatrixCKT == "aveDiag")
  {
    if(is.null(blockStructure))
    {
      stop(paste("blockStructure not specified, when typeMatrixCKT = ", typeMatrixCKT))
    }
    if( all.equal( sort(unname(unlist(blockStructure))) , 1:d ) != TRUE )
    {
      stop("blockStructure must be a partition.")
    }
    totalGroups = length(blockStructure)
    estimate = array(data = 1 , dim = c(totalGroups , totalGroups , nz))
    for (g1 in 2:totalGroups)
    {
      for (g2 in 1:(g1-1))
      {

        diagSize = min(length(blockStructure[[g1]]),
                       length(blockStructure[[g2]]) )
        matrixBlockCKT = matrix(NA, nrow = diagSize , ncol = nz)
        for (j in 1:diagSize)
        {
          vectorX1 = dataMatrix[ , blockStructure[[g1]][j] ]
          vectorX2 = dataMatrix[ , blockStructure[[g2]][j] ]
          matrixSigns = computeMatrixSignPairs(vectorX1 = vectorX1,
                                               vectorX2 = vectorX2,
                                               typeEstCKT = typeEstCKT)
          matrixBlockCKT[j,] = apply(X = arrayWeights , MARGIN = 3,
                                     FUN = function(x)
                                     {return(sum(x * matrixSigns))})
        }
        blockCKT = apply(X = matrixBlockCKT , MARGIN = 2 , mean)
        estimate[g1 , g2 ,] = blockCKT
        estimate[g2 , g1 ,] = blockCKT

      }
    }
  }


  else if(typeMatrixCKT == "aveAll")
  {
    if(is.null(blockStructure))
    {
      stop(paste("blockStructure not specified, when typeMatrixCKT = ", typeMatrixCKT))
    }
    else if( all.equal( sort(unname(unlist(blockStructure))) , 1:d ) != TRUE )
    {
      stop("blockStructure must be a partition.")
    }
    totalGroups = length(blockStructure)
    estimate = array(data = 1 , dim = c(totalGroups , totalGroups , nz))
    for (g1 in 2:totalGroups)
    {
      for (g2 in 1:(g1-1))
      {
        arrayBlockCKT = array(NA, dim = c(length(blockStructure[[g1]]),
                                          length(blockStructure[[g2]]), nz) )
        for (j1 in 1:length(blockStructure[[g1]]) )
        {
          for (j2 in 1:length(blockStructure[[g2]]) )
          {
            vectorX1 = dataMatrix[ , blockStructure[[g1]][j1] ]
            vectorX2 = dataMatrix[ , blockStructure[[g2]][j2] ]
            matrixSigns = computeMatrixSignPairs(vectorX1 = vectorX1,
                                                 vectorX2 = vectorX2,
                                                 typeEstCKT = typeEstCKT)
            arrayBlockCKT[j1,j2,] = apply(X = arrayWeights , MARGIN = 3,
                                          FUN = function(x)
                                          {return(sum(x * matrixSigns))})
          }
        }
        blockCKT = apply(X = arrayBlockCKT , MARGIN = 3 , mean)
        estimate[g1,g2,] = blockCKT
        estimate[g2,g1,] = blockCKT

      }
    }
  }

  else
  {
    stop("typeMatrixCKT must be: 'all', 'aveDiag' or 'aveAll'.")
  }

  return(estimate)
}


