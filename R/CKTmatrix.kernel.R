

#' Estimate the conditional Kendall's tau matrix
#' at different conditioning points
#'
#' Assume that we are interested in a random vector \eqn{(X, Z)},
#' where \eqn{X} is of dimension \eqn{d > 2} and \eqn{Z} is of dimension \eqn{1}.
#' We want to estimate the dependence across the elements of the conditioned vector \eqn{X}
#' given \eqn{Z=z}.
#' This function takes in parameter observations of \eqn{(X,Z)}
#' and returns kernel-based estimators of \deqn{\tau_{i,j | Z=zk}}
#' which is the conditional Kendall's tau between \eqn{X_i} and \eqn{X_j}
#' given to \eqn{Z=zk}, for every conditioning point \eqn{zk} in \code{gridZ}.
#'
#' If the conditional Kendall's tau matrix has a block structure,
#' then improved estimation is possible by averaging over the kernel-based estimators of
#' pairwise conditional Kendall's taus.
#' Groups of variables composing the same blocks can be defined
#' using the parameter \code{blockStructure}, and the averaging can be set on using
#' the parameter \code{averaging=all}, or \code{averaging=diag}
#' for faster estimation by averaging only over diagonal elements of each block.
#'
#' @param dataMatrix a matrix of size \code{(n,d)} containing \code{n} observations of a
#' \code{d}-dimensional random vector \eqn{X}.
#'
#' @param observedZ vector of observed points of a conditioning variable \eqn{Z}.
#' It must have the same length as the number of rows of \code{dataMatrix}.
#'
#' @param h bandwidth. It can be a real, in this case the same \code{h}
#' will be used for every element of \code{gridZ}.
#' If \code{h} is a vector then its elements are recycled to match the length of
#' \code{gridZ}.
#'
#' @param gridZ points at which the conditional Kendall's tau is computed.
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#' @param kernel.name name of the kernel used for smoothing.
#' Possible choices are: \code{"Gaussian"} (Gaussian kernel)
#' and \code{"Epa"} (Epanechnikov kernel).
#'
#' @param averaging type of averaging used for fast estimation.
#' Possible choices are \itemize{
#'   \item \code{no}: no averaging;
#'   \item \code{all}: averaging all Kendall's taus in each block;
#'   \item \code{diag}: averaging along diagonal blocks elements.
#' }
#'
#' @param blockStructure list of vectors.
#' Each vector corresponds to one group of variables
#' and contains the indexes of the variables that belongs to this group.
#' \code{blockStructure} must be a partition of \code{1:d},
#' where \code{d} is the number of columns in \code{dataMatrix}.
#'
#'
#' @return array with dimensions depending on \code{averaging}:
#' \itemize{
#'   \item If \code{averaging = "no"}:
#'   it returns an array of dimensions \code{(n, n, length(gridZ))},
#'   containing the estimated conditional Kendall's tau matrix given \eqn{Z = z}.
#'   Here, \code{n} is the number of rows in \code{dataMatrix}.
#'
#'   \item If \code{averaging = "all"} or \code{"diag"}:
#'   it returns an array of dimensions
#'   \code{(length(blockStructure), length(blockStructure), length(gridZ))},
#'   containing the block estimates of the conditional Kendall's tau given \eqn{Z = z}
#'   with ones on the diagonal.
#' }
#'
#' @seealso \code{\link{CKT.kernel}} for kernel-based estimation of conditional Kendall's tau
#' between two variables (i.e. the equivalent of this function
#' when \eqn{X} is bivariate and \code{d=2}).
#' \code{ElliptCopulas::\link[ElliptCopulas]{KTMatrixEst}()} for the fast estimation
#' of Kendall's tau matrix in the unconditional case (i.e., without Z and without smoothing).
#'
#' @examples
#'
#' # Data simulation
#' n = 100
#' Z = runif(n)
#' d = 5
#' CKT_11 = 0.8
#' CKT_22 = 0.9
#' CKT_12 = 0.1 + 0.5 * cos(pi * Z)
#' data_X = matrix(nrow = n, ncol = d)
#' for (i in 1:n){
#'   CKT_matrix = matrix(data =
#'     c(  1      , CKT_11   , CKT_11   , CKT_12[i], CKT_12[i] ,
#'       CKT_11   ,   1      , CKT_11   , CKT_12[i], CKT_12[i] ,
#'       CKT_11   , CKT_11   ,    1     , CKT_12[i], CKT_12[i] ,
#'       CKT_12[i], CKT_12[i], CKT_12[i],   1      , CKT_22    ,
#'       CKT_12[i], CKT_12[i], CKT_12[i], CKT_22   ,   1
#'       ) ,
#'      nrow = 5, ncol = 5)
#'   sigma = sin(pi * CKT_matrix/2)
#'   data_X[i, ] = mvtnorm::rmvnorm(n = 1, sigma = sigma)
#' }
#' plot(as.data.frame.matrix(data_X))
#'
#' # Estimation of CKT matrix
#' h = 1.06 * sd(Z) * n^{-1/5}
#' gridZ = c(0.2, 0.8)
#' estMatrixAll <- CKTmatrix.kernel(
#'   dataMatrix = data_X, observedZ = Z, gridZ = gridZ, h = h)
#' # Averaging estimator
#' estMatrixAve <- CKTmatrix.kernel(
#'   dataMatrix = data_X, observedZ = Z, gridZ = gridZ,
#'   averaging = "diag", blockStructure = list(1:3,4:5), h = h)
#'
#' # The estimated CKT matrix conditionally to Z=0.2 is:
#' estMatrixAll[ , , 1]
#' # Using the averaging estimator,
#' # the estimated CKT between the first group (variables 1 to 3)
#' # and the second group (variables 4 and 5) is
#' estMatrixAve[1, 2, 1]
#'
#' # True value (of CKT between variables in block 1 and 2 given Z = 0.2):
#' 0.1 + 0.5 * cos(pi * 0.2)
#'
#'
#' @author Rutger van der Spek, Alexis Derumigny
#'
#' @export
#'
CKTmatrix.kernel <- function(dataMatrix, observedZ, gridZ,
                             averaging = "no", blockStructure = NULL,
                             h, kernel.name = "Epa",
                             typeEstCKT = 4)
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


  if(averaging == "no")
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
  } else if(averaging == "diag")
  {
    if(is.null(blockStructure))
    {
      stop(paste("blockStructure not specified, when averaging = ", averaging))
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
  } else if (averaging == "all")
  {
    if(is.null(blockStructure))
    {
      stop(paste("blockStructure not specified, when averaging = ", averaging))
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
  } else
  {
    stop("'averaging' must be one of: 'no', 'all' or 'diag'.")
  }

  return(estimate)
}


