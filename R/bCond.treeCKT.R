
#' Construct a binary tree for the modeling the conditional Kendall's tau
#'
#' This function takes in parameter two matrices of observations:
#' the first one contains the observations of \code{XI} (the conditioned variables)
#' and the second on contains the observations of \code{XJ} (the conditioning variables).
#' The goal of this procedure is to find which of the variables in \code{XJ}
#' have important influence on the dependence between the components of \code{XI},
#' (measured by the Kendall's tau).
#'
#' The object return by this function is a binary tree. Each leaf of this tree
#' correspond to one event (or, equivalently, one subset of \eqn{R^{dim(XJ)}}),
#' and the conditional Kendall's tau conditionally to it.
#'
#'
#' @param XI matrix of size n*p of observations of the conditioned variables.
#' @param XJ matrix of size n*(d-p) containing observations of the conditioning vector.
#'
#' @param minCut minimum difference in probabilities that is necessary to cut.
#' @param minProb minimum probability of being in one of the node.
#' @param minSize minimum number of observations in each node.
#'   This is an alternative to minProb and has priority over it.
#'
#' @param nPoints_xJ number of points in the grid that are considered
#' when choosing the point for splitting the tree.
#'
#' @param type.quantile way of computing the quantiles,
#' see \code{stats::\link[stats]{quantile}()}.
#'
#' @param verbose control the text output of the procedure.
#' If \code{verbose = 0}, suppress all output.
#' If \code{verbose = 2}, the progress of the computation
#' is printed during the computation.
#'
#'
#' @return the estimated tree using the data `XI, XJ`.
#'
#'
#' @references Derumigny, A., Fermanian, J. D., & Min, A. (2020).
#' Testing for equality between conditional copulas
#' given discretized conditioning events.
#' ArXiv preprint \href{https://arxiv.org/abs/2008.09498}{arxiv:2008.09498}.
#'
#'
#' @seealso \code{\link{bCond.simpA.CKT}} for a test of the simplifying assumption
#' that all these conditional Kendall's tau are equal.
#'
#' \code{\link{treeCKT2matrixInd}} for converting this tree to a matrix of indicators
#' of each event. \code{\link{matrixInd2matrixCKT}} for getting the matrix of estimated
#' conditional Kendall's taus for each event.
#'
#' \code{\link{CKT.estimate}} for the estimation of
#' pointwise conditional Kendall's tau,
#' i.e. assuming a continuous conditioning variable \eqn{Z}.
#'
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
#' result = bCond.treeCKT(XI = XI, XJ = XJ, minSize = 50, verbose = 2)
#'
#' # Number of observations in the first two children
#' print(length(data.tree::GetAttribute(result$children[[1]], "condObs")))
#' print(length(data.tree::GetAttribute(result$children[[2]], "condObs")))
#'
#'
#' @export
#'
bCond.treeCKT <- function(XI, XJ,
                          minCut = 0, minProb = 0.01, minSize = minProb * length(XI),
                          nPoints_xJ = 10, type.quantile = 7,
                          verbose = 2)
{
  # Doing various checks
  if(nrow(XI) != NROW(XJ)){
    stop("XI and XJ must have the same number of observations.")
  }
  if (anyNA(XI) || anyNA(XJ)){
    stop("XI and XJ must not contain missing values.")
  }
  if (ncol(XI) < 2){
    stop("XI must have at least two columns")
  }
  if (NCOL(XJ) == 1){
    XJ = matrix(XJ, ncol = 1)
  } else {
    XJ = as.matrix(XJ)
  }
  if (minSize < 4){
    warning("minSize < 4 not supported. Defaulting to minSize = 4")
    minSize = 4
  }

  # We add names on the variables if necessary
  if (is.null(colnames(XI))){
    namesXI <- paste0("X", 1:ncol(XI))
  } else {
    namesXI <- colnames(XI)
  }
  if (is.null(colnames(XJ))){
    namesXJ <- paste0("X", ncol(XI) + 1:ncol(XJ))
  } else {
    namesXJ <- colnames(XJ)
  }

  result = CutCKTMult(
    XI = XI, XJ = XJ, namesXI = namesXI, namesXJ = namesXJ,
    node = data.tree::Node$new("Init"),
    minCut = minCut, minProb = minProb, minSize = minSize,
    nPoints_xJ = nPoints_xJ, type.quantile = type.quantile,
    verbose = verbose, verbose.space = "")

  return (result)
}


#' Utility function for building the recursive tree
#' for conditional Kendall's tau of discretized conditioning events
#'
#' @param XI matrix of size n*p of observations of the conditioned variables.
#' @param XJ matrix of size n*(d-p) containing observations of the conditioning vector.
#' @param node current node of the tree ;
#'   by default, the function creates a new tree at each call.
#'
#' @param minCut minimum difference in probabilities that is necessary to cut.
#' @param minProb minimum probability of being in one of the node.
#' @param minSize minimum number of observations in each node.
#'   This is an alternative to minProb and has priority over it.
#'
#' @param nPoints_xJ number of points in the grid that are considered
#' when choosing the point for splitting the tree.
#' @param type.quantile way of computing the quantiles,
#' see \code{\link[stats::quantile]{stats::quantile()}}.
#' @param verbose control the text output of the procedure.
#' If \code{verbose = 0}, suppress all output.
#' If \code{verbose = 2}, the progress of the computation tree
#' is printed during the computation.
#' @param verbose.space used for padding space when printing the tree
#' (used only if \code{verbose = 2}).
#'
#' @noRd
#'
CutCKTMult <- function(XI , XJ, namesXI, namesXJ,
                       node = data.tree::Node$new("Init"),
                       minCut = 0, minProb = 0.01, minSize = minProb * length(XI),
                       nPoints_xJ = 10, type.quantile = 7,
                       verbose = 2, verbose.space = "")
{
  # Number of conditioned variables
  sizeI = ncol(XI)
  # Number of conditioning variables
  sizeJ = ncol(XJ)

  # We first compute the Kendall's tau at this point of the tree
  node$CKT = rep(NA, sizeI * (sizeI - 1 )/2 )
  index_pair = 1
  for (iVar in 1:(sizeI-1)){
    for (jVar in (iVar+1):sizeI){
      node$CKT[index_pair] = wdm::wdm(XI[,iVar], XI[,jVar], method = "kendall")
      names(node$CKT)[index_pair] = paste0(namesXI[c(iVar, jVar)], collapse = "_")
      index_pair = index_pair + 1
    }
  }
  # If this is the root, we indicate it in the list of observations
  if (data.tree::isRoot(node)){node$condObs = 1:length(XI[,1])}
  # In any case, we add the information about the size of the conditional dataset
  node$size = length(node$condObs)


  # Computation of the criterion and optimization -------------------------

  # Creation of arrays to store the conditional Kendall's taus lower and greater
  #
  # The value arrayCKTl[i,j,k] corresponds to
  # the conditional Kendall's tau of the i-th couple of conditioned variable
  # conditionally to the j-th conditioning variable being lower
  # than its quantile at level k/(nPoints_xJ+1).
  # Conversely, arrayCKTl[i,j,k] corresponds to the conditioning
  # with respect to the "greater than" event.

  arrayCKTl = array(dim = c(sizeI * (sizeI - 1 )/2 , sizeJ , nPoints_xJ) )
  arrayCKTg = array(dim = c(sizeI * (sizeI - 1 )/2 , sizeJ , nPoints_xJ) )

  # For each couple of conditioned variables
  index_pair = 1
  for (iVar1 in 1:(sizeI-1)){
    for (iVar2 in (iVar1+1):sizeI){
      # for each conditioning variable
      for (jVar in (1:sizeJ)) {
        # For each quantile that we consider
        for (k_xj in 1:nPoints_xJ) {
          xj = stats::quantile(x = as.numeric(XJ[,jVar]),
                               probs = k_xj/(nPoints_xJ+1) , type = type.quantile)
          whichij = which(XJ[,jVar] <= xj)

          arrayCKTl[index_pair, jVar, k_xj] =
            wdm::wdm(XI[whichij, iVar1], XI[whichij, iVar2], method = "kendall")

          arrayCKTg[index_pair, jVar, k_xj] =
            wdm::wdm(XI[-whichij, iVar1], XI[-whichij, iVar2], method = "kendall")
        }
      }
      # We switch to a new pair of conditioned variables
      index_pair = index_pair + 1
    }
  }

  # Our criteria is the difference between the two conditional Kendall's tau
  arrayCrit = abs(arrayCKTl - arrayCKTg)

  # We choose one that maximize the criteria.
  # `solution` is a 3 dimensional vector of indexes determining
  # the position of the optimal criterion in the 3-dimensional array arrayCrit
  solution = arrayInd(which.max(arrayCrit) , dim(arrayCrit))[1,]

  nameXIstar = names(node$CKT)[solution[1]]
  jstar = solution[2]
  nameJstar = namesXJ[jstar]
  xjstar = stats::quantile(x = XJ[,jstar], probs = solution[3]/(nPoints_xJ+1) ,
                           type = type.quantile)
  whichl = which(XJ[,jstar] <= xjstar)
  whichg = which(XJ[,jstar] > xjstar)


  # Creation of new leaves and end of the function ----------------------------------------


  # If one of our stopping criteria is true
  if ( length(whichl) < minSize ||
       length(whichg) < minSize ||
       arrayCrit[solution[1], solution[2], solution[3]] < minCut){
    # Then we don't do the cut
    if (data.tree::isRoot(node) & verbose > 0) {warning("No significant partionning found.")}
    return (node)
  } else {
    # We do the cut on variable j at point xjstar

    # Creation of the two children node
    nameChildl = paste0(nameXIstar, " ; ", nameJstar, " ; ", prettyNum(xjstar), " -")
    nameChildg = paste0(nameXIstar, " ; ", nameJstar, " ; ", prettyNum(xjstar), " +")

    nodel = node$AddChild(nameChildl)
    nodeg = node$AddChild(nameChildg)
    nodel$namejstar = nameJstar
    nodeg$namejstar = nameJstar
    nodel$jstar = jstar
    nodeg$jstar = jstar
    nodel$xjstar = xjstar
    nodeg$xjstar = xjstar
    nodel$sign = -1
    nodeg$sign = +1
    nodel$condObs = node$condObs[whichl]
    nodeg$condObs = node$condObs[whichg]

    # We call now CutCKTMult on each of the children of the tree
    if (verbose > 1) {
      cat(paste0(verbose.space, "|-- ", nodel$name,
                 "  ", prettyNum(arrayCKTl[solution[1], solution[2], solution[3]]), "\n"))
    }
    CutCKTMult(XI = XI[whichl, ], XJ = data.frame(XJ[whichl,]),
               namesXI = namesXI, namesXJ = namesXJ,
               node = nodel,
               minCut = minCut, minProb = minProb, minSize = minSize,
               nPoints_xJ = nPoints_xJ, type.quantile = type.quantile,
               verbose = verbose, verbose.space = paste0(verbose.space, "|   "))
    if (verbose > 1) {
      cat(paste0(verbose.space, "o-- ", nodeg$name,
                 "  ", prettyNum(arrayCKTg[solution[1], solution[2], solution[3]]), "\n"))
    }
    CutCKTMult(XI = XI[whichg, ], XJ = data.frame(XJ[whichg,]),
               namesXI = namesXI, namesXJ = namesXJ,
               node = nodeg,
               minCut = minCut, minProb = minProb, minSize = minSize,
               nPoints_xJ = nPoints_xJ, type.quantile = type.quantile,
               verbose = verbose, verbose.space = paste0(verbose.space, "    "))
  }

  return(node)
}

