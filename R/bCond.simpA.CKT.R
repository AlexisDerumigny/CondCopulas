

#' Function for testing the simplifying assumption with
#' data-driven box-type conditioning events
#'
#' This function takes in parameter the matrix of (observations) of the
#' conditioned variables and either \code{matrixInd}, a matrix of indicator variables
#' describing which events occur for which observations
#'
#' @param matrixInd a matrix of indexes of size (n, N.boxes) describing
#' for each observation i to which box ( = event) it belongs.
#'
#' If it is \code{NULL}, then a tree will be estimated to provide relevant boxes
#' (by using \code{\link{bCond.treeCKT}()})
#' and then converting to a \code{matrixInd} by \code{\link{treeCKT2matrixInd}()}.
#'
#' @inheritParams bCond.treeCKT
#'
#' @param methodTree method for constructing the tree
#' \itemize{
#' \item \code{doSplit} some part of the data is used for constructing the tree
#'   and the other part for constructing the test statistic
#'   using the boxes defined by the estimated tree.
#'   The share of the data used for construction the tree is controlled by
#'   the parameter \code{propTree}.
#' \item \code{noSplit} all of the data is used for
#'   both the tree and the test statistic on it.
#'   Note that p-values obtained by this method have an upward bias
#'   due to the lack of independence between these two steps.
#' }
#' Only used if \code{matrixInd} is not provided.
#'
#' @param propTree share of observations used to build the tree
#' (the rest of the observations are used for the computation of the p-value).
#' Only used if \code{matrixInd} is not provided.
#'
#' @param methodPvalue method for computing the p-value \itemize{
#' \item \code{covMatrix} by computation of the covariance matrix of
#'   the random vector \eqn{(\tau_{i,k|X_J \in A_j}, 1\leq,i,k\leq p, 1\leq j \leq m)}.
#' \item \code{bootNP} by the usual non-parametric bootstrap
#' \item \code{bootInd} by the independent bootstrap
#' }
#'
#' @param nBootstrap number of bootstrap replications
#' (Only used if \code{methodPvalue} is not \code{covMatrix}).
#'
#'
#' @return a list with the following components \itemize{
#' \item \code{p.value} the estimated p-value.
#' \item \code{stat} the test statistic.
#' \item \code{treeCKT} the estimated tree if \code{matrixInd} is not provided.
#' \item \code{vec_statB} the vector of bootstrapped statistics
#' if \code{methodPvalue} is not \code{covMatrix}.
#' }
#'
#' @references Derumigny, A., Fermanian, J. D., & Min, A. (2022).
#' Testing for equality between conditional copulas
#' given discretized conditioning events.
#' Canadian Journal of Statistics.
#' \doi{10.1002/cjs.11742}
#'
#' Derumigny, A., & Fermanian, J. D. (2022)
#' Conditional empirical copula processes and generalized dependence measures
#' Electronic Journal of Statistics, 16(2), 5692-5719.
#' \doi{10.1214/22-EJS2075}
#'
#'
#' @seealso \code{\link{bCond.simpA.param}} for a test of this simplifying assumption
#' in a parametric framework.
#'
#' \code{\link{bCond.treeCKT}} provides the binary tree that is used in this function
#' (if \code{matrixInd} is not provided).
#'
#' Tests of the simplifying assumption for conditional copulas with a continuous
#' conditioning variable:
#' \itemize{
#'   \item \code{\link{simpA.NP}} in a nonparametric setting
#'   \item \code{\link{simpA.param}} in a (semi)parametric setting,
#'   where the conditional copula belongs to a parametric family,
#'   but the conditional margins are estimated arbitrarily through
#'   kernel smoothing
#'   \item \code{\link{simpA.kendallReg}}: test based on the constancy of
#'   conditional Kendall's tau
#' }
#'
#'
#' @author Alexis Derumigny, Jean-David Fermanian and Aleksey Min
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
#' result = bCond.simpA.CKT(XI = XI, XJ = XJ, minSize = 10, verbose = 2,
#'                          methodTree = "doSplit", nBootstrap = 4)
#' print(result$p.value)
#' result2 = bCond.simpA.CKT(XI = XI, XJ = XJ, minSize = 10, verbose = 2,
#'                           methodTree = "noSplit", nBootstrap = 4)
#' print(result2$p.value)
#'
#'
#' @export
#'
bCond.simpA.CKT <- function(XI, XJ = NULL, matrixInd = NULL,
                            minCut = 0,
                            minProb = 0.01, minSize = minProb * nrow(XI),
                            nPoints_xJ = 10, type.quantile = 7,
                            verbose = 2,
                            methodTree = "doSplit", propTree = 0.5,
                            methodPvalue = "bootNP", nBootstrap = 100)
{
  n = nrow(XI)

  if (is.null(matrixInd))
  {
    if (is.null(XJ)){
      stop("Exactly one of `XJ` or `matrixInd` should not be `NULL`. ",
           "Currently, both are `NULL`.")
    }
    if (nrow(XI) != nrow(XJ)){
      stop("XI and XJ should have the same number of rows, ",
           "equal to the number of observations in the dataset.")
    }

    # First part: construction of the tree -----------------------

    if (methodTree == "doSplit") {
      # We select some part of the observations for building the tree
      sampleForTree = sample(1:n, size = propTree * n)
      sampleForStat = (1:n)[-sampleForTree]
      XItree = XI[sampleForTree, ]
      XJtree = data.frame(XJ[sampleForTree, ])

      treeCKT = bCond.treeCKT(
        XI = XItree, XJ = XJtree,
        minCut = minCut, minSize = minSize, minProb = minProb,
        nPoints_xJ = nPoints_xJ, type.quantile = type.quantile,
        verbose = verbose)

      XIstat = XI[sampleForStat, ]
      XJstat = data.frame(XJ[sampleForStat, ])
    } else {
      sampleForStat = 1:n
      XIstat = XI
      XJstat = NULL

      treeCKT = bCond.treeCKT(
        XI = XI, XJ = XJ,
        minCut = minCut, minSize = minSize, minProb = minProb,
        nPoints_xJ = nPoints_xJ, type.quantile = type.quantile,
        verbose = verbose)
    }

    # If we cannot find any relevant partition, we stop immediately
    if (treeCKT$height == 1)
    {
      return(list("p.value" = "",
                  "stat" = "",
                  "matrix_hat_CKT" = as.matrix(treeCKT$CKT, ncol = 1),
                  "treeCKT" = treeCKT ))
    }

    # Second part: construction of the matrix of indices -----------------------------
    matrixInd = treeCKT2matrixInd(
      estimatedTree = treeCKT, newDataXJ = XJstat)

    # Third part: construction of the matrix of estimated CKT ------------------------

    if (methodTree == "doSplit")
    {
      matrix_hat_CKT = matrixInd2matrixCKT(
        matrixInd = matrixInd, newDataXI = XIstat)
    } else {
      matrix_hat_CKT = treeCKT2matrixCKT(
        estimatedTree = treeCKT, newDataXI = XIstat, newDataXJ = XJstat)
    }
  } else {
    if (!is.null(XJ)){
      stop("Exactly one of `XJ` or `matrixInd` should not be `NULL`. ",
           "Currently, none are `NULL`.")
    }
    if (nrow(matrixInd) != nrow(XJ)){
      stop("matrixInd and XJ should have the same number of rows, ",
           "equal to the number of observations in the dataset.")
    }

    # If a matrixInd has been given by the user
    # We just need to compute the matrix of estimated CKT
    matrix_hat_CKT = matrixInd2matrixCKT(
      matrixInd = matrixInd, newDataXI = XI)

    # We initialize other variables that will be used for the computation of the p-values.
    XIstat = XI
    treeCKT = NULL
  }

  # Final part: computation of the p-value ---------------------------------------
  # using either the estimated covariance matrix
  # or using bootstrap techniques

  result = bCond.simpA.CKT.computepvalue(
    XIstat = XIstat, matrixInd = matrixInd, matrix_hat_CKT = matrix_hat_CKT,
    methodPvalue = methodPvalue, nBootstrap = nBootstrap)

  result$treeCKT = treeCKT

  return (result)
}


#' Computation of the p-value for tests of the simplifying assumption
#' based on the constancy of the conditional Kendall's tau across boxes
#'
#' @param XIstat matrix of observations of the conditioned variables.
#' @param matrixInd matrix of indicators of each of the events.
#' @param methodPvalue method for computing the p-values.
#' @param nBootstrap number of boostrap replications.
#'
#' @noRd
#'
bCond.simpA.CKT.computepvalue <- function(
  XIstat, matrixInd, matrix_hat_CKT,
  methodPvalue, nBootstrap = 100)
{
  nSample = nrow(XIstat)   # should be equal to nrow(matrixInd)

  switch (
    methodPvalue,

    "covMatrix" = {
      # computation of diverse elements for the covariance matrix
      statCovMat = fun.mult.norm.prop.4(
        xI = XIstat ,
        ind = matrixInd ,
        hat_CKT = matrix_hat_CKT ,
        tau.cond.true = 0)

      # computation of the test statistic and the p-value
      test = fun.Wald.Stat.gen.dim(R.Object = statCovMat)

      # We return the elements as a list

      return(list("p.value" = test$p.value,
                  "stat" = test$stat,
                  "covMatrix" = test$covMatrix,
                  "matrix_hat_CKT" = matrix_hat_CKT))
    },

    "bootNP" = {
      true_stat = fun.Maximum.Stat.gen.dim(hat_CKT = matrix_hat_CKT, tau.cond.true = 0)

      vec_statB = rep(NA, nBootstrap)

      for (iBootstrap in 1:nBootstrap){
        # We construct the new sample
        sample_st <- sample(x = 1:nSample, size = nSample, replace = TRUE)

        # We compute the corresponding Kendall's taus
        matrix_hat_CKT_st <- matrixInd2matrixCKT(
          matrixInd = matrixInd[sample_st, ],
          newDataXI = XIstat[sample_st, ])

        vec_statB[iBootstrap] <- fun.Maximum.Stat.gen.dim(hat_CKT = matrix_hat_CKT_st,
                                                          tau.cond.true = matrix_hat_CKT)
      }
      p.value = mean(as.numeric(true_stat < vec_statB), na.rm = TRUE)

      return(list("p.value" = p.value,
                  "stat" = true_stat,
                  "vec_statB" = vec_statB,
                  "matrix_hat_CKT" = matrix_hat_CKT))
    },

    "bootInd" = {
      true_stat = fun.Maximum.Stat.gen.dim(hat_CKT = matrix_hat_CKT, tau.cond.true = 0)

      vec_statB = rep(NA, nBootstrap)

      for (iBootstrap in 1:nBootstrap){
        # We construct the new sample
        sample_st_I <- sample(x = 1:nSample, size = nSample, replace = TRUE)
        sample_st_J <- sample(x = 1:nSample, size = nSample, replace = TRUE)

        # We compute the corresponding Kendall's taus
        matrix_hat_CKT_st <- matrixInd2matrixCKT(
          matrixInd = matrixInd[sample_st_J, ],
          newDataXI = XIstat[sample_st_I, ])

        vec_statB[iBootstrap] <- fun.Maximum.Stat.gen.dim(
          hat_CKT = matrix_hat_CKT_st, tau.cond.true = 0)
      }
      p.value = mean(as.numeric(true_stat < vec_statB))

      return(list("p.value" = p.value,
                  "stat" = true_stat,
                  "vec_statB" = vec_statB,
                  "matrix_hat_CKT" = matrix_hat_CKT ))

    },

    {stop("Unknown bootstrap method.")}
  )

}

