

#' Semiparametric testing of the simplifying assumption
#'
#' @param X1 observed vector of the first conditionned variable
#' @param X2 observed vector of the second conditionned variable
#' @param X3 observed vector of the conditioning variable
#'
#' @param testStat name of the test statistic to be used
#' @param typeBoot the type of bootstrap to be used
#' @param nBootstrap number of bootstrap replications
#' @param h the bandwidth used for kernel smoothing
#' @param kernel.name the name of the kernel
#' @param truncVal the value of truncation for the integral,
#' i.e. the integrals are computed from \code{truncVal} to \code{1-truncVal}
#' instead of from 0 to 1.
#' @param numericalInt parameters to be given to
#' \code{statmod::\link[statmod]{gauss.quad}}, including the number of
#' quadrature points and the type of interpolation.
#'
#' @param family the chosen family of copulas
#'   (see the documentation of the class \code{VineCopula::\link[VineCopula]{BiCop}()}
#'   for the available families).
#'
#' @return a list containing
#' \itemize{
#'     \item \code{true_stat}: the value of the test statistic computed on the whole sample
#'     \item \code{vect_statB}: a vector of length \code{nBootstrap} containing the bootstrapped
#'     test statistics.
#'     \item \code{p_val}: the p-value of the test.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 800
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1], mean = Z)
#' X2 = qnorm(simCopula[,2], mean = - Z)
#'
#' result <- simpA.param(
#'    X1 = X1, X2 = X2, X3 = Z, family = 1,
#'    h = 0.03, kernel.name = "Epanechnikov")
#' print(result$p_val)
#'
#' \donttest{
#' set.seed(1)
#' N = 500
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 0.8
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1], mean = Z)
#' X2 = qnorm(simCopula[,2], mean = - Z)
#'
#' result <- simpA.param(
#'    X1 = X1, X2 = X2, X3 = Z, family = 1,
#'    h = 0.08, kernel.name = "Epanechnikov", nBootstrap = 30)
#' print(result$p_val)
#' }
#'
#' @export
#'
simpA.param <- function(
  X1, X2, X3, family, testStat = "T2c", typeBoot = "boot.NP", h,
  nBootstrap = 100,
  kernel.name = "Epanechnikov", truncVal = h,
  numericalInt = list(kind = "legendre", nGrid = 10) )
{
  stopifnot(length(X1) == length(X2))
  stopifnot(length(X1) == length(X3))
  n <- length(X1)

  nGrid = numericalInt$nGrid
  grid <- statmod::gauss.quad(n = nGrid, kind = numericalInt$kind)
  # Change of range to be on [h, 1-h]
  grid$nodes <- grid$nodes * (1/2 - h) + 1/2
  grid$weights <- grid$weights / 2

  switch(
    typeBoot,

    "boot.NP" = {
      FUN_boot <- boot.NP
      nStar <- 1
    },

    "boot.pseudoInd" = {
      FUN_boot <- boot.pseudoInd
      nStar <- 2
    },

    "boot.pseudoInd.sameX3" = {
      FUN_boot <- boot.pseudoInd.sameX3
      nStar <- 2
    },

    "boot.pseudoNP" = {
      FUN_boot <- boot.pseudoNP
      nStar <- 1
    },

    "boot.cond" = {
      FUN_boot <- boot.cond
      nStar <- 1
    },

    "boot.paramInd" = {
      FUN_boot <- boot.paramInd
      nStar <- 2
    },

    "boot.paramCond" = {
      FUN_boot <- boot.paramCond
      nStar <- 1
    },

    stop("Unknown bootstrap method name.")
  )

  switch(
    testStat,

    "T2c" = {

      FUN_trueStat <- testStat_T2c
      if (nStar == 1){
        FUN_stat_st <- testStat_T2c_boot1st
      } else if (nStar == 2){
        FUN_stat_st <- testStat_T2c_boot2st
      }
    },

    stop("Unknown test statistic name. Possible choice are: ",
         "'T2c'.")
  )

  env <- environment()
  FUN_boot(env = env, FUN_trueStat = FUN_trueStat, FUN_stat_st = FUN_stat_st)

  return ( list(true_stat = env$true_stat ,
                vect_statB = env$vect_statB ,
                p_val = env$p_val) )
}

