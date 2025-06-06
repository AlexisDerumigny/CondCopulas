

#' Semiparametric testing of the simplifying assumption
#'
#' This function tests the ``simplifying assumption'' that a conditional
#' copula \deqn{C_{1,2|3}(u_1, u_2 | X_3 = x_3)} does not depend on the
#' value of the conditioning variable \eqn{x_3} in a semiparametric setting,
#' where the conditional copula is of the form
#' \deqn{C_{1,2|3}(u_1, u_2 | X_3 = x_3) = C_{\theta(x_3)}(u_1,u_2),}
#' for all \eqn{0 <= u_1, u_2 <= 1} and all \eqn{x_3}.
#' Here, \eqn{(C_\theta)} is a known family of copula and \eqn{\theta(x_3)}
#' is an unknown conditional dependence parameter.
#' In this setting, the simplifying assumption can be rewritten as
#' \strong{``\eqn{\theta(x_3)} does not depend on \eqn{x_3}, i.e. is a constant
#' function of \eqn{x_3}''}.
#'
#' @param X1 vector of \code{n} observations of the first conditioned variable
#'
#' @param X2 vector of \code{n} observations of the second conditioned variable
#'
#' @param X3 vector of \code{n} observations  of the conditioning variable
#'
#' @param testStat name of the test statistic to be used.
#' The only choice implemented yet is \code{'T2c'}.
#'
#' @param typeBoot the type of bootstrap to be used.
#' (see Derumigny and Fermanian, 2017, p.165).
#' Possible values are
#' \itemize{
#'    \item \code{"boot.NP"}: usual (Efron's) non-parametric bootstrap
#'
#'    \item \code{"boot.pseudoInd"}: pseudo-independent bootstrap
#'
#'    \item \code{"boot.pseudoInd.sameX3"}: pseudo-independent bootstrap
#'    without resampling on \eqn{X_3}
#'
#'    \item \code{"boot.pseudoNP"}: pseudo-non-parametric bootstrap
#'
#'    \item \code{"boot.cond"}: conditional bootstrap
#'
#'    \item \code{"boot.paramInd"}: parametric independent bootstrap
#'
#'    \item \code{"boot.paramCond"}: parametric conditional bootstrap
#' }
#'
#'
#' @param nBootstrap number of bootstrap replications
#'
#' @param h the bandwidth used for kernel smoothing
#'
#' @param kernel.name the name of the kernel
#'
#' @param truncVal the value of truncation for the integral,
#' i.e. the integrals are computed from \code{truncVal} to \code{1-truncVal}
#' instead of from 0 to 1.
#' Note that \code{truncVal} must be in the interval \eqn{[0, 0.5)},
#' i.e. \eqn{0} is allowed but not \eqn{0.5}.
#'
#' The default is \code{truncVal = NULL}, which actually means that
#' \code{truncVal = h} if \code{h < 0.5} and \code{truncVal = 0} else.
#'
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
#'     \item \code{true_stat}: the value of the test statistic
#'     computed on the whole sample
#'     \item \code{vect_statB}: a vector of length \code{nBootstrap}
#'     containing the bootstrapped test statistics.
#'     \item \code{p_val}: the p-value of the test.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#' \doi{10.1515/demo-2017-0011}
#'
#' @seealso Other tests of the simplifying assumption:
#' \itemize{
#'   \item \code{\link{simpA.NP}} in a nonparametric setting
#'   \item \code{\link{simpA.kendallReg}}: test based on the constancy of
#'   conditional Kendall's tau
#'
#'   \item the counterparts of these tests in the discrete conditioning setting:
#'   \code{\link{bCond.simpA.CKT}}
#'   (test based on conditional Kendall's tau)
#'   \code{\link{bCond.simpA.param}}
#'   (test assuming a parametric form for the conditional copula)
#' }
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 500
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1], mean = Z)
#' X2 = qnorm(simCopula[,2], mean = - Z)
#'
#' result <- simpA.param(
#'    X1 = X1, X2 = X2, X3 = Z, family = 1,
#'    h = 0.03, kernel.name = "Epanechnikov", nBootstrap = 5)
#' print(result$p_val)
#' # In practice, it is recommended to use at least nBootstrap = 100
#' # with nBootstrap = 200 being a good choice.
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
#'    h = 0.08, kernel.name = "Epanechnikov", nBootstrap = 5)
#' print(result$p_val)
#' }
#'
#' @export
#'
simpA.param <- function(
  X1, X2, X3, family, testStat = "T2c", typeBoot = "boot.NP", h,
  nBootstrap = 100,
  kernel.name = "Epanechnikov", truncVal = NULL,
  numericalInt = list(kind = "legendre", nGrid = 10) )
{
  .checkSame_nobs_X1X2X3(X1, X2, X3)
  .checkUnivX1X2X3(X1, X2, X3)

  if (is.null(truncVal)){
    if (h < 0.5){
      truncVal = h
    } else {
      truncVal = 0
    }
  } else {
    if (truncVal < 0 || truncVal >= 0.5){
      stop(errorCondition(
        message = paste0("'truncVal' must be between 0 and 1/2 (strictly). ",
                         "Here it is: ", truncVal),
        class = "InvalidInputError"
      ))
    }
  }

  n <- length(X1)

  nGrid = numericalInt$nGrid
  grid <- statmod::gauss.quad(n = nGrid, kind = numericalInt$kind)
  # Change of range to be on [truncVal , 1 - truncVal]
  grid$nodes <- grid$nodes * (1/2 - truncVal) + 1/2
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

    stop("Unknown test statistic name. Possible choice is: ",
         "'T2c' (only implemented for the moment).")
  )

  env <- environment()
  FUN_boot(env = env, FUN_trueStat = FUN_trueStat, FUN_stat_st = FUN_stat_st)

  return ( list(true_stat = env$true_stat ,
                vect_statB = env$vect_statB ,
                p_val = env$p_val) )
}

