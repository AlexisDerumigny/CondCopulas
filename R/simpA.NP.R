
#' Nonparametric testing of the simplifying assumption
#'
#' This function tests the ``simplifying assumption'' that a conditional
#' copula \deqn{C_{1,2|3}(u_1, u_2 | X_3 = x_3)} does not depend on the
#' value of the conditioning variable \eqn{x_3} in a nonparametric setting,
#' where the conditional copula is estimated by kernel smoothing.
#'
#' @param X1 vector of \code{n} observations of the first conditioned variable
#'
#' @param X2 vector of \code{n} observations of the second conditioned variable
#'
#' @param X3 vector of \code{n} observations  of the conditioning variable
#'
#' @param testStat name of the test statistic to be used.
#' Possible values are
#' \itemize{
#'    \item \code{T1_CvM_Cs3}: Equation (3) of (Derumigny & Fermanian, 2017) with
#'    the simplified copula estimated by Equation (6) and the weight
#'    \eqn{w(u_1, u_2, u_3) = \hat{F}_1(u_1) \hat{F}_2(u_2) \hat{F}_3(u_3)}.
#'
#'    \item \code{T1_CvM_Cs4}: Equation (3) of (Derumigny & Fermanian, 2017) with
#'    the simplified copula estimated by Equation (7) and the weight
#'    \eqn{w(u_1, u_2, u_3) = \hat{F}_1(u_1) \hat{F}_2(u_2) \hat{F}_3(u_3)}.
#'
#'    \item \code{T1_KS_Cs3}: Equation (4) of (Derumigny & Fermanian, 2017) with
#'    the simplified copula estimated by Equation (6).
#'
#'    \item \code{T1_KS_Cs4}: Equation (4) of (Derumigny & Fermanian, 2017) with
#'    the simplified copula estimated by Equation (7).
#'
#'    \item \code{tilde_T0_CvM}: Equation (10) of (Derumigny & Fermanian, 2017).
#'
#'    \item \code{tilde_T0_KS}: Equation (9) of (Derumigny & Fermanian, 2017).
#'
#'    \item \code{I_chi}: Equation (13) of (Derumigny & Fermanian, 2017).
#'
#'    \item \code{I_2n}: Equation (15) of (Derumigny & Fermanian, 2017).
#' }
#'
#' @param typeBoot the type of bootstrap to be used
#' (see Derumigny and Fermanian, 2017, p.165).
#' Possible values are
#' \itemize{
#'    \item \code{boot.NP}: usual (Efron's) non-parametric bootstrap
#'
#'    \item \code{boot.pseudoInd}: pseudo-independent bootstrap
#'
#'    \item \code{boot.pseudoInd.sameX3}: pseudo-independent bootstrap
#'    without resampling on \eqn{X_3}
#'
#'    \item \code{boot.pseudoNP}: pseudo-non-parametric bootstrap
#'
#'    \item \code{boot.cond}: conditional bootstrap
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
#'
#' @param numericalInt parameters to be given to
#' \code{statmod::\link[statmod]{gauss.quad}}, including the number of
#' quadrature points and the type of interpolation.
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
#'   \item \code{\link{simpA.param}} in a (semi)parametric setting,
#'   where the conditional copula belongs to a parametric family,
#'   but the conditional margins are estimated arbitrarily through
#'   kernel smoothing
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
#' result <- simpA.NP(
#'    X1 = X1, X2 = X2, X3 = Z,
#'    testStat = "I_chi", typeBoot = "boot.pseudoInd",
#'    h = 0.03, kernel.name = "Epanechnikov", nBootstrap = 10)
#'
#' # In practice, it is recommended to use at least nBootstrap = 100
#' # with nBootstrap = 200 being a good choice.
#'
#' print(result$p_val)
#'
#' set.seed(1)
#' N = 500
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 0.8
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1], mean = Z)
#' X2 = qnorm(simCopula[,2], mean = - Z)
#'
#' result <- simpA.NP(
#'    X1 = X1, X2 = X2, X3 = Z,
#'    testStat = "I_chi", typeBoot = "boot.pseudoInd",
#'    h = 0.08, kernel.name = "Epanechnikov", nBootstrap = 10)
#' print(result$p_val)
#'
#' @export
#'
simpA.NP <- function(
  X1, X2, X3, testStat, typeBoot = "bootNP", h,
  nBootstrap = 100,
  kernel.name = "Epanechnikov", truncVal = h,
  numericalInt = list(kind = "legendre", nGrid = 10))
{
  stopifnot(length(X1) == length(X2))
  stopifnot(length(X1) == length(X3))
  n <- length(X1)

  if (testStat %in% c("T1_CvM_Cs3", "T1_CvM_Cs4", "tilde_T0_CvM",
                      "T1_KS_Cs3", "T1_KS_Cs4", "tilde_T0_KS")){
    nGrid = numericalInt$nGrid
    grid <- statmod::gauss.quad(n = nGrid, kind = numericalInt$kind)
    # Change of range to be on [0,1]
    grid$nodes <- grid$nodes * (1/2 - h) + 1/2
    grid$weights <- grid$weights / 2
  }

  if (testStat %in% c("T1_CvM_Cs3", "T1_CvM_Cs4", "tilde_T0_CvM",
                      "T1_KS_Cs3", "T1_KS_Cs4", "tilde_T0_KS") &
      typeBoot %in% c("boot.pseudoInd", "boot.pseudoInd.sameX3",
                      "boot.pseudoNP", "boot.cond")){
    stop("For this test statistic, only the nonparametric bootstrap ",
         "is available.")
  }

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

    stop("Unknown bootstrap method name.")
  )

  switch(
    testStat,

    "T1_CvM_Cs3" = {

      FUN_trueStat <- testStat_T1_CvM_Cs3
      if (nStar == 1){
        FUN_stat_st <- testStat_T1_CvM_Cs3_boot1st
      }
    },

    "T1_CvM_Cs4" = {

      FUN_trueStat <- testStat_T1_CvM_Cs4
      if (nStar == 1){
        FUN_stat_st <- testStat_T1_CvM_Cs4_boot1st
      }
    },

    "tilde_T0_CvM" = {

      FUN_trueStat <- testStat_tilde_T0_CvM
      if (nStar == 1){
        FUN_stat_st <- testStat_tilde_T0_CvM_boot1st
      }
    },

    "T1_KS_Cs3" = {

      FUN_trueStat <- testStat_T1_KS_Cs3
      if (nStar == 1){
        FUN_stat_st <- testStat_T1_KS_Cs3_boot1st
      }
    },

    "T1_KS_Cs4" = {

      FUN_trueStat <- testStat_T1_KS_Cs4
      if (nStar == 1){
        FUN_stat_st <- testStat_T1_KS_Cs4_boot1st
      }
    },

    "tilde_T0_KS" = {

      FUN_trueStat <- testStat_tilde_T0_KS
      if (nStar == 1){
        FUN_stat_st <- testStat_tilde_T0_KS_boot1st
      }
    },

    "I_chi" = {

      FUN_trueStat <- testStat_Ichi
      if (nStar == 1){
        FUN_stat_st <- testStat_Ichi_boot1st
      } else if (nStar == 2){
        FUN_stat_st <- testStat_Ichi_boot2st
      }
    },

    "I_2n" = {

      FUN_trueStat <- testStat_I2n
      if (nStar == 1){
        FUN_stat_st <- testStat_I2n_boot1st
      } else if (nStar == 2){
        FUN_stat_st <- testStat_I2n_boot2st
      }
    },

    stop("Unknown test statistic name. Possible choices are: ",
         "'T1_CvM_Cs3', 'T1_CvM_Cs4', 'tilde_T0_CvM', ",
         "'T1_KS_Cs3', 'T1_KS_Cs4', 'tilde_T0_KS', ",
         "'I_chi', 'I_2n'. (see documentation)")
  )

  env <- environment()
  FUN_boot(env = env, FUN_trueStat = FUN_trueStat, FUN_stat_st = FUN_stat_st)

  return ( list(true_stat = env$true_stat ,
                vect_statB = env$vect_statB ,
                p_val = env$p_val) )
}


