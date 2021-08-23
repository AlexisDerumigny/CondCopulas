


#' Nonparametric testing of the simplifying assumption
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
#' result <- simpA.NP(
#'    X1 = X1, X2 = X2, X3 = Z,
#'    testStat = "I_chi", typeBoot = "boot.pseudoInd",
#'    h = 0.03, kernel.name = "Epanechnikov")
#' print(result$p_val)
#'
#' set.seed(1)
#' N = 800
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
#'    h = 0.08, kernel.name = "Epanechnikov")
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
         "'I_chi', 'I_2n'.")
  )

  env <- environment()
  FUN_boot(env = env, FUN_trueStat = FUN_trueStat, FUN_stat_st = FUN_stat_st)

  return ( list(true_stat = env$true_stat ,
                vect_statB = env$vect_statB ,
                p_val = env$p_val) )
}


