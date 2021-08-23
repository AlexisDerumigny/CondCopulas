
#' Test of the assumption that a conditional copulas does not vary
#' through a list of discrete conditioning events
#'
#'
#' @param X1 vector of \code{n} observations of the first conditioned variable.
#' @param X2 vector of \code{n} observations of the second conditioned variable.
#' @param partition matrix of size \code{n * p},
#' where \code{p} is the number of conditioning events that are considered.
#' partition[i,j] should be the indicator of whether the \code{i}-th observation
#' belongs or not to the \code{j}-th conditioning event.
#'
#' @param family family of parametric copulas used
#' @param testStat test statistic used. Possible choices are
#' \itemize{
#'   \item \code{T2c_par} \eqn{\sum_{box} (\theta_0 - \theta(box))^2}
#'   \item \code{T2c_tau} Same as above, except that the copula family is now parametrized
#'   by its Kendall's tau instead of its natural parameter.
#' }
#' @param typeBoot
#' @param nBootstrap
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
#' n = 800
#' Z = stats::runif(n = n)
#' CKT = 0.2 * as.numeric(Z <= 0.3) +
#'   0.5 * as.numeric(Z > 0.3 & Z <= 0.5) +
#'   - 0.8 * as.numeric(Z > 0.5)
#' simCopula = VineCopula::BiCopSim(N = n,
#'   par = VineCopula::BiCopTau2Par(CKT, family = 1), family = 1)
#' X1 = simCopula[,1]
#' X2 = simCopula[,2]
#' partition = cbind(Z <= 0.3, Z > 0.3 & Z <= 0.5, Z > 0.5)
#'
#' result = bCond.simpA.param(X1 = X1, X2 = X2,
#'   partition = partition, family = 1, typeBoot = "boot.paramCond")
#' print(result$p_val)
#'
#' n = 800
#' Z = stats::runif(n = n)
#' CKT = 0.1
#' simCopula = VineCopula::BiCopSim(N = n,
#'   par = VineCopula::BiCopTau2Par(CKT, family = 1), family = 1)
#' X1 = simCopula[,1]
#' X2 = simCopula[,2]
#' partition = cbind(Z <= 0.3, Z > 0.3 & Z <= 0.5, Z > 0.5)
#'
#' result = bCond.simpA.param(X1 = X1, X2 = X2,
#'   partition = partition, family = 1, typeBoot = "boot.paramCond")
#' print(result$p_val)
#'
#'
#' @export
#'
bCond.simpA.param <- function(
  X1, X2, partition, family, testStat = "T2c_tau", typeBoot = "boot.NP",
  nBootstrap = 100)
{
  if (length(X1) != length(X2)){stop("X1 and X2 should be of the same length.")}
  if (length(X1) != nrow(partition)){
    stop("X1 should have the same length as the number of rows in 'partition'")
  }
  n <- length(X1)

  if (family == 2){
    method = "itau"; family_est <- 1
  } else {
    method = "mle"; family_est <- family
    }

  switch(
    typeBoot,

    "boot.NP" = {
      FUN_boot <- boot.NP.bCond
      nStar <- 1
    },

    "boot.paramInd" = {
      FUN_boot <- boot.paramInd.bCond
      nStar <- 2
    },

    "boot.paramCond" = {
      FUN_boot <- boot.paramCond.bCond
      nStar <- 1
    },

    stop("Unknown bootstrap method name.")
  )

  switch(
    testStat,

    "T2c_tau" = {

      parametrization = "tau"
      FUN_trueStat <- testStat_bT2c
      if (nStar == 1){
        FUN_stat_st <- testStat_bT2c_boot1st
      } else if (nStar == 2){
        FUN_stat_st <- testStat_bT2c_boot2st
      }
    },

    "T2c_par" = {

      parametrization = "par"
      FUN_trueStat <- testStat_bT2c
      if (nStar == 1){
        FUN_stat_st <- testStat_bT2c_boot1st
      } else if (nStar == 2){
        FUN_stat_st <- testStat_bT2c_boot2st
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



testStat_bT2c <- function(env){
  condPobs <- bCond.pobs(X1 = env$X1, X2 = env$X2, partition = env$partition)

  # Estimation of the simplified (conditional) parameter
  env$theta_0 = VineCopula::BiCopEst(u1 = condPobs[,1], u2 = condPobs[,2],
                                     family = env$family_est , method = env$method)$par

  env$theta_boxes = bCond.estParamCopula(U1 = condPobs[,1], U2 = condPobs[,2],
                                         family = env$family_est, partition = env$partition)

  if (parametrization == "par"){
    env$true_stat = sum((env$theta_0 - env$theta_boxes)^2)
  } else if (parametrization == "tau"){
    env$true_stat = sum(( VineCopula::BiCopPar2Tau(env$theta_0, family = family)
                          - VineCopula::BiCopPar2Tau(env$theta_boxes, family = family) )^2)
  } else {
    stop("Unknown parametrization. Possible parametrizations are 'tau' and 'par'.")
  }

}

testStat_bT2c_boot1st <- function(env){
  condPobs_st <- bCond.pobs(X1 = env$X1_st, X2 = env$X2_st, partition = env$partition_st)

  # Estimation of the simplified (conditional) parameter
  env$theta_0_st = VineCopula::BiCopEst(u1 = condPobs_st[,1], u2 = condPobs_st[,2],
                                    family = env$family_est , method = env$method)$par

  env$theta_boxes_st = bCond.estParamCopula(U1 = condPobs_st[,1], U2 = condPobs_st[,2],
                                        family = env$family_est, partition = env$partition_st)

  if (parametrization == "par"){
    env$stat_st = sum((env$theta_boxes_st - env$theta_boxes - env$theta_0_st + env$theta_0)^2)
  } else if (parametrization == "tau"){
    env$stat_st = sum(( VineCopula::BiCopPar2Tau(env$theta_boxes_st, family = family)
                        - VineCopula::BiCopPar2Tau(env$theta_boxes, family = family)
                        - VineCopula::BiCopPar2Tau(env$theta_0_st, family = family)
                        + VineCopula::BiCopPar2Tau(env$theta_0, family = family))^2)
  } else {
    stop("Unknown parametrization. Possible parametrizations are 'tau' and 'par'.")
  }

}

testStat_bT2c_boot2st <- function(env){
  condPobs_st <- bCond.pobs(X1 = env$X1_st, X2 = env$X2_st, partition = env$partition_st)

  # Estimation of the simplified (conditional) parameter
  env$theta_0_st = VineCopula::BiCopEst(u1 = condPobs_st[,1], u2 = condPobs_st[,2],
                                    family = env$family_est , method = env$method)$par

  env$theta_boxes_st = bCond.estParamCopula(U1 = condPobs_st[,1], U2 = condPobs_st[,2],
                                        family = env$family_est, partition = env$partition_st)

  if (parametrization == "par"){
    env$stat_st = sum( (env$theta_boxes_st - env$theta_0_st)^2 )
  } else if (parametrization == "tau"){
    env$stat_st = sum( ( VineCopula::BiCopPar2Tau(env$theta_boxes_st, family = family)
                         - VineCopula::BiCopPar2Tau(env$theta_0_st, family = family) )^2)
  } else {
    stop("Unknown parametrization. Possible parametrizations are 'tau' and 'par'.")
  }

}

