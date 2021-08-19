
# Usual non-parametric bootstrap ====================================================

#' Usual non-parametric bootstrap testing
#'
#' @param env the environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item nBootstrap: the number of bootstrap replications
#'    \item other variables necessary for the computation of the test statistic
#' }
#' @param FUN_trueStat the function for computing the (true) test statistic
#' @param FUN_stat_st the function for computing the bootstrapped test statistic.
#' This corresponds to the test statistic with *one star*.
#'
#' @noRd
#'
boot.NP <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    # Resampling to create the bootstrapped sample
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$X1_st = env$X1[permutation]
    env$X2_st = env$X2[permutation]
    env$X3_st = env$X3[permutation]

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


# Bootstrap schemes for tests of simpA based on the independence prop ================

#' Pseudo-independent bootstrap
#'
#'
#' @noRd
#'
boot.pseudoInd <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$existZU_st = TRUE

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    # Resampling to create the bootstrapped sample
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$Z1_st = env$Z1[permutation]
    env$Z2_st = env$Z2[permutation]

    permutation2 = as.integer(stats::runif(env$n, 1, env$n))
    env$U3_st = env$U3[permutation2]

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


#' Pseudo-independent bootstrap with same conditioning variables
#'
#'
#' @noRd
#'
boot.pseudoInd.sameX3 <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$existZU_st = TRUE

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    # Resampling to create the bootstrapped sample
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$Z1_st = env$Z1[permutation]
    env$Z2_st = env$Z2[permutation]

    # Keeping the same observations of the conditioning variable
    env$U3_st = env$U3

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


#' Pseudo-nonparametric bootstrap
#'
#'
#' @noRd
#'
boot.pseudoNP <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$existZU_st = TRUE

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    # Resampling to create the bootstrapped sample
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$Z1_st = env$Z1[permutation]
    env$Z2_st = env$Z2[permutation]
    env$U3_st = env$U3[permutation]

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


#' Conditional bootstrap
#'
#'
#' @noRd
#'
boot.cond <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$existZU_st = TRUE

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    # Resampling to create the bootstrapped sample
    permutation2 = as.integer(stats::runif(env$n, 1, env$n))
    env$U3_st = env$U3[permutation2]
    env$Z1_st = rep(NA, env$n)
    env$Z2_st = rep(NA, env$n)
    for (i in 1:env$n)
    {
      permutation_i = sample(x = 1:env$n, size = 1, prob = env$matrixK3[permutation2[i],])
      env$Z1_st[i] = env$Z1[permutation_i]
      env$Z2_st[i] = env$Z2[permutation_i]
    }

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


# Bootstrap schemes for parametric tests of simpA ================


#' Parametric independent bootstrap
#'
#'
#' @noRd
#'
boot.paramInd <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$existZU_J_st = TRUE

  env$vect_statB = rep(NA , env$nBootstrap)

  if (env$family == 2) {
    for (iBootstrap in 1:env$nBootstrap) {
      # Resampling to create the bootstrapped sample
      simCopule_st = VineCopula::BiCopSim(N = env$n , family = env$family,
                                          par = env$theta_0, par2 = 4)
      env$Z1_J_st = simCopule_st[,1]
      env$Z2_J_st = simCopule_st[,2]

      permutation2 = as.integer(stats::runif(env$n, 1, env$n))
      env$U3_st = env$U3[permutation2]

      FUN_stat_st(env)
      env$vect_statB[iBootstrap] <- env$stat_st
    }
  } else {
    for (iBootstrap in 1:env$nBootstrap) {
      # Resampling to create the bootstrapped sample
      simCopule_st = VineCopula::BiCopSim(N = env$n , family = env$family,
                                          par = env$theta_0)
      env$Z1_J_st = simCopule_st[,1]
      env$Z2_J_st = simCopule_st[,2]

      permutation2 = as.integer(stats::runif(env$n, 1, env$n))
      env$U3_st = env$U3[permutation2]

      FUN_stat_st(env)
      env$vect_statB[iBootstrap] <- env$stat_st
    }
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


#' Parametric conditional bootstrap
#'
#'
#' @noRd
#'
boot.paramCond <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$existZU_st = TRUE

  env$vect_statB = rep(NA , env$nBootstrap)

  if (env$family == 2) {

    theta_xJ_n = estimateParCondCopula_ZIJ(
        Z1_J = env$Z1_J, Z2_J = env$Z2_J, U3 = env$U3,
        newU3 = env$U3, family = 1, h = env$h, method = "itau")

    for (iBootstrap in 1:env$nBootstrap) {
      # Resampling to create the bootstrapped sample
      env$Z1_st = rep(NA, env$n)
      env$Z2_st = rep(NA, env$n)
      env$X3_st = rep(NA, env$n)
      for (i in 1:env$n)
      {
        i_X3 = as.integer(stats::runif(1,1,env$n))
        simCopule_st = VineCopula::BiCopSim(N=1 , family = env$family,
                                              par = env$theta_xJ_n[i_X3], par2 = 4)
        env$Z1_J_st[i] = simCopule_st[1]
        env$Z2_J_st[i] = simCopule_st[2]
        env$X3_st[i] = env$X3[i_X3]
      }
      ecdf3_st = stats::ecdf(env$X3_st)
      env$U3_st = ecdf3_st(env$X3_st)
    }
  } else if (env$family == 3) {
    theta_xJ_n = estimateParCondCopula_ZIJ(
      Z1_J = env$Z1_J, Z2_J = env$Z2_J, U3 = env$U3,
      newU3 = env$U3, family = env$family, h = env$h, method = "mle")

    for (iBootstrap in 1:env$nBootstrap) {
      # Resampling to create the bootstrapped sample
      env$Z1_J_st = rep(NA, env$n)
      env$Z2_J_st = rep(NA, env$n)
      env$X3_st = rep(NA, env$n)
      for (i in 1:env$n)
      {
        i_X3 = as.integer(stats::runif(1,1,env$n))
        simCopule_st = VineCopula::BiCopSim(N=1 , family = env$family,
                                            par = min(env$theta_xJ_n[i_X3],100))
        env$Z1_J_st[i] = simCopule_st[1]
        env$Z2_J_st[i] = simCopule_st[2]
        env$X3_st[i] = env$X3[i_X3]
      }
      ecdf3_st = stats::ecdf(env$X3_st)
      env$U3_st = ecdf3_st(env$X3_st)
    }
  } else {
    theta_xJ_n = estimateParCondCopula_ZIJ(
      Z1_J = env$Z1_J, Z2_J = env$Z2_J, U3 = env$U3,
      newU3 = env$U3, family = env$family, h = env$h, method = "mle")

    for (iBootstrap in 1:env$nBootstrap) {
      # Resampling to create the bootstrapped sample
      env$Z1_J_st = rep(NA, env$n)
      env$Z2_J_st = rep(NA, env$n)
      env$X3_st = rep(NA, env$n)
      for (i in 1:env$n)
      {
        i_X3 = as.integer(stats::runif(1,1,env$n))
        simCopule_st = VineCopula::BiCopSim(N=1 , family = env$family,
                                            par = env$theta_xJ_n[i_X3])
        env$Z1_J_st[i] = simCopule_st[1]
        env$Z2_J_st[i] = simCopule_st[2]
        env$X3_st[i] = env$X3[i_X3]
      }
      ecdf3_st = stats::ecdf(env$X3_st)
      env$U3_st = ecdf3_st(env$X3_st)
    }
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


