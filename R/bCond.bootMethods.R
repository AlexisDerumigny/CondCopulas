

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
boot.NP.bCond <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    # Resampling to create the bootstrapped sample
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$X1_st = env$X1[permutation]
    env$X2_st = env$X2[permutation]
    env$partition_st = env$partition[permutation, ]

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


boot.paramInd.bCond <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    if (env$family == 2) {
      simCopule_st = VineCopula::BiCopSim(N = env$n , family = env$family,
                                          par = env$theta_0, par2 = 4)
    }
    else {
      simCopule_st = VineCopula::BiCopSim(N = env$n , family = env$family, par = env$theta_0)
    }
    env$X1_st = simCopule_st[,1]
    env$X2_st = simCopule_st[,2]

    # Independent resampling of the conditional events
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$partition_st = env$partition[permutation, ]

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}


boot.paramCond.bCond <- function(env, FUN_trueStat, FUN_stat_st)
{
  FUN_trueStat(env)

  env$vect_statB = rep(NA , env$nBootstrap)
  for (iBootstrap in 1:env$nBootstrap)
  {
    permutation = as.integer(stats::runif(env$n, 1, env$n))
    env$partition_st = env$partition[permutation, ]

    simCopule_st = matrix(nrow = env$n, ncol = 2)
    for (box in 1:ncol(env$partition))
    {
      nobs_box = length(which(env$partition_st[,box]))
      where_box = as.logical(env$partition_st[,box])
      if (env$family == 2) {
        simCopule_st[where_box,] = VineCopula::BiCopSim(
          N = nobs_box, family = env$family, par = env$theta_boxes[box], par2 = 4)
      } else {
        simCopule_st[where_box,] = VineCopula::BiCopSim(
          N = nobs_box, family = env$family, par = env$theta_boxes[box])
      }
    }

    env$X1_st = simCopule_st[,1]
    env$X2_st = simCopule_st[,2]

    FUN_stat_st(env)
    env$vect_statB[iBootstrap] <- env$stat_st
  }

  env$p_val = NA
  try(env$p_val <-  1-stats::ecdf(env$vect_statB)(env$true_stat) , silent = TRUE)
}

