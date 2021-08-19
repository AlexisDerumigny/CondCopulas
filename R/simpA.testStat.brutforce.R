
#' Computation of the marginal conditional pseudo-observations
#'
#' @noRd
#'
estimationOfZ_I_J <- function (U1, U2, U3, kernel, h)
{
  # Computation of the kernel
  matrixK3 = computeKernelMatrix(observedX = U3, newX = U3, kernel = kernel, h = h)

  # Computation of conditional cumulative distribution functions
  Z1 = estimateCondCDF_vec(observedX1 = U1, newX1 = U1, matrixK3 = matrixK3)

  Z2 = estimateCondCDF_vec(observedX1 = U2, newX1 = U2, matrixK3 = matrixK3)

  return (list(Z1 = Z1, Z2 = Z2, matrixK3 = matrixK3))
}


## Cramer von Mises test statistics ===================================================


#' Compute the test statistic T1_CvM
#' using the estimator Cs_3
#' (i.e. empirical averaging of the estimated conditional copula)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_CvM_Cs3 <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Estimation of the conditional copula
  env$array_C_IJ = estimateNPCondCopula(
    observedX1 = env$U1, observedX2 = env$U2, observedX3 = env$U3,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by averaging
  env$mat_C_sIJ = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ[i,j] = mean(env$array_C_IJ[i,j,])
    }
  }

  # Computation of the test statistic
  env$true_stat = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$true_stat = env$true_stat +
        env$grid$weights[i] * env$grid$weights[j] *
        mean(env$grid$weights * (env$array_C_IJ[i,j,] - env$mat_C_sIJ[i,j])^2)
    }
  }
}


#' Compute the bootstrapped version of the test statistic T1_CvM
#' using the estimator Cs_3
#' (i.e. empirical averaging of the estimated conditional copula)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the bootstrapped data vectors X1_st, X2_st, X3_st
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_CvM_Cs3_boot1st <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1_st = stats::ecdf(env$X1_st)
  ecdf2_st = stats::ecdf(env$X2_st)
  ecdf3_st = stats::ecdf(env$X3_st)

  env$U1_st = ecdf1_st(env$X1_st)
  env$U2_st = ecdf2_st(env$X2_st)
  env$U3_st = ecdf3_st(env$X3_st)

  # Estimation of the conditional copula
  env$array_C_IJ_st = estimateNPCondCopula(
    observedX1 = env$U1_st, observedX2 = env$U2_st, observedX3 = env$U3_st,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by averaging
  env$mat_C_sIJ_st = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ_st[i,j] = mean(env$array_C_IJ_st[i,j,])
    }
  }

  # Computation of the bootstrapped test statistic
  env$stat_st = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$stat_st = env$stat_st +
        env$grid$weights[i] * env$grid$weights[j] *
        mean(env$grid$weights * (env$array_C_IJ_st[i,j,] - env$array_C_IJ[i,j,]
                                 - env$mat_C_sIJ_st[i,j] + env$mat_C_sIJ[i,j])^2)
    }
  }
}


#' Compute the test statistic T1_CvM
#' using the estimator Cs_4
#' (i.e. empirical cdf of the marginal conditional pseudo-observations)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_CvM_Cs4 <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Estimation of the conditional copula
  env$array_C_IJ = estimateNPCondCopula(
    observedX1 = env$U1, observedX2 = env$U2, observedX3 = env$U3,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by the
  # empirical cdf of the marginal conditional pseudo-observations
  env$resultZ = estimationOfZ_I_J(U1 = env$U1, U2 = env$U2,
                                  U3 = env$U3, kernel = env$kernel.name, h = env$h)
  env$ecdf_cond1 = env$resultZ$Z1
  env$ecdf_cond2 = env$resultZ$Z2

  C_sIJ <- function(u1, u2) {
    indicator = (env$ecdf_cond1 <= u1) & (env$ecdf_cond2 <= u2)
    result = (length(which(indicator)))/env$n
    return (result)
  }
  env$mat_C_sIJ = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ[i,j] = C_sIJ(env$grid$nodes[i] , env$grid$nodes[j])
    }
  }

  # Computation of the test statistic
  env$true_stat = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$true_stat = env$true_stat +
        env$grid$weights[i] * env$grid$weights[j] *
        mean(env$grid$weights * (env$array_C_IJ[i,j,] - env$mat_C_sIJ[i,j])^2)
    }
  }
}



#' Compute the bootstrapped version of the test statistic T1_CvM
#' using the estimator Cs_4
#' (i.e. empirical cdf of the marginal conditional pseudo-observations)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the bootstrapped data vectors X1_st, X2_st, X3_st
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_CvM_Cs4_boot1st <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1_st = stats::ecdf(env$X1_st)
  ecdf2_st = stats::ecdf(env$X2_st)
  ecdf3_st = stats::ecdf(env$X3_st)

  env$U1_st = ecdf1_st(env$X1_st)
  env$U2_st = ecdf2_st(env$X2_st)
  env$U3_st = ecdf3_st(env$X3_st)

  # Estimation of the conditional copula
  env$array_C_IJ_st = estimateNPCondCopula(
    observedX1 = env$U1_st, observedX2 = env$U2_st, observedX3 = env$U3_st,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by the
  # empirical cdf of the marginal conditional pseudo-observations
  env$resultZ_st = estimationOfZ_I_J(U1 = env$U1_st, U2 = env$U2_st,
                                     U3 = env$U3_st, kernel = env$kernel.name, h = env$h)
  env$ecdf_cond1_st = env$resultZ_st$Z1
  env$ecdf_cond2_st = env$resultZ_st$Z2

  C_sIJ_st <- function(u1, u2) {
    indicator = (env$ecdf_cond1_st <= u1) & (env$ecdf_cond2_st <= u2)
    result = (length(which(indicator))) / env$n
    return (result)
  }
  env$mat_C_sIJ_st = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ_st[i,j] = C_sIJ_st(env$grid$nodes[i] , env$grid$nodes[j])
    }
  }


  # Computation of the test statistic
  env$stat_st = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$stat_st = env$stat_st +
        env$grid$weights[i] * env$grid$weights[j] *
        mean(env$grid$weights * (env$array_C_IJ_st[i,j,] - env$array_C_IJ[i,j,]
                                 - env$mat_C_sIJ_st[i,j] + env$mat_C_sIJ[i,j])^2)
    }
  }
}


#' Compute the test statistic tilde_T0
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_tilde_T0_CvM <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Estimation of the conditional copula
  env$array_C_IJ = estimateNPCondCopula(
    observedX1 = env$U1, observedX2 = env$U2, observedX3 = env$U3,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Computation of the test statistic
  env$true_stat = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$true_stat = env$true_stat +
        env$grid$weights[i] * env$grid$weights[j] *
        mean(outer(env$grid$weights, env$grid$weights) *
               (env$array_C_IJ[,,i] - env$array_C_IJ[,,j])^2)


    }
  }
}


#' Compute the bootstrapped version of the test statistic tilde_T0_CvM
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the bootstrapped data vectors X1_st, X2_st, X3_st
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_tilde_T0_CvM_boot1st <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1_st = stats::ecdf(env$X1_st)
  ecdf2_st = stats::ecdf(env$X2_st)
  ecdf3_st = stats::ecdf(env$X3_st)

  env$U1_st = ecdf1_st(env$X1_st)
  env$U2_st = ecdf2_st(env$X2_st)
  env$U3_st = ecdf3_st(env$X3_st)

  # Estimation of the conditional copula
  env$array_C_IJ_st = estimateNPCondCopula(
    observedX1 = env$U1_st, observedX2 = env$U2_st, observedX3 = env$U3_st,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Computation of the test statistic
  env$stat_st = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$stat_st = env$stat_st +
        env$grid$weights[i] * env$grid$weights[j] *
        mean(outer(env$grid$weights, env$grid$weights) *
               (env$array_C_IJ_st[,,i] - env$array_C_IJ[,,i] +
                  env$array_C_IJ[,,j] - env$array_C_IJ_st[,,j])^2)


    }
  }
}



# Kolmogorov-Smirnov test statistics ====================================


#' Compute the test statistic T1_KS
#' using the estimator Cs_3
#' (i.e. empirical averaging of the estimated conditional copula)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_KS_Cs3 <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Estimation of the conditional copula
  env$array_C_IJ = estimateNPCondCopula(
    observedX1 = env$U1, observedX2 = env$U2, observedX3 = env$U3,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by averaging
  env$mat_C_sIJ = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ[i,j] = mean(env$array_C_IJ[i,j,])
    }
  }

  # Computation of the test statistic
  env$true_stat = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$true_stat = max(env$true_stat ,
                          abs(env$array_C_IJ[i,j,] - env$mat_C_sIJ[i,j]) )
    }
  }
}


#' Compute the bootstrapped version of the test statistic T1_KS
#' using the estimator Cs_3
#' (i.e. empirical averaging of the estimated conditional copula)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the bootstrapped data vectors X1_st, X2_st, X3_st
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_KS_Cs3_boot1st <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1_st = stats::ecdf(env$X1_st)
  ecdf2_st = stats::ecdf(env$X2_st)
  ecdf3_st = stats::ecdf(env$X3_st)

  env$U1_st = ecdf1_st(env$X1_st)
  env$U2_st = ecdf2_st(env$X2_st)
  env$U3_st = ecdf3_st(env$X3_st)

  # Estimation of the conditional copula
  env$array_C_IJ_st = estimateNPCondCopula(
    observedX1 = env$U1_st, observedX2 = env$U2_st, observedX3 = env$U3_st,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by averaging
  env$mat_C_sIJ_st = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ_st[i,j] = mean(env$array_C_IJ_st[i,j,])
    }
  }

  # Computation of the bootstrapped test statistic
  env$stat_st = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$stat_st = max(env$stat_st,
                        abs(env$array_C_IJ_st[i,j,] - env$array_C_IJ[i,j,]
                                 - env$mat_C_sIJ_st[i,j] + env$mat_C_sIJ[i,j]) )
    }
  }
}


#' Compute the test statistic T1_KS
#' using the estimator Cs_4
#' (i.e. empirical cdf of the marginal conditional pseudo-observations)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_KS_Cs4 <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Estimation of the conditional copula
  env$array_C_IJ = estimateNPCondCopula(
    observedX1 = env$U1, observedX2 = env$U2, observedX3 = env$U3,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by the
  # empirical cdf of the marginal conditional pseudo-observations
  env$resultZ = estimationOfZ_I_J(U1 = env$U1, U2 = env$U2,
                                  U3 = env$U3, kernel = env$kernel.name, h = env$h)
  env$ecdf_cond1 = env$resultZ$Z1
  env$ecdf_cond2 = env$resultZ$Z2

  C_sIJ <- function(u1, u2) {
    indicator = (env$ecdf_cond1 <= u1) & (env$ecdf_cond2 <= u2)
    result = (length(which(indicator)))/env$n
    return (result)
  }
  env$mat_C_sIJ = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ[i,j] = C_sIJ(env$grid$nodes[i] , env$grid$nodes[j])
    }
  }

  # Computation of the test statistic
  env$true_stat = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$true_stat = max(env$true_stat ,
                          abs(env$array_C_IJ[i,j,] - env$mat_C_sIJ[i,j]) )
    }
  }
}



#' Compute the bootstrapped version of the test statistic T1_KS
#' using the estimator Cs_4
#' (i.e. empirical cdf of the marginal conditional pseudo-observations)
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the bootstrapped data vectors X1_st, X2_st, X3_st
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_T1_KS_Cs4_boot1st <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1_st = stats::ecdf(env$X1_st)
  ecdf2_st = stats::ecdf(env$X2_st)
  ecdf3_st = stats::ecdf(env$X3_st)

  env$U1_st = ecdf1_st(env$X1_st)
  env$U2_st = ecdf2_st(env$X2_st)
  env$U3_st = ecdf3_st(env$X3_st)

  # Estimation of the conditional copula
  env$array_C_IJ_st = estimateNPCondCopula(
    observedX1 = env$U1_st, observedX2 = env$U2_st, observedX3 = env$U3_st,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Estimation of the simplified copula by the
  # empirical cdf of the marginal conditional pseudo-observations
  env$resultZ_st = estimationOfZ_I_J(U1 = env$U1_st, U2 = env$U2_st,
                                     U3 = env$U3_st, kernel = env$kernel.name, h = env$h)
  env$ecdf_cond1_st = env$resultZ_st$Z1
  env$ecdf_cond2_st = env$resultZ_st$Z2

  C_sIJ_st <- function(u1, u2) {
    indicator = (env$ecdf_cond1_st <= u1) & (env$ecdf_cond2_st <= u2)
    result = (length(which(indicator))) / env$n
    return (result)
  }
  env$mat_C_sIJ_st = matrix(nrow = env$nGrid , ncol = env$nGrid)
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$mat_C_sIJ_st[i,j] = C_sIJ_st(env$grid$nodes[i] , env$grid$nodes[j])
    }
  }


  # Computation of the test statistic
  env$stat_st = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$stat_st = max(env$stat_st,
                        abs(env$array_C_IJ_st[i,j,] - env$array_C_IJ[i,j,]
                                 - env$mat_C_sIJ_st[i,j] + env$mat_C_sIJ[i,j]) )
    }
  }
}


#' Compute the test statistic tilde_T0_KS
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the data vectors X1, X2, X3
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_tilde_T0_KS <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Estimation of the conditional copula
  env$array_C_IJ = estimateNPCondCopula(
    observedX1 = env$U1, observedX2 = env$U2, observedX3 = env$U3,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Computation of the test statistic
  env$true_stat = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$true_stat = max(env$true_stat ,
                          abs( env$array_C_IJ[,,i] - env$array_C_IJ[,,j]) )


    }
  }
}


#' Compute the bootstrapped version of the test statistic tilde_T0_KS
#'
#' @param env an environment which should contain (at least)
#' \itemize{
#'    \item the bootstrapped data vectors X1_st, X2_st, X3_st
#'    \item \code{grid} list of nodes and weights
#'    at which the conditional copula should be estimated.
#'    \item \code{nGrid} the number of nodes
#' }
#'
#' @noRd
testStat_tilde_T0_KS_boot1st <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1_st = stats::ecdf(env$X1_st)
  ecdf2_st = stats::ecdf(env$X2_st)
  ecdf3_st = stats::ecdf(env$X3_st)

  env$U1_st = ecdf1_st(env$X1_st)
  env$U2_st = ecdf2_st(env$X2_st)
  env$U3_st = ecdf3_st(env$X3_st)

  # Estimation of the conditional copula
  env$array_C_IJ_st = estimateNPCondCopula(
    observedX1 = env$U1_st, observedX2 = env$U2_st, observedX3 = env$U3_st,
    U1_ = env$grid$nodes, U2_ = env$grid$nodes, newX3 = env$grid$nodes,
    kernel = env$kernel.name, h = env$h)

  # Computation of the test statistic
  env$stat_st = 0
  for (i in 1:env$nGrid)
  {
    for (j in 1:env$nGrid)
    {
      env$stat_st = max( env$stat_st ,
                         abs(env$array_C_IJ_st[,,i] - env$array_C_IJ[,,i] +
                               env$array_C_IJ[,,j] - env$array_C_IJ_st[,,j]) )
    }
  }
}



