

testStat_T2c <- function(env)
{
  # Computation of the pseudos-observations
  ecdf1 = stats::ecdf(env$X1)
  ecdf2 = stats::ecdf(env$X2)
  ecdf3 = stats::ecdf(env$X3)

  env$U1 = ecdf1(env$X1)
  env$U2 = ecdf2(env$X2)
  env$U3 = ecdf3(env$X3)

  # Computation of Z
  env$resultZ = estimationOfZ_I_J(U1 = env$U1, U2 = env$U2, U3 = env$U3,
                                  kernel = env$kernel.name, h = env$h)

  env$Z1 = env$resultZ$Z1
  env$Z2 = env$resultZ$Z2
  env$matrixK3 = env$resultZ$matrixK3

  # Estimation by conditional CMLE
  if (env$family == 2) {
    env$theta_0 = VineCopula::BiCopEst(u1 = env$Z1, u2 = env$Z2,
                                   family = 1 , method = "itau")$par

    env$theta_xJ = estimateParCondCopula_ZIJ(Z1_J = env$Z1, Z2_J = env$Z2,
                                             U3 = env$U3, newU3 = env$grid$nodes,
                                             family = 1 , method = "itau",
                                             h = env$h)
  } else {
    env$theta_0 = VineCopula::BiCopEst(u1 = env$Z1, u2 = env$Z2,
                                       family = env$family , method = "mle")$par

    env$theta_xJ = estimateParCondCopula_ZIJ(Z1_J = env$Z1, Z2_J = env$Z2,
                                             U3 = env$U3, newU3 = env$grid$nodes,
                                             family = env$family , method = "mle",
                                             h = env$h)
  }

  env$true_stat = sum( env$grid$weights * (env$theta_0 - env$theta_xJ)^2)
}


testStat_T2c_boot1st <- function(env)
{
  if (is.null(env$existZU_st)) {
    # Computation of the pseudos-observations
    ecdf1_st = stats::ecdf(env$X1_st)
    ecdf2_st = stats::ecdf(env$X2_st)
    ecdf3_st = stats::ecdf(env$X3_st)

    env$U1_st = ecdf1_st(env$X1_st)
    env$U2_st = ecdf2_st(env$X2_st)
    env$U3_st = ecdf3_st(env$X3_st)

    # Computation of Z
    env$resultZ_st = estimationOfZ_I_J(U1 = env$U1_st, U2 = env$U2_st, U3 = env$U3_st,
                                       kernel = env$kernel.name, h = env$h)

    env$Z1_st = env$resultZ_st$Z1
    env$Z2_st = env$resultZ_st$Z2
  }

  # Estimation by conditional CMLE
  if (env$family == 2) {
    env$theta_0_st = VineCopula::BiCopEst(u1 = env$Z1_st, u2 = env$Z2_st,
                                       family = 1 , method = "itau")$par

    env$theta_xJ_st = estimateParCondCopula_ZIJ(Z1_J = env$Z1_st, Z2_J = env$Z2_st,
                                             U3 = env$U3_st, newU3 = env$grid$nodes,
                                             family = 1 , method = "itau",
                                             h = env$h)
  } else {
    env$theta_0_st = VineCopula::BiCopEst(u1 = env$Z1_st, u2 = env$Z2_st,
                                       family = env$family , method = "mle")$par

    env$theta_xJ_st = estimateParCondCopula_ZIJ(Z1_J = env$Z1_st, Z2_J = env$Z2_st,
                                             U3 = env$U3_st, newU3 = env$grid$nodes,
                                             family = env$family , method = "mle",
                                             h = env$h)
  }

  env$stat_st = sum(env$grid$weights *
                      (env$theta_xJ_st - env$theta_xJ - env$theta_0_st + env$theta_0)^2)
}


testStat_T2c_boot2st <- function(env)
{
  if (is.null(env$existZU_st)) {
    # Computation of the pseudos-observations
    ecdf1_st = stats::ecdf(env$X1_st)
    ecdf2_st = stats::ecdf(env$X2_st)
    ecdf3_st = stats::ecdf(env$X3_st)

    env$U1_st = ecdf1_st(env$X1_st)
    env$U2_st = ecdf2_st(env$X2_st)
    env$U3_st = ecdf3_st(env$X3_st)

    # Computation of Z
    env$resultZ_st = estimationOfZ_I_J(U1 = env$U1_st, U2 = env$U2_st, U3 = env$U3_st,
                                       kernel = env$kernel.name, h = env$h)

    env$Z1_st = env$resultZ_st$Z1
    env$Z2_st = env$resultZ_st$Z2
  }

  # Estimation by conditional CMLE
  if (env$family == 2) {
    env$theta_0_st = VineCopula::BiCopEst(u1 = env$Z1_st, u2 = env$Z2_st,
                                          family = 1 , method = "itau")$par

    env$theta_xJ_st = estimateParCondCopula_ZIJ(Z1_J = env$Z1_st, Z2_J = env$Z2_st,
                                                U3 = env$U3_st, newU3 = env$grid$nodes,
                                                family = 1 , method = "itau",
                                                h = env$h)
  } else {
    env$theta_0_st = VineCopula::BiCopEst(u1 = env$Z1_st, u2 = env$Z2_st,
                                          family = env$family , method = "mle")$par

    env$theta_xJ_st = estimateParCondCopula_ZIJ(Z1_J = env$Z1_st, Z2_J = env$Z2_st,
                                                U3 = env$U3_st, newU3 = env$grid$nodes,
                                                family = env$family , method = "mle",
                                                h = env$h)
  }

  env$stat_st = sum(env$grid$weights *
                      (env$theta_xJ_st - env$theta_0_st)^2)
}
