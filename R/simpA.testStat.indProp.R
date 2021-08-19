

# Utils ===================================================================


computationOf_G_chi_ZU <- function (Z1, Z2, U3, n)
{

  # Computation of GI,J(x1 , x2, x3)

  fonctionG_IJ <- function (x1, x2, x3)
  {
    return (mean(as.numeric((Z1 <= x1) & (Z2 <= x2) & (U3 <= x3) )))
  }

  # Computation of GI,J(Bk, Ai)

  listG = array(dim = c(5 , 5 , 5))
  partition = seq(from = 0 , to = 1 , length.out = 6)


  for (k1 in (1:5))
  {
    for (k2 in 1:5)
    {
      for (k3 in 1:5)
      {
        listG[k1, k2, k3] =
          fonctionG_IJ(partition[k1+1] , partition[k2+1] , partition[k3+1])

        - fonctionG_IJ(partition[k1  ] , partition[k2+1] , partition[k3+1])
        - fonctionG_IJ(partition[k1+1] , partition[k2  ] , partition[k3+1])
        - fonctionG_IJ(partition[k1+1] , partition[k2+1] , partition[k3  ])

        + fonctionG_IJ(partition[k1+1] , partition[k2  ] , partition[k3  ])
        + fonctionG_IJ(partition[k1  ] , partition[k2+1] , partition[k3  ])
        + fonctionG_IJ(partition[k1  ] , partition[k2  ] , partition[k3+1])

        - fonctionG_IJ(partition[k1] , partition[k2] , partition[k3])
      }
    }
  }

  # Computation of GI,J(Bk, R)

  listG_BkR = matrix(nrow = 5, ncol = 5)

  for (k1 in 1:5)
  {
    for (k2 in 1:5)
    {
      listG_BkR[k1 , k2] =
        fonctionG_IJ(partition[k1+1] , partition[k2+1] , 1)

      - fonctionG_IJ(partition[k1  ] , partition[k2+1] , 1)
      - fonctionG_IJ(partition[k1+1] , partition[k2  ] , 1)
      - fonctionG_IJ(partition[k1+1] , partition[k2+1] , 0)

      + fonctionG_IJ(partition[k1+1] , partition[k2  ] , 0)
      + fonctionG_IJ(partition[k1  ] , partition[k2+1] , 0)
      + fonctionG_IJ(partition[k1  ] , partition[k2  ] , 1)

      - fonctionG_IJ(partition[k1] , partition[k2] , 0)
    }
  }

  # Computation of GI,J(R^2, Al)

  listG_R2Al = rep(NA , n)

  for (k3 in 1:5)
  {

    listG_R2Al[k3] =
      fonctionG_IJ(1               , 1               , partition[k3+1])

    - fonctionG_IJ(0               , 1               , partition[k3+1])
    - fonctionG_IJ(1               , 0               , partition[k3+1])
    - fonctionG_IJ(1               , 1               , partition[k3  ])

    + fonctionG_IJ(1               , 0               , partition[k3  ])
    + fonctionG_IJ(0               , 1               , partition[k3  ])
    + fonctionG_IJ(0               , 0               , partition[k3+1])

    - fonctionG_IJ(0               , 0               , partition[k3])

  }

  # Computation of G indep

  listG_indep = array(dim = c(5 , 5 , 5))

  for (k1 in (1:5))
  {
    for (k2 in 1:5)
    {
      for (k3 in 1:5)
      {
        listG_indep[k1, k2, k3] = listG_BkR[k1 , k2] * listG_R2Al[k3]
      }
    }
  }

  return (list(G=listG , G_indep=listG_indep))
}


# Test statistics =========================================================


testStat_Ichi <- function(env)
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

  # Computation of G
  env$resultG = computationOf_G_chi_ZU(Z1 = env$Z1, Z2 = env$Z2,
                                       U3 = env$U3, n = env$n)
  env$listG = env$resultG$G
  env$listG_indep = env$resultG$G_indep
  env$true_stat = sum((env$listG - env$listG_indep)^2 / env$listG_indep)
}


testStat_Ichi_boot1st <- function(env)
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

  # Computation of G
  env$resultG_st = computationOf_G_chi_ZU(Z1 = env$Z1_st, Z2 = env$Z2_st,
                                          U3 = env$U3_st, n = env$n)
  env$listG_st = env$resultG_st$G
  env$listG_indep_st = env$resultG_st$G_indep
  env$stat_st = sum((env$listG_st - env$listG
                     - env$listG_indep_st + env$listG_indep)^2 / env$listG_indep_st)
}


testStat_Ichi_boot2st <- function(env)
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

  # Computation of G
  env$resultG_st = computationOf_G_chi_ZU(Z1 = env$Z1_st, Z2 = env$Z2_st,
                                          U3 = env$U3_st, n = env$n)
  env$listG_st = env$resultG_st$G
  env$listG_indep_st = env$resultG_st$G_indep
  env$stat_st = sum((env$listG_st - env$listG_indep_st)^2 / env$listG_indep_st)
}



testStat_I2n <- function(env)
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

  # Computation of the test statistic
  env$array_C_IJ = array(dim = c(env$nGrid, env$nGrid, env$nGrid))
  for (k in 1:env$nGrid)
  {
    for (i in 1:env$nGrid)
    {
      for (j in 1:env$nGrid)
      {
        env$C_IJ_ijk = mean(as.numeric(env$Z1_J <= env$grid$nodes[i] &
                                         env$Z2_J <= env$grid$nodes[j] &
                                         env$U3 <= env$grid$nodes[k]))
        env$C_I_given_J_ij = mean(as.numeric(env$Z1_J <= env$grid$nodes[i] &
                                               env$Z2_J <= env$grid$nodes[j]))
        env$C_J_k = mean(as.numeric(env$U3 <= env$grid$nodes[k]))

        env$array_C_IJ[i,j,k] = env$grid$weights[i] *
          env$grid$weights[j] * env$grid$weights[k] *
          (env$C_IJ_ijk - env$C_I_given_J_ij * env$C_J_k)^2
      }
    }
  }

  env$true_stat = sum(env$array_C_IJ)
}


testStat_I2n_boot1st <- function(env)
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

  # Computation of G
  env$array_C_IJ_st = array(dim = c(env$nGrid, env$nGrid, env$nGrid))
  for (k in 1:env$nGrid)
  {
    for (i in 1:env$nGrid)
    {
      for (j in 1:env$nGrid)
      {
        env$C_IJ_ijk_st = mean(as.numeric(env$Z1_J_st <= env$grid$nodes[i]
                                      & env$Z2_J_st <= env$grid$nodes[j]
                                      & env$U3_st <= env$grid$nodes[k]))

        env$C_I_given_J_ij_st = mean(as.numeric(env$Z1_J_st <= env$grid$nodes[i]
                                            & env$Z2_J_st <= env$grid$nodes[j]))

        env$C_J_k_st = mean(as.numeric(env$U3_st <= env$grid$nodes[k]))

        env$array_C_IJ_st[i,j,k] = env$grid$weights[i] *
          env$grid$weights[j] * env$grid$weights[k] *
          (env$C_IJ_ijk_st
           - env$C_IJ_ijk
           - env$C_I_given_J_ij_st * env$C_J_k_st
           + env$C_I_given_J_ij * env$C_J_k)^2
      }
    }
  }

  env$stat_st = sum(env$array_C_IJ_st)
}


testStat_I2n_boot2st <- function(env)
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

  # Computation of G
  env$array_C_IJ_st = array(dim = c(env$nGrid, env$nGrid, env$nGrid))
  for (k in 1:env$nGrid)
  {
    for (i in 1:env$nGrid)
    {
      for (j in 1:env$nGrid)
      {
        env$C_IJ_ijk_st = mean(as.numeric(env$Z1_J_st <= env$grid$nodes[i]
                                          & env$Z2_J_st <= env$grid$nodes[j]
                                          & env$U3_st <= env$grid$nodes[k]))

        env$C_I_given_J_ij_st = mean(as.numeric(env$Z1_J_st <= env$grid$nodes[i]
                                                & env$Z2_J_st <= env$grid$nodes[j]))

        env$C_J_k_st = mean(as.numeric(env$U3_st <= env$grid$nodes[k]))

        env$array_C_IJ_st[i,j,k] = env$grid$weights[i] *
          env$grid$weights[j] * env$grid$weights[k] *
          (env$C_IJ_ijk_st
           - env$C_I_given_J_ij_st * env$C_J_k_st)^2
      }
    }
  }

  env$stat_st = sum(env$array_C_IJ_st)
}


