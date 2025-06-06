

test_that("the truncation works", {
  set.seed(1)
  N = 500
  Z = rnorm(n = N, mean = 5, sd = 2)
  conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
  simCopula = VineCopula::BiCopSim(N=N , family = 1,
                                   par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
  X1 = qnorm(simCopula[,1], mean = Z)
  X2 = qnorm(simCopula[,2], mean = - Z)

  result_1 <- simpA.param(
    X1 = X1, X2 = X2, X3 = Z, family = 1,
    h = 2, kernel.name = "Epanechnikov", nBootstrap = 1, truncVal = NULL)

  result_2 <- simpA.param(
    X1 = X1, X2 = X2, X3 = Z, family = 1,
    h = 2, kernel.name = "Epanechnikov", nBootstrap = 1, truncVal = 0)

  expect_identical(result_1$true_stat, result_2$true_stat)

  expect_error(
    {simpA.param(
      X1 = X1, X2 = X2, X3 = Z, family = 1,
      h = 2, kernel.name = "Epanechnikov", nBootstrap = 1, truncVal = -1)},

    class = "InvalidInputError"
  )

  expect_error(
    {simpA.param(
      X1 = X1, X2 = X2, X3 = Z, family = 1,
      h = 2, kernel.name = "Epanechnikov", nBootstrap = 1, truncVal = 0.5)},

    class = "InvalidInputError"
  )

})
