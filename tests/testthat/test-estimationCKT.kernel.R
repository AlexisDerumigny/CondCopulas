test_that("CKT.kernel gives the same results with inputs that are matrices with 1 column", {

  set.seed(1)
  N = 50
  Z = rnorm(n = N)
  X1 = rnorm(n = N)
  X2 = rnorm(n = N)

  newZ = seq(-10, 10, by = 0.1)

  estimatedCKT_kernel <- CKT.kernel(
     X1 = X1, X2 = X2, Z = Z,
     newZ = newZ, h = 10, kernel.name = "Epa")

  estimatedCKT_kernel_matX1 <- CKT.kernel(
    X1 = matrix(X1, ncol = 1), X2 = X2, Z = Z,
    newZ = newZ, h = 10, kernel.name = "Epa")

  estimatedCKT_kernel_matX2 <- CKT.kernel(
    X1 = X1, X2 = matrix(X2, ncol = 1), Z = Z,
    newZ = newZ, h = 10, kernel.name = "Epa")

  estimatedCKT_kernel_matZ <- CKT.kernel(
    X1 = X1, X2 = X2, Z = matrix(Z, ncol = 1),
    newZ = newZ, h = 10, kernel.name = "Epa")

  estimatedCKT_kernel_dfZ <- CKT.kernel(
    X1 = X1, X2 = X2, Z = data.frame(Z = Z),
    newZ = newZ, h = 10, kernel.name = "Epa")

  expect_identical(estimatedCKT_kernel_matX1, estimatedCKT_kernel)
  expect_identical(estimatedCKT_kernel_matX1, estimatedCKT_kernel)
  expect_identical(estimatedCKT_kernel_matZ, estimatedCKT_kernel)
  expect_identical(estimatedCKT_kernel_dfZ, estimatedCKT_kernel)
})


test_that("old interface to CKT.kernel works", {

  set.seed(1)
  N = 50
  Z = rnorm(n = N)
  X1 = rnorm(n = N)
  X2 = rnorm(n = N)

  newZ = seq(-10, 10, by = 0.1)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 10, kernel.name = "Epa")

  estimatedCKT_kernel_old <- CKT.kernel(
    observedX1 = X1, observedX2 = X2, observedZ = Z,
    newZ = newZ, h = 10, kernel.name = "Epa")

  expect_identical(estimatedCKT_kernel_old, estimatedCKT_kernel)
})


test_that("CKT.kernel returns an error in case of wrong dimensions", {

  N = 10
  Z1 = rnorm(n = N, mean = 5, sd = 2)
  Z2 = rnorm(n = N, mean = 5, sd = 2)
  conditionalTau = -0.9 + 1.8 * pnorm(Z1 - Z2, mean = 2, sd = 2)
  simCopula = VineCopula::BiCopSim(N = N , family = 1,
      par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  newZ = expand.grid(Z1 = seq(2,10,by = 0.5),
                     Z2 = seq(2,10,by = 1))
  expect_error(
    {CKT.kernel(
      X1 = X1, X2 = X2, Z = Z1,
      newZ = newZ, h = 1, kernel.name = "Epa")},
    class = "WrongDimensionError"
  )

})


