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


test_that("CKT.kernel works without error for the different types of estimators in dim 1", {

  set.seed(1)
  N = 50
  Z = rnorm(n = N)
  X1 = rnorm(n = N)
  X2 = rnorm(n = N)

  newZ = seq(-10, 10, by = 1)

  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = "wdm")
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 1)
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 2)
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 3)
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 4)
  })

  # Expect no warning nor errors with Gaussian kernel

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = "wdm")

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 1)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 2)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 3)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 4)
})



test_that("CKT.kernel works without error for the different types of estimators in dim 2", {

  set.seed(1)
  N = 50
  Z = cbind(rnorm(n = N), rnorm(n = N))
  X1 = rnorm(n = N)
  X2 = rnorm(n = N)

  newZ = expand.grid(Z1 = seq(-3, 3, by = 0.5),
                     Z2 = seq(-3, 3, by = 1))

  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = "wdm")
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 1)
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 2)
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 3)
  })
  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = 4)
  })

  # Expect no warning nor errors with Gaussian kernel

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = "wdm")

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 1)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 2)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 3)

  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Gaussian", typeEstCKT = 4)
})



test_that("CKT.kernel gives warning in presence of NAs", {

  set.seed(1)
  N = 50
  Z = cbind(rnorm(n = N), rnorm(n = N))
  X1 = rnorm(n = N)
  X2 = rnorm(n = N)

  newZ = expand.grid(Z1 = seq(-3, 3, by = 0.5),
                     Z2 = seq(-3, 3, by = 1))

  newZ[1,1] = NA
  X1[2] = NA

  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = "wdm")
  })

  # Warnings can be disabled with the corresponding option
  estimatedCKT_kernel <- CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = newZ, h = 1, kernel.name = "Epa", typeEstCKT = "wdm", warnNA = FALSE)
})



test_that("Cross-validation works for CKT.kernel", {

  set.seed(1)
  N = 50
  Z = rnorm(n = N)
  X1 = rnorm(n = N)
  X2 = rnorm(n = N)

  newZ = seq(-2, 2, by = 0.5)


  expect_no_error({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = c(0.1, 1, 10), kernel.name = "Epa", methodCV = "Kfolds")
  })

  expect_no_error({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = c(0.01, 0.1, 1, 10), kernel.name = "Epa", methodCV = "leave-one-out")
  })

  expect_warning({
    estimatedCKT_kernel <- CKT.kernel(
      X1 = X1, X2 = X2, Z = Z,
      newZ = newZ, h = c(0.01, 0.1, 1, 10), kernel.name = "Epa", methodCV = "leave-one-out",
      nPairs = Inf)
  })
})

