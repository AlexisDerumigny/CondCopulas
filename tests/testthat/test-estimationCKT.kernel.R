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

  expect_identical(estimatedCKT_kernel_matX1, estimatedCKT_kernel)
  expect_identical(estimatedCKT_kernel_matX1, estimatedCKT_kernel)
  expect_identical(estimatedCKT_kernel_matZ, estimatedCKT_kernel)
})
