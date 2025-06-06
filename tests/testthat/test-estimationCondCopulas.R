

test_that("estimateNPCondCopula gives suitable warning for h too small", {

  # We simulate from a conditional copula
  N = 500
  X3 = rnorm(n = N, mean = 5, sd = 2)
  conditionalTau = 0.9 * pnorm(X3, mean = 5, sd = 2)
  simCopula = VineCopula::BiCopSim(N=N , family = 3,
                                   par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  # We do the estimation
  grid = c(0.2, 0.4, 0.6, 0.8)

  # Warning when 2 out of 4 points are too far away
  capturedWarning = tryCatch({
    estimateNPCondCopula(
      X1 = X1, X2 = X2, X3 = X3,
      U1_ = grid, U2_ = grid, newX3 = c(2, 5, 7, 100, 200),
      kernel = "Gaussian", h = 0.1)},
    warning = function(w){w}
  )

  expect_true(inherits(capturedWarning, "ZeroWeights_KernelWarning"))
  expect_identical(capturedWarning$problematicPoints, c(100, 200))

  # Warning when 2 out of 2 points are too far away
  capturedWarning = tryCatch({
    estimateNPCondCopula(
      X1 = X1, X2 = X2, X3 = X3,
      U1_ = grid, U2_ = grid, newX3 = c(2, 5, 7, 100, 200),
      kernel = "Gaussian", h = 0.1)},
    warning = function(w){w}
  )

  expect_true(inherits(capturedWarning, "ZeroWeights_KernelWarning"))
  expect_identical(capturedWarning$problematicPoints, c(100, 200))


  # Test whether the presence of the warning (because of the problematic points)
  # does not influence the estimation for the non-problematic points

  arrayEst_noproblem = estimateNPCondCopula(
    X1 = X1, X2 = X2, X3 = X3,
    U1_ = grid, U2_ = grid, newX3 = c(2, 5, 7),
    kernel = "Gaussian", h = 0.1)

  suppressWarnings(
    {arrayEst <- estimateNPCondCopula(
      X1 = X1, X2 = X2, X3 = X3,
      U1_ = grid, U2_ = grid, newX3 = c(2, 5, 7, 1000),
      kernel = "Gaussian", h = 0.1)},

    classes = "ZeroWeights_KernelWarning"
  )

  expect_true(all(is.na(arrayEst[, , 4])))

  expect_identical(arrayEst[, , -4], arrayEst_noproblem)
})



