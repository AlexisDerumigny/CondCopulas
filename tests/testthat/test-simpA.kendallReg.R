
test_that("simpA.kendallReg does not work for multivariate X1 or X2", {

  X1 = matrix(data = c(1,1,1,1), ncol = 2)
  X2 = c(1,1)
  Z = c(1,2)

  expect_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 1,
      listPhi = list( function(x){return(x)} ),
      lambda = 0)
  }, class = "WrongDimensionError")

  X1 = c(1,1)
  X2 = matrix(data = c(1,1,1,1), ncol = 2)
  Z = c(1,2)

  expect_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 1,
      listPhi = list( function(x){return(x)} ),
      lambda = 0)
  }, class = "WrongDimensionError")
})

test_that("simpA.kendallReg does not work for X1, X2 of different lengths", {

  X1 = c(1,1,3)
  X2 = c(1,1)
  Z = c(1,2)

  expect_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 1,
      listPhi = list( function(x){return(x)} ),
      lambda = 0)
  }, class = "DifferentLengthsError")

  X1 = c(1,1)
  X2 = c(1,1,3)
  Z = c(1,2)

  expect_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 1,
      listPhi = list( function(x){return(x)} ),
      lambda = 0)
  }, class = "DifferentLengthsError")

  X1 = c(1,1)
  X2 = c(1,1)
  Z = c(1,2,3)

  expect_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 1,
      listPhi = list( function(x){return(x)} ),
      lambda = 0)
  }, class = "DifferentLengthsError")
})

test_that("simpA.kendallReg works if only one function 'phi' is given", {
  set.seed(1)
  # We simulate from a conditional copula  set.seed(1)
  N = 300
  Z = runif(n = N, min = 0, max = 1)
  conditionalTau = -0.9 + 1.8 * Z

  simCopula = VineCopula::BiCopSim(N=N , family = 1,
      par = VineCopula::BiCopTau2Par(1 , conditionalTau ))

  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  expect_no_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 0.03,
      listPhi = list( function(x){return(x)} ) )
  })

  expect_equal(result$p_val, 0, tolerance = 0.02)

  coefs = coef(result)
  expect_true(!is.null(coefs))

  varcov = vcov(result)
  expect_true(!is.null(varcov))
})

test_that("simpA.kendallReg works if SA is true", {
  set.seed(1)
  # We simulate from a conditional copula  set.seed(1)
  N = 300
  Z = runif(n = N, min = 0, max = 1)
  conditionalTau = 0.5

  simCopula = VineCopula::BiCopSim(N=N , family = 1,
                                   par = VineCopula::BiCopTau2Par(1 , conditionalTau ))

  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  expect_no_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 0.03,
      listPhi = list( function(x){return(x)} ) )
  })

  expect_gt(result$p_val, 0.05)


  coefs = coef(result)
  expect_true(!is.null(coefs))

  varcov = vcov(result)
  expect_true(!is.null(varcov))
})

test_that("simpA.kendallReg works if only two functions 'phi' are given", {
  set.seed(1)
  # We simulate from a conditional copula  set.seed(1)
  N = 300
  Z = runif(n = N, min = 0, max = 1)
  conditionalTau = -0.9 + 1.8 * Z

  simCopula = VineCopula::BiCopSim(N=N , family = 1,
                                   par = VineCopula::BiCopTau2Par(1 , conditionalTau ))

  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  expect_no_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 0.03,
      listPhi = list( function(x){return(x)} ,
                      function(x){return(x^2)} ),
      lambda = 0)
  })

  expect_equal(result$p_val, 0, tolerance = 0.02)
  result$coef
  result$varCov

  print(result)
  plot(result)

  result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 0.03,
      listPhi = list( function(x){return(x)} ,
                      function(x){return(x^2)} ),
      lambda = 0.1)
  # lambda too high to get reliable p-values...
})

test_that("simpA.kendallReg works with many phi", {
  set.seed(1)
  # We simulate from a conditional copula  set.seed(1)
  N = 300
  Z = runif(n = N, min = 0, max = 1)
  conditionalTau = -0.9 + 1.8 * Z

  simCopula = VineCopula::BiCopSim(N = N , family = 1,
                                   par = VineCopula::BiCopTau2Par(1 , conditionalTau ))

  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  result = simpA.kendallReg(
     X1, X2, Z, h_kernel = 0.03,
     listPhi = list(
       function(x){return(x)},
       function(x){return(cos(10 * x))},
       function(x){return(sin(10 * x))},
       function(x){return(as.numeric(x <= 0.4))},
       function(x){return(as.numeric(x <= 0.6))}) )

  result = simpA.kendallReg(
    X1, X2, Z, h_kernel = 0.03,
    listPhi = list(
      function(x){return(x)}) )

  expect_equal(result$p_val, 0, tolerance = 0.02)

  # We simulate from a simplified conditional copula
  set.seed(1)
  N = 300
  Z = runif(n = N, min = 0, max = 1)
  conditionalTau = -0.3
  simCopula = VineCopula::BiCopSim(N=N , family = 1,
      par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  result = simpA.kendallReg(
     X1, X2, Z, h_kernel = 0.03,
     listPhi = list(
       function(x){return(x)},
       function(x){return(cos(10 * x))},
       function(x){return(sin(10 * x))},
       function(x){return(as.numeric(x <= 0.4))},
       function(x){return(as.numeric(x <= 0.6))}) )
  print(result)
  plot(result)

  expect_gt(result$p_val, 0.05)

  coefs = coef(result)
  expect_true(!is.null(coefs))

  varcov = vcov(result)
  expect_true(!is.null(varcov))
})

test_that("simpA.kendallReg works if only two functions 'phi' are given and SA is true", {
  set.seed(1)
  # We simulate from a conditional copula  set.seed(1)
  N = 300
  Z = runif(n = N, min = 0, max = 1)
  conditionalTau = 0.5

  simCopula = VineCopula::BiCopSim(N=N , family = 1,
                                   par = VineCopula::BiCopTau2Par(1 , conditionalTau ))

  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  expect_no_error({
    result = simpA.kendallReg(
      X1, X2, Z, h_kernel = 0.03,
      listPhi = list( function(x){return(x)} ,
                      function(x){return(x^2)} ),
      lambda = 0)
  })

  expect_gt(result$p_val, 0.05)

  # plot(result_$vectorZToEstimate, result_$vector_hat_CKT_NP^2, col = "blue", type = "l")
  # lines(result_$vectorZToEstimate, result_$resultWn$Gn_zipr, type = "l")
})
