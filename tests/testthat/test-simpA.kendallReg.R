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

  result$coef
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

  print(result$p_val)
})
