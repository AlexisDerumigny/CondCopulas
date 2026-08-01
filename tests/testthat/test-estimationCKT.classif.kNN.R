test_that("multiplication works", {
  # We simulate from a conditional copula
  set.seed(1)
  N = 800
  Z = rnorm(n = N, mean = 5, sd = 2)
  conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
  simCopula = VineCopula::BiCopSim(N=N , family = 1,
      par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
  X1 = qnorm(simCopula[,1])
  X2 = qnorm(simCopula[,2])

  newZ = seq(2,10,by = 0.1)
  datasetP = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
  expect_no_error({
    estimatedCKT_knn <- CKT.predict.kNN(
      datasetPairs = datasetP,
      newZ = matrix(newZ,ncol = 1),
      number_nn = c(50,80),
      partition = NULL)
  })
})
