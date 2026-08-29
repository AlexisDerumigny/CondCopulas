test_that("computationOf_G_chi_ZU returns box probabilities", {
  # one point in the first box: mass 1 in cell [1,1,1], 0 elsewhere
  expect_equal(computationOf_G_chi_ZU(0.1, 0.1, 0.1, 1)$G, array(c(1, rep(0, 124)), c(5, 5, 5)))

  # random data: must match direct box counts, and G_indep the product of marginals
  set.seed(123)
  Z <- matrix(runif(150), 50)
  res <- computationOf_G_chi_ZU(Z[, 1], Z[, 2], Z[, 3], 50)
  G <- unclass(unname(table(lapply(1:3, function(j) cut(Z[, j], 0:5 / 5, include.lowest = TRUE))))) / 50
  expect_equal(res$G, G)
  expect_equal(res$G_indep, outer(apply(G, 1:2, sum), apply(G, 3, sum)))
})
