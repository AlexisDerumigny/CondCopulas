
test_that("measures_nonsimplifyingness_NP does not work for non-valid input", {

  X1 = c(1,2)
  X2 = 3
  Z = 4

  expect_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1) },
    class = "DifferentLengthsError")

  X1 = matrix(c(1,2), ncol = 2)

  expect_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1) },
    class = "WrongDimensionError")

  X1 = 2
  Z = matrix(c(1,2), ncol = 2)

  expect_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1) },
    class = "WrongDimensionError")

})
