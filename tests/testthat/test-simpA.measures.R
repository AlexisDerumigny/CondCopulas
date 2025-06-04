
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

  Z = 2

  # We now give invalid measures

  expect_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1, measures = c()) },
    class = "ZeroLengthError")

  expect_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1, measures = "aaaa") },
    class = "UnknownMeasureNameError")

})


test_that("measures_nonsimplifyingness_NP works for valid input", {

  X1 = rnorm(100)
  X2 = rnorm(100)
  Z = rnorm(100)

  # It works by default
  expect_no_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 0.1) })

  # It works for one of the measures
  expect_no_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1,
                                   measures = "T1_CvM_Cs3") })

  # It works for two measures
  expect_no_error({
    measures_nonsimplifyingness_NP(X1 = X1, X2 = X2, Z = Z, h = 1,
                                   measures = c("T1_CvM_Cs3", "T1_CvM_Cs4")) })
})

