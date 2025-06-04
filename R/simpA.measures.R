

#' Compute a measure of non-simplifyingness based on non-parametric estimation
#' of the conditional copula
#'
#' @param X1,X2 vector of \code{n} observations of the conditioned variables
#'
#' @param Z vector of \code{n} observations of the conditioning variable
#'
#'
#'
measures_nonsimplifyingness_NP <- function(
    X1, X2, Z, h,
    measures = "all",
    kernel.name = "Epanechnikov", truncVal = h,
    numericalInt = list(kind = "legendre", nGrid = 10))
{
  .checkSame_nobs_X1X2Z(X1, X2, Z)
  .checkUnivX1X2Z(X1, X2, Z)


}

