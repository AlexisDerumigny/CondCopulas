

#' Compute a measure of non-simplifyingness based on non-parametric estimation
#' of the conditional copula
#'
#' @param X1,X2 vector of \code{n} observations of the conditioned variables
#'
#' @param Z vector of \code{n} observations of the conditioning variable
#'
#'
#' @export
#'
measures_nonsimplifyingness_NP <- function(
    X1, X2, Z, h,
    measures = "all",
    kernel.name = "Epanechnikov", truncVal = h,
    numericalInt = list(kind = "legendre", nGrid = 10))
{
  .checkSame_nobs_X1X2Z(X1, X2, Z)
  n = length(X1)

  .checkUnivX1X2Z(X1, X2, Z)

  nGrid = numericalInt$nGrid
  grid <- statmod::gauss.quad(n = nGrid, kind = numericalInt$kind)
  # Change of range to be on [0,1]
  grid$nodes <- grid$nodes * (1/2 - h) + 1/2
  grid$weights <- grid$weights / 2

  possible_measures = c("T1_CvM_Cs3", "T1_CvM_Cs4", "tilde_T0_CvM",
                        "T1_KS_Cs3", "T1_KS_Cs4", "tilde_T0_KS")

  if (length(measures) == 0){
    stop(errorCondition(
      message = paste0("'measures' should not be of length 0."),
      class = "ZeroLengthError") )
  } else if (length(measures) == 1 && measures == "all"){
    measures = possible_measures
  } else {
    which_bad = which(!(measures %in% possible_measures))
    if(length(which_bad) > 0){
      stop(errorCondition(
        message = paste0(
          "Unknwon measure(s): ",
          "'", paste0(measures[which_bad], collapse = "', '"), "'.",
          "\n",
          "Possible measures are: ",
          "'", paste0(possible_measures, collapse = "', '"), "'."),
        class = "UnknownMeasureNameError" ) )
    }
  }

  result = data.frame(
    name = measures,
    value = NA_real_
  )

  env <- environment()

  # Necessary to make the other functions work...
  X3 = Z

  for (i in 1:length(measures)){
    result$value[i] <- switch(
      measures[i],

      "T1_CvM_Cs3" = {testStat_T1_CvM_Cs3(env); true_stat},

      "T1_CvM_Cs4" = {testStat_T1_CvM_Cs4(env); true_stat},

      "tilde_T0_CvM" = {testStat_tilde_T0_CvM(env); true_stat},

      "T1_KS_Cs3" = {testStat_T1_KS_Cs3(env); true_stat},

      "T1_KS_Cs4" = {testStat_T1_KS_Cs4(env); true_stat},

      "tilde_T0_KS" = {testStat_tilde_T0_KS(env); true_stat}
    )
  }

  return (result)
}

