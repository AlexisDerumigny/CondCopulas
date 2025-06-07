

#' Compute a measure of non-simplifyingness based on non-parametric estimation
#' of the conditional copula
#'
#' @param X1,X2 vector of \code{n} observations of the conditioned variables
#'
#' @param Z vector of \code{n} observations of the conditioning variable
#'
#' @param measures choices of measures of non-simplifyingness to be computed.
#' \code{measures = "all"} includes all available non-parametric measures of
#' non-simplifyingness. Otherwise, \code{measures} must be a character vector
#' and a subset of \code{c("T1_CvM_Cs3", "T1_CvM_Cs4", "tilde_T0_CvM", }
#' \code{"T1_KS_Cs3", "T1_KS_Cs4", "tilde_T0_KS")}.
#'
#' @param h the bandwidth used for kernel smoothing
#'
#' @param kernel.name the name of the kernel
#'
#' @param truncVal the value of truncation for the integral,
#' i.e. the integrals are computed from \code{truncVal} to \code{1-truncVal}
#' instead of from 0 to 1.
#' Note that \code{truncVal} must be in the interval \eqn{[0, 0.5)},
#' i.e. \eqn{0} is allowed but not \eqn{0.5}.
#'
#' The default is \code{truncVal = NULL}, which actually means that
#' \code{truncVal = h} if \code{h < 0.5} and \code{truncVal = 0} else.
#'
#' @param numericalInt parameters to be given to
#' \code{statmod::\link[statmod]{gauss.quad}}, including the number of
#' quadrature points and the type of interpolation.
#'
#' @param verbose option used for debugging. If \code{verbose = 0}, the function
#' is silent. Higher values of \code{verbose} give more explicit details on the
#' computations.
#'
#'
#' @export
#'
measures_nonsimplifyingness_NP <- function(
    X1, X2, Z, h,
    measures = "all",
    kernel.name = "Epanechnikov", truncVal = NULL,
    numericalInt = list(kind = "legendre", nGrid = 10),
    verbose = 0)
{
  .checkSame_nobs_X1X2Z(X1, X2, Z)
  n = length(X1)

  .checkUnivX1X2Z(X1, X2, Z)
  if (length(truncVal) > 1){
    stop("'truncVal' must be of length 1 or NULL")
  }

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

  result = expand.grid(measure = measures,
                       h = h,
                       stringsAsFactors = FALSE)

  # We now determine the truncation value based on the user-specified `truncVal`
  # or on the `h` if `truncVal` is missing.
  if (is.null(truncVal)){
    result$truncVal = ifelse(result$h < 0.5, yes = result$h, no = 0)

  } else {
    if (truncVal < 0 || truncVal >= 0.5){
      stop(errorCondition(
        message = paste0("'truncVal' must be in the interval [0, 0.5). ",
                         "Here it is: ", truncVal),
        class = "InvalidInputError"
      ))
    }
    result$truncVal = truncVal
  }

  result$value = NA_real_

  env <- environment()

  # Necessary to make the other functions work...
  X3 = Z

  for (i in 1:nrow(result)){
    h = result$h[i]
    truncVal = result$truncVal[i]
    if (verbose > 0){
      cat("h = ", h, "; truncVal = ", truncVal, "; ")
    }

    nGrid = numericalInt$nGrid
    grid <- statmod::gauss.quad(n = nGrid, kind = numericalInt$kind)
    # Change of range to be on [truncVal , 1 - truncVal]
    grid$nodes <- grid$nodes * (1/2 - truncVal) + 1/2

    # FIXME: the weights should depend on `truncVal`
    grid$weights <- grid$weights / 2

    result$value[i] <- switch(
      result$measure[i],

      "T1_CvM_Cs3" = {testStat_T1_CvM_Cs3(env); env$true_stat},

      "T1_CvM_Cs4" = {testStat_T1_CvM_Cs4(env); env$true_stat},

      "tilde_T0_CvM" = {testStat_tilde_T0_CvM(env); env$true_stat},

      "T1_KS_Cs3" = {testStat_T1_KS_Cs3(env); env$true_stat},

      "T1_KS_Cs4" = {testStat_T1_KS_Cs4(env); env$true_stat},

      "tilde_T0_KS" = {testStat_tilde_T0_KS(env); env$true_stat}
    )
    if (verbose > 0){
      cat(result$measure[i], " = ", result$value[i], "\n")
    }
  }

  return (result)
}

