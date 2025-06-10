#' Standard Error of Model Coefficients
#'
#' Generic function to calculate the standard error of estimated coefficients
#' from a model object.
#'
#' @param object A model object.
#' @param ... Additional arguments passed to methods.
#'
#' @export
se <- function(object, ...) {
  UseMethod("se")
}
