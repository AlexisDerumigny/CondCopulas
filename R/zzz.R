
.checkSame_nobs_X1X2Z <- function(X1, X2, Z)
{
  if ((NROW(X1) != NROW(X2)) || (NROW(X1) != NROW(Z))){
    stop(errorCondition(
      message = paste0("X1, X2 and Z must have the same number of observations. ",
                       "Here they are respectively: ",
                       NROW(X1), ", ", NROW(X2), ", ", NROW(Z)),
      class = "DifferentLengthsError") )
  }
}

.checkSame_nobs_X1X2X3 <- function(X1, X2, X3)
{
  if ((NROW(X1) != NROW(X2)) || (NROW(X1) != NROW(X3))){
    stop(errorCondition(
      message = paste0("X1, X2 and X3 must have the same number of observations. ",
                       "Here they are respectively: ",
                       NROW(X1), ", ", NROW(X2), ", ", NROW(X3)),
      class = "DifferentLengthsError") )
  }
}

.checkUnivX1X2 <- function(X1, X2)
{
  if (NCOL(X1) > 1){
    stop(errorCondition(
      message = paste0("X1 must be univariate. Here it has ",
                       NCOL(X1), " columns"),
      class = "WrongDimensionError") )
  }
  if (NCOL(X2) > 1){
    stop(errorCondition(
      message = paste0("X2 must be univariate. Here it has ",
                       NCOL(X2), " columns"),
      class = "WrongDimensionError") )
  }
}

.checkUnivX1X2X3 <- function(X1, X2, X3)
{
  if (NCOL(X1) > 1){
    stop(errorCondition(
      message = paste0("X1 must be univariate. Here it has ",
                       NCOL(X1), " columns"),
      class = "WrongDimensionError") )
  }
  if (NCOL(X2) > 1){
    stop(errorCondition(
      message = paste0("X2 must be univariate. Here it has ",
                       NCOL(X2), " columns"),
      class = "WrongDimensionError") )
  }
  if (NCOL(X3) > 1){
    stop(errorCondition(
      message = paste0("X3 must be univariate. Here it has ",
                       NCOL(X3), " columns"),
      class = "WrongDimensionError") )
  }
}

.checkUnivX1X2Z <- function(X1, X2, Z)
{
  if (NCOL(X1) > 1){
    stop(errorCondition(
      message = paste0("X1 must be univariate. Here it has ",
                       NCOL(X1), " columns"),
      class = "WrongDimensionError") )
  }
  if (NCOL(X2) > 1){
    stop(errorCondition(
      message = paste0("X2 must be univariate. Here it has ",
                       NCOL(X2), " columns"),
      class = "WrongDimensionError") )
  }
  if (NCOL(Z) > 1){
    stop(errorCondition(
      message = paste0("Z must be univariate. Here it has ",
                       NCOL(Z), " columns"),
      class = "WrongDimensionError") )
  }
}

.observedX1X2_to_X1X2 <- function(env)
{
  if (is.null(env$X1)){
    if (is.null(env$observedX1)){
      stop("X1 must be non-null.")
    } else {
      env$X1 = env$observedX1
    }
  }
  if (is.null(env$X2)){
    if (is.null(env$observedX2)){
      stop("X2 must be non-null.")
    } else {
      env$X2 = env$observedX2
    }
  }
}

.observedZ_to_Z <- function(env)
{
  if (is.null(env$Z)){
    if (is.null(env$observedZ)){
      stop("Z must be non-null.")
    } else {
      env$Z = env$observedZ
    }
  }
}

.observedX3_to_X3 <- function(env)
{
  if (is.null(env$X3)){
    if (is.null(env$observedX3)){
      stop("X3 must be non-null.")
    } else {
      env$X3 = env$observedX3
    }
  }
}


#' Constructor for warning conditions of the package
#'
#' @noRd
CondCopulas_warning_condition_base <- function(message, subclass = NULL, call = sys.call(-1), ...) {
  # warningCondition() automatically adds 'warning' and 'condition' to the class
  return (
    warningCondition(
      message = message,
      class = c(subclass, "CondCopulasWarning"), # We add a base warning class
      call = call,
      ... # Allows for additional custom fields
    )
  )
}
