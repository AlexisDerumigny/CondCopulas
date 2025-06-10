
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

.checkSame_ncols_Z_newZ <- function(Z, newZ, name_Z, name_newZ){
  if (NCOL(Z) != NCOL(newZ)){
    stop(errorCondition(
      message = paste0(name_Z, " and ", name_newZ ,
                       " must have the same number of columns ",
                       "(so the same number of conditioning variables). ",
                       "However, here ", name_Z, " has ", NCOL(Z), " columns ",
                       "while ", name_newZ, " has ", NCOL(newZ), " columns." ),
      class = "WrongDimensionError")
    )
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

.check_MatrixSignPairs <- function(matrixSignsPairs)
{
  if (nrow(matrixSignsPairs) != ncol(matrixSignsPairs)){
    stop(errorCondition(
      message = paste0("'matrixSignsPairs' must be a square matrix. ",
                       "Here, its dimensions are: ",
                       nrow(matrixSignsPairs), " rows and ",
                       ncol(matrixSignsPairs), " columns."),
      class = "WrongDimensionError")
    )
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

#' Check whether an object Z is either a matrix or a data.frame or a vector
#' with only numeric components
#'
#' @returns a vector or a matrix with at least 2 columns, with the same content
#' as Z. This is guaranteed to be of type `numeric`.
#'
#' @noRd
.ensure_Z_numeric_vector_or_matrix <- function(Z, nameZ){
  if (is.vector(Z)){
    if(!is.numeric(Z)){
      stop(errorCondition(
        message = paste0("If ", nameZ, " is a vector, it should be numeric. ",
                         "Here, ", nameZ, " is of class ", class(Z) ,".") ,
        class = "NonNumericInputError"
      ))
    }
  } else if (inherits(Z, "data.frame")){
    Z = as.matrix.data.frame(Z)
    if(!is.numeric(Z)){
      stop(errorCondition(
        message = paste0(nameZ, " should be composed of numeric values. ",
                         "Here, ", nameZ, " is of storage mode ", mode(Z) ,".") ,
        class = "NonNumericInputError"
      ))
    }
  } else if (!inherits(Z, "matrix")){
    stop(errorCondition(
      message = paste0(nameZ, " should be a numeric matrix or vector.",
                       "Here, ", nameZ, " is of class ", class(Z) ,".") ,
      class = "NonNumericInputError"
    ))
  }

  if (NCOL(Z) == 1){
    Z = as.numeric(Z)
  }

  return (Z)
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
