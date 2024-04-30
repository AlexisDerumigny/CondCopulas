
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
