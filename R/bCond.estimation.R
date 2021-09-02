
#' Computing the pseudo-observations in case of discrete
#' conditioning events
#'
#' @param X1 vector of \code{n} observations of the first conditioned variable.
#' @param X2 vector of \code{n} observations of the second conditioned variable.
#' @param partition matrix of size \code{n * p},
#' where \code{p} is the number of conditioning events that are considered.
#' partition[i,j] should be the indicator of whether the \code{i}-th observation
#' belongs or not to the \code{j}-th conditioning event.
#'
#' @return a matrix containing the conditional pseudo-observations.
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#'
#' @examples
#' n = 800
#' Z = stats::runif(n = n)
#' CKT = 0.2 * as.numeric(Z <= 0.3) +
#'   0.5 * as.numeric(Z > 0.3 & Z <= 0.5) +
#'   - 0.8 * as.numeric(Z > 0.5)
#' simCopula = VineCopula::BiCopSim(N = n,
#'   par = VineCopula::BiCopTau2Par(CKT, family = 1), family = 1)
#' X1 = simCopula[,1]
#' X2 = simCopula[,2]
#' partition = cbind(Z <= 0.3, Z > 0.3 & Z <= 0.5, Z > 0.5)
#' condPseudoObs = bCond.pobs(X1 = X1, X2 = X2, partition = partition)
#'
#' @export
#'
bCond.pobs <- function(X1, X2, partition)
{
  if (length(X1) != length(X2)){stop("X1 and X2 should be of the same length.")}
  if (length(X1) != nrow(partition)){
    stop("X1 should have the same length as the number of rows in 'partition'")
  }
  n = nrow(partition)
  p = ncol(partition)
  V1_A = rep(NA , n)
  V2_A = rep(NA , n)

  for (box in 1:p)
  {
    indexesBoxes = which(partition[,box])
    V1_A[indexesBoxes] = stats::ecdf(X1[indexesBoxes]) (X1[indexesBoxes])
    V2_A[indexesBoxes] = stats::ecdf(X2[indexesBoxes]) (X2[indexesBoxes])
  }

  return (cbind(V1_A , V2_A))
}


#' Estimation of the conditional parameters of a parametric conditional
#' copula with discrete conditioning events.
#'
#'
#' This function uses canonical maximum likelihood to estimate
#'
#' This function is currently implemented only for one-parameter families
#' of conditional copulas and for the Student family with fixed degree of freedom.
#'
#' @param U1 vector of \code{n} conditional pseudo-observations of the first conditioned variable
#' @param U2 vector of \code{n} conditional pseudo-observations of the second conditioned variable
#' @param family the family of conditional copulas used.
#' Can be a number or a vector of size p
#' @param partition matrix of size \code{n * p},
#' where \code{p} is the number of conditioning events that are considered.
#' partition[i,j] should be the indicator of whether the \code{i}-th observation
#' belongs or not to the \code{j}-th conditioning event
#'
#' @return a list of size \code{p} containing the \code{p} conditional copulas
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#'
#' @examples
#' n = 800
#' Z = stats::runif(n = n)
#' CKT = 0.2 * as.numeric(Z <= 0.3) +
#'   0.5 * as.numeric(Z > 0.3 & Z <= 0.5) +
#'   - 0.8 * as.numeric(Z > 0.5)
#' simCopula = VineCopula::BiCopSim(N = n,
#'   par = VineCopula::BiCopTau2Par(CKT, family = 1), family = 1)
#' X1 = simCopula[,1]
#' X2 = simCopula[,2]
#' partition = cbind(Z <= 0.3, Z > 0.3 & Z <= 0.5, Z > 0.5)
#' condPseudoObs = bCond.pobs(X1 = X1, X2 = X2, partition = partition)
#'
#' estimatedCondCopulas = bCond.estParamCopula(
#'   U1 = condPseudoObs[,1], U2 = condPseudoObs[,2],
#'   family = 1, partition = partition)
#' print(estimatedCondCopulas)
#' # Comparison with the true conditional parameters: 0.2, 0.5, -0.8.
#'
#'
#' @export
#'
bCond.estParamCopula <- function(U1, U2, family, partition)
{
  if (length(U1) != length(U2)){stop("U1 and U2 should be of the same length.")}
  if (length(U1) != nrow(partition)){
    stop("U1 should have the same length as the number of rows in 'partition'")
  }
  p = ncol(partition)
  family = rep(family, length.out = p)
  copulas_boxes = as.list(rep(NA , p))

  for (box in 1:p)
  {
    if (family[box] != 2) {
      copulas_boxes[[box]] =
        try(VineCopula::BiCopEst(u1 = U1[which(partition[,box])] ,
                                 u2 = U2[which(partition[,box])] ,
                                 family = family[box], method = "mle") , silent = TRUE)
    } else if (family[box] == 2) {
      copulas_boxes[[box]] =
        try(VineCopula::BiCopEst(u1 = U1[which(partition[,box])] ,
                                 u2 = U2[which(partition[,box])] ,
                                 family = 1, method = "itau")$par , silent = TRUE)
    }
  }

  return (copulas_boxes)
}

