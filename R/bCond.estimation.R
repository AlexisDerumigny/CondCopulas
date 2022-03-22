
#' Computing the pseudo-observations in case of discrete
#' conditioning events
#'
#' @description
#' Let \eqn{A_1, ..., A_p} be \eqn{p} events forming a partition of
#' a probability space and \eqn{X_1, ..., X_d} be \eqn{d} random variables.
#' Assume that we observe \eqn{n} i.i.d. replications of \eqn{(X_1, ..., X_d)},
#' and that for each \eqn{i=1, ..., d},
#' \deqn{V_{i,j|A} = F_{X_j | A_k}(X_{i,j} | A_k),}
#' we also know which of the \eqn{A_k} was realized.
#' This function computes the pseudo-observations
#' where \eqn{k} is such that the event \eqn{A_k}
#' is realized for the \eqn{i}-th observation.
#'
#' @param X matrix of size \code{n * d} observations of conditioned variables.
#'
#' @param partition matrix of size \code{n * p},
#' where \code{p} is the number of conditioning events that are considered.
#' partition[i,k] should be the indicator of whether the \code{i}-th observation
#' belongs or not to the \code{k}-th conditioning event.
#'
#' @return a matrix of size \code{n * d}
#' containing the conditional pseudo-observations \eqn{V_{i,j|A}}.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2017).
#' About tests of the “simplifying” assumption for conditional copulas.
#' Dependence Modeling, 5(1), 154-197.
#'
#' @seealso \code{\link{bCond.estParamCopula}} for the estimation
#' of a (conditional) parametric copula model in this framework.
#'
#' \code{\link{bCond.treeCKT}} that provides a binary tree
#' based on conditional Kendall's tau
#' and that can be used to derive relevant conditioning events.
#'
#'
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
#' condPseudoObs = bCond.pobs(X = cbind(X1, X2),
#'                            partition = partition)
#'
#' @export
#'
bCond.pobs <- function(X, partition)
{
  if (nrow(X) != nrow(partition)){
    stop("The number of rows in 'partition' and in `X` ",
         "must be equal.")
  }

  n = nrow(partition)
  p = ncol(partition)
  d = ncol(X)

  matrix_VjA = matrix(nrow = n, ncol = d)

  for (box in 1:p)
  {
    indexesBoxes = which(partition[,box])
    for (j in 1:d)
    {
      matrix_VjA[indexesBoxes, j] = stats::ecdf(X[indexesBoxes,j]) (X[indexesBoxes,j])
    }
  }

  return (matrix_VjA)
}


#' Estimation of the conditional parameters of a parametric conditional
#' copula with discrete conditioning events.
#'
#'
#' By Sklar's theorem, any conditional distribution function
#' can be written as
#' \deqn{F_{1,2|A}(x_1, x_2) = c_{1,2|A}(F_{1|A}(x_1), F_{2,A}(x_2)),}
#' where \eqn{A} is an event and
#' \eqn{c_{1,2|A}} is a copula depending on the event \eqn{A}.
#' In this function, we assume that we have a partition \eqn{A_1,... A_p}
#' of the probability space, and that for each \eqn{k=1,...,p},
#' the conditional copula is parametric according to the following model
#' \deqn{c_{1,2|Ak} = c_{\theta(Ak)},}
#' for some parameter \eqn{\theta(Ak)} depending on the realized event \eqn{Ak}.
#' This function uses canonical maximum likelihood to estimate
#' \eqn{\theta(Ak)} and the corresponding copulas \eqn{c_{1,2|Ak}}.
#'
#'
#' @param U1 vector of \code{n} conditional pseudo-observations
#' of the first conditioned variable.
#'
#' @param U2 vector of \code{n} conditional pseudo-observations
#' of the second conditioned variable.
#'
#' @param family the family of conditional copulas
#' used for each conditioning event \eqn{A_k}. If not of length \eqn{p},
#' it is recycled to match the number of events \eqn{p}.
#'
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
#' @seealso \code{\link{bCond.pobs}} for the computation
#' of (conditional) pseudo-observations in this framework.
#'
#' \code{\link{bCond.simpA.param}} for a test of the simplifying assumption
#' that all these conditional copulas are equal
#' (assuming they all belong to the same parametric family).
#' \code{\link{bCond.simpA.CKT}} for a test of the simplifying assumption
#' that all these conditional copulas are equal,
#' based on the equality of conditional Kendall's tau.
#'
#'
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
#' condPseudoObs = bCond.pobs(X = cbind(X1, X2), partition = partition)
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

