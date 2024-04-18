

#' Test of the simplifying assumption using the constancy
#' of conditional Kendall's tau
#'
#' This function computes Kendall's regression, a regression-like
#' model for conditional Kendall's tau. More precisely, it fits the model
#' \deqn{\Lambda(\tau_{X_1, X_2 | Z = z}) = \sum_{j=1}^{p'} \beta_j \psi_j(z),}
#' where \eqn{\tau_{X_1, X_2 | Z = z}} is the conditional Kendall's tau
#' between \eqn{X_1} and \eqn{X_2} conditionally to \eqn{Z=z},
#' \eqn{\Lambda} is a function from \eqn{]-1, 1]} to \eqn{R},
#' \eqn{(\beta_1, \dots, \beta_p)} are unknown coefficients to be estimated
#' and \eqn{\psi_1, \dots, \psi_{p'})} are a dictionary of functions.
#' Then, this function tests the assumption
#' \deqn{\beta_2 = \beta_3 = ... = \beta_{p'} = 0,}
#' where the coefficient corresponding to the intercept is removed.
#'
#' @param X1 vector of observations of the first conditioned variable
#'
#' @param X2 vector of observations of the second conditioned variable
#'
#' @param Z vector of observations of the conditioning variable
#'
#' @param vectorZToEstimate vector containing the points \eqn{Z'_i}
#' to be used at which the conditional Kendall's tau should be estimated.
#'
#' @param listPhi the list of transformations \eqn{phi} to be used.
#'
#' @param typeEstCKT the type of estimation of the kernel-based estimation
#' of conditional Kendall's tau.
#'
#' @param h_kernel the bandwidth used for the kernel-based estimations.
#'
#' @param Lambda the function to be applied on conditional Kendall's tau.
#' By default, the identity function is used.
#'
#' @param Lambda_deriv the derivative of the function \code{Lambda}.
#'
#' @param lambda the penalization parameter used for Kendall's regression.
#' By default, cross-validation is used to find the best value of \code{lambda}
#' if \code{length(listPhi) > 1}. Otherwise \code{lambda = 0} is used.
#'
#' @param h_lambda bandwidth used for the smooth cross-validation
#' in order to get a value for \code{lambda}.
#'
#' @param Kfolds_lambda the number of subsets used for the cross-validation
#' in order to get a value for \code{lambda}.
#'
#' @param l_norm type of norm used for selection of the optimal lambda
#' by cross-validation. \code{l_norm=1} corresponds to the sum of
#' absolute values of differences between predicted and estimated
#' conditional Kendall's tau while \code{l_norm=2} corresponds to
#' the sum of squares of differences.
#'
#'
#' @return a list containing
#' \itemize{
#'     \item \code{statWn}: the value of the test statistic.
#'     \item \code{p_val}: the p-value of the test.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2020).
#' On Kendall’s regression.
#' Journal of Multivariate Analysis, 178, 104610.
#' (page 7)
#' \doi{10.1016/j.jmva.2020.104610}
#'
#' @seealso The function to fit Kendall's regression:
#' \code{\link{CKT.kendallReg.fit}}.
#'
#' Other tests of the simplifying assumption:
#' \itemize{
#'   \item \code{\link{simpA.NP}} in a nonparametric setting
#'   \item \code{\link{simpA.param}} in a (semi)parametric setting,
#'   where the conditional copula belongs to a parametric family,
#'   but the conditional margins are estimated arbitrarily through
#'   kernel smoothing
#'
#'   \item the counterparts of these tests in the discrete conditioning setting:
#'   \code{\link{bCond.simpA.CKT}}
#'   (test based on conditional Kendall's tau)
#'   \code{\link{bCond.simpA.param}}
#'   (test assuming a parametric form for the conditional copula)
#' }
#'
#'
#' @examples
#'
#' \donttest{
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 300
#' Z = runif(n = N, min = 0, max = 1)
#' conditionalTau = -0.9 + 1.8 * Z
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' result = simpA.kendallReg(
#'    X1, X2, Z, h_kernel = 0.03,
#'    listPhi = list(
#'      function(x){return(x)}, function(x){return(x^2)},
#'      function(x){return(x^3)}, function(x){return(x^4)},
#'      function(x){return(x^5)},
#'      function(x){return(cos(x))}, function(x){return(sin(x))},
#'      function(x){return(as.numeric(x <= 0.4))},
#'      function(x){return(as.numeric(x <= 0.6))}) )
#' print(result$p_val)
#'
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 300
#' Z = runif(n = N, min = 0, max = 1)
#' conditionalTau = -0.3
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' result = simpA.kendallReg(
#'    X1, X2, Z, h_kernel = 0.03,
#'    listPhi = list(
#'      function(x){return(x)}, function(x){return(x^2)},
#'      function(x){return(x^3)}, function(x){return(x^4)},
#'      function(x){return(x^5)},
#'      function(x){return(cos(x))}, function(x){return(sin(x))},
#'      function(x){return(as.numeric(x <= 0.4))},
#'      function(x){return(as.numeric(x <= 0.6))}) )
#' print(result$p_val)
#' }
#'
#' @export
#'
simpA.kendallReg <- function(
    X1, X2, Z,
    vectorZToEstimate = NULL,
    listPhi = list(function(x){return(x)}, function(x){return(x^2)},
                   function(x){return(x^3)}),
    typeEstCKT = 4,
    h_kernel,
    Lambda = function(x){return(x)}, Lambda_deriv = function(x){return(1)},
    lambda = NULL, h_lambda = h_kernel,
    Kfolds_lambda = 5, l_norm = 1
)
{
  if (is.null(vectorZToEstimate)){
    nprime = 100
    vectorZToEstimate = computeVectorZToEstimate(
      vecteurZrealised = Z, length.out = nprime,
      percentage = 0.95)
  }

  matrixSignsPairs = computeMatrixSignPairs(X1, X2, typeEstCKT = typeEstCKT)

  # Computation of the \hat \tau_{1,2 | \Z = \Z'_i}
  vectorEstimate_1step = CKT.kernel.univariate(
    matrixSignsPairs = matrixSignsPairs ,
    observedZ = Z ,
    h = h_kernel ,
    ZToEstimate = vectorZToEstimate ,
    kernel.name = "Epa" ,
    typeEstCKT = typeEstCKT )

  LambdaCKT = Lambda(vectorEstimate_1step)

  whichFinite = which( is.finite(LambdaCKT))
  if (is.null(whichFinite)) {
    stop("No kernel estimation successful. ",
         "Maybe h_kernel is too small?")
  }

  # To compute \hat \beta,
  # we first prepare the design matrix and the penalization parameter lambda

  # Computation of the design matrix \Zb by applying each phi to Z'_i
  designMatrixZ = sapply(listPhi,
                         function(phi) {phi(vectorZToEstimate)},
                         simplify = "array")

  # Choice of lambda
  if (length(listPhi) == 1){
    if (!is.null(lambda)){
      if (lambda != 0){
        warning(paste0("There is only one function 'phi' provided. ",
                       "Therefore the penalization parameter lambda is set to 0."))
      }
    }
    lambda = 0
  } else {

    # Choice of lambda by cross-validation if it is `NULL`
    if (is.null(lambda)){
      # Computation of \lambda chosen by cross-validation
      resultCV = CKT.KendallReg.LambdaCV(
        observedX1 = X1, observedX2 = X2,
        observedZ = Z, ZToEstimate = vectorZToEstimate,
        designMatrixZ = designMatrixZ,
        typeEstCKT = typeEstCKT,
        h_lambda = h_lambda, Lambda = Lambda, kernel.name = "Epa",
        Kfolds_lambda = Kfolds_lambda, l_norm = l_norm)

      lambda = resultCV$lambdaCV
    }
  }

  # Estimation
  if (lambda == 0){
    # Estimation by least-squares
    designMatrix_withIntercept = cbind(1 , designMatrixZ)
    colnames(d)[1] <- "Intercept"

    reg = stats::lm.fit(x = designMatrix_withIntercept[whichFinite, ],
                        y = LambdaCKT[whichFinite])
    vector_hat_beta = reg$coefficient
  } else {
    # Penalized estimation
    reg = glmnet::glmnet(x = designMatrixZ[whichFinite, ],
                         y = LambdaCKT[whichFinite],
                         family = "gaussian")

    vector_hat_beta = stats::coef(reg, s = lambda)
  }

  # Using Wn
  progressBar = TRUE #### TODO : manage the progressbars

  resultWn = computeWn(vectorZ = Z,
                       vectorZToEstimate = vectorZToEstimate[whichFinite],
                       vector_hat_CKT_NP = vectorEstimate_1step[whichFinite],
                       vector_hat_beta = vector_hat_beta,
                       matrixSignsPairs = matrixSignsPairs,
                       inputMatrix = designMatrixZ[whichFinite, , drop = FALSE],
                       h = h_kernel,
                       kernel.name = "Epa",
                       intK2 = 3/5,
                       Lambda_deriv = Lambda_deriv,
                       progressBar = progressBar)

  df = length(whichFinite)
  pval_Wn = 1 - stats::pchisq(as.numeric(resultWn$W_n), df = df)

  return (list(p_val = pval_Wn, statWn = resultWn$W_n, df = df,
               coef = vector_hat_beta,
               resultWn = resultWn))
}



# Compute the Wald-type test statistic related to
# the hypothesis $\beta = 0$ against $\beta != 0$
# vectorZ is the vector of Z in the database of length n
# vectorZToEstimate is the vector of z'_i
# vector_hat_CKT_NP is the vector of the nonparametric estimators of the CKT at every z'_i
# vector_hat_beta is the vector of coefficients, of size p'
# matrixSignsPairs is a matrix of size n * n
# inputMatrix is the inputMatrix of size n' * p'
computeWn = function(
    vectorZ, vectorZToEstimate, vector_hat_CKT_NP, vector_hat_beta,
    matrixSignsPairs, inputMatrix, h, kernel.name, intK2, Lambda_deriv,
    progressBar)
{
  nprime = length(vectorZToEstimate)
  n = length(vectorZ)
  pprime = length(inputMatrix[1,])

  if (length(vector_hat_CKT_NP) != nprime) {
    stop("vector_hat_CKT_NP and vectorZToEstimate ",
         "do not have the same length : ",
         length(vector_hat_CKT_NP), " and ", nprime)
  } else if ( (dim(matrixSignsPairs)[1] != n) || (dim(matrixSignsPairs)[2] != n) ) {
    stop("Wrong dimensions for matrixSignsPairs : ",
         paste0(dim(matrixSignsPairs), collapse =" "), " instead of ", n, " ", n)
  } else if ( length(inputMatrix[,1]) !=  nprime) {
    stop("Wrong nrow for inputMatrix : ",
         length(inputMatrix[,1]), " instead of ", nprime)
  }

  # 1 - Computation of G_n(z'_i) for all i
  vectorZToEstimate_arr = array(vectorZToEstimate)
  matrixSignsPairsSymmetrized = (matrixSignsPairs + t(matrixSignsPairs)) / 2
  if (progressBar){
    Gn_zipr = pbapply::pbapply(
      X = vectorZToEstimate_arr, MARGIN = 1,
      FUN = function(pointZ) {compute_vect_Gn_zipr(
        pointZ, vectorZ, h,
        kernel.name, matrixSignsPairsSymmetrized) } )
  } else {
    Gn_zipr = apply(
      X = vectorZToEstimate_arr, MARGIN = 1,
      FUN = function(pointZ) {compute_vect_Gn_zipr(
        pointZ, vectorZ, h,
        kernel.name, matrixSignsPairsSymmetrized) } )
  }

  # 2 - Computation of H_(i,i) under the hypothesis that all z'_i are distinct
  vect_H_ii = rep(NA, nprime)
  for (iprime in 1:nprime) {
    pointZ = vectorZToEstimate[iprime]
    listKh = computeWeights.univariate(
      vectorZ = vectorZ, h = h, pointZ = pointZ,
      kernel.name = kernel.name, normalization = FALSE)
    estimator_fZ = mean(listKh)
    vect_H_ii[iprime] = 4 * (intK2 / estimator_fZ) *
      (Lambda_deriv(pointZ))^2 * (Gn_zipr[iprime] - (vector_hat_CKT_NP[iprime])^2)
  }

  # 3 - Computation of Sigma_npr
  matrix_Sigma_npr = matrix(rep(0, pprime^2), nrow = pprime, ncol = pprime)
  for (iprime in 1:nprime) {
    matrix_Sigma_npr = matrix_Sigma_npr +
      ( t( t( inputMatrix[iprime, ] ) ) %*% t( inputMatrix[iprime, ] ) )
  }
  inv_matrix_Sigma_npr = solve(matrix_Sigma_npr)

  # 4 - Computation of V_n
  matrix_Vn = matrix(rep(0, pprime^2), nrow = pprime, ncol = pprime)
  for (iprime in 1:nprime) {
    matrix_Vn = matrix_Vn + vect_H_ii[iprime] *
      ( t( t( inputMatrix[iprime, ] ) ) %*% t( inputMatrix[iprime, ] ) )
  }
  matrix_Vn_completed =
    inv_matrix_Sigma_npr %*% matrix_Vn %*% inv_matrix_Sigma_npr

  # 5 - Computation of W_n
  W_n = n * h * t(vector_hat_beta[-1]) %*%
    matrix_Vn_completed %*% t( t(vector_hat_beta[-1]) )

  return (list(W_n = W_n, Gn_zipr = Gn_zipr,
               vect_H_ii = vect_H_ii,
               matrix_Sigma_npr = matrix_Sigma_npr,
               matrix_Vn = matrix_Vn,
               matrix_Vn_completed = matrix_Vn_completed))
}


# Compute a series of points to estimate
# vecteurZrealised is a vector of observations of Z
# length.out is the length of the output
# percentage is the percentage of "coverage"
# of the interval from which the new Z are given.
computeVectorZToEstimate = function(vecteurZrealised,
                                    length.out, percentage)
{
  minZ = min(vecteurZrealised)
  maxZ = max(vecteurZrealised)
  fracExcluded = (1-percentage)/2 # on the left and on the right of the interval
  newMinZ = minZ + fracExcluded*(maxZ-minZ)
  newMaxZ = maxZ - fracExcluded*(maxZ-minZ)
  return (seq(newMinZ, newMaxZ, length.out = length.out))
}


# This function returns Gn(z'_i), z'_i = pointZ
# which is an estimator of the conditional expectation
# E(g(X_1, X_2) g(X_1, X_3) | Z_1 = Z_2 = Z_3 = z'_i]
compute_vect_Gn_zipr = function(pointZ, vectorZ, h,
                                kernel.name, matrixSignsPairsSymmetrized)
{
  n = length(vectorZ)

  listWeights = computeWeights.univariate(
    vectorZ = vectorZ, h = h, pointZ = pointZ,
    kernel.name = kernel.name, normalization = TRUE)
  matrixWeights = outer(listWeights, listWeights)
  result = 0
  for (k in 1:n) {
    result = result + listWeights[k] *
      (matrixSignsPairsSymmetrized[k,-k] %*%
         matrixWeights[-k,-k] %*%
         matrixSignsPairsSymmetrized[-k,k])
  }
  return(result)

  # for (iprime in 1:nprime)
  # {
  # pointZ = vectorZToEstimate[iprime]
  # listWeights = computeWeights(vectorZ, h, pointZ, kernel.name)
  # matrixWeights = outer(listWeights, listWeights)
  # Gn_zipr[iprime] = 0
  # for (k in 1:n)
  # {
  # Gn_zipr[iprime] = Gn_zipr[iprime] +
  # listWeights[k] * (matrixSignsPairsSymmetrized[k,-k] %*%
  # matrixWeights[-k,-k] %*% matrixSignsPairsSymmetrized[-k,k])
  # }
  # }
}


