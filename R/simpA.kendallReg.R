

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
#' @param Lambda_inv the inverse function of \code{Lambda}.
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
#' @param x an \code{S3} object of class \code{simpA_kendallReg_test}.
#' @param ylim graphical parameter, see \link{plot}
#'
#'
#' @return \code{simpA.kendallReg} returns an \code{S3} object of
#' class \code{simpA_kendallReg_test}, containing
#' \itemize{
#'     \item \code{statWn}: the value of the test statistic.
#'     \item \code{p_val}: the p-value of the test.
#' }
#'
#' \code{plot.simpA_kendallReg_test} returns (invisibly) a matrix with columns
#' \code{z}, \code{est_CKT_NP}, \code{asympt_se_np}, \code{est_CKT_NP_q025},
#' \code{est_CKT_NP_q975}, \code{est_CKT_reg}, \code{asympt_se_reg},
#' \code{est_CKT_reg_q025}, code{est_CKT_reg_q975}.
#' The first column correspond to the grid of values of z. The next 4 columns
#' are the NP (kernel-based) estimator of conditional Kendall's tau, with its
#' standard error, and lower/upper confidence bands. The last 4 columns are the
#' equivalents for the estimator based on Kendall's regression.
#'
#' \code{plot.simpA_kendallReg_test} plots the kernel-based estimator and its
#' confidence band (in red), and the estimator based on Kendall's regression
#' and its confidence band (in blue).
#'
#' Usually the confidence band for Kendall's regression is much tighter than the
#' pure non-parametric counterpart. This is because the parametric model is
#' sparser and the corresponding estimator converges faster (even without
#' penalization).
#'
#' \code{print.simpA_kendallReg_test} has no return values and is only called
#' for its side effects.
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
#' # We simulate from a non-simplified conditional copula
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
#'   X1, X2, Z, h_kernel = 0.03,
#'   listPhi = list(z = function(z){return(z)} ) )
#' print(result)
#' plot(result)
#'
#' result_morePhi = simpA.kendallReg(
#'    X1, X2, Z, h_kernel = 0.03,
#'    listPhi = list(
#'      z = function(z){return(z)},
#'      cos10z = function(z){return(cos(10 * z))},
#'      sin10z = function(z){return(sin(10 * z))},
#'      `1(z <= 0.4)` = function(z){return(as.numeric(z <= 0.4))},
#'      `1(z <= 0.6)` = function(z){return(as.numeric(z <= 0.6))}) )
#' print(result_morePhi)
#' plot(result_morePhi)
#'
#' # We simulate from a simplified conditional copula
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
#'      z = function(z){return(z)},
#'      cos10z = function(z){return(cos(10 * z))},
#'      sin10z = function(z){return(sin(10 * z))},
#'      `1(z <= 0.4)` = function(z){return(as.numeric(z <= 0.4))},
#'      `1(z <= 0.6)` = function(z){return(as.numeric(z <= 0.6))}) )
#' print(result)
#' plot(result)
#' }
#'
#' @export
#'
simpA.kendallReg <- function(
    X1, X2, Z,
    vectorZToEstimate = NULL,
    listPhi = list(z = function(z){return(z)}),
    typeEstCKT = 4,
    h_kernel,
    Lambda = function(x){return(x)},
    Lambda_deriv = function(x){return(1)},
    Lambda_inv = function(x){return(x)},
    lambda = NULL, h_lambda = h_kernel,
    Kfolds_lambda = 5, l_norm = 1
)
{
  if ((NROW(X1) != NROW(X2)) || (NROW(X1) != NROW(Z))){
    stop(errorCondition(
      message = paste0("X1, X2 and Z should have the same number of observations. ",
                       "Here they are respectively: ",
                       NROW(X1), ", ", NROW(X2), ", ", NROW(Z)),
      class = "DifferentLengthsError") )
  }
  if (NCOL(X1) > 1){
    stop(errorCondition(
      message = paste0("X1 should be univariate. Here it has ",
                       NCOL(X1), " columns"),
      class = "WrongDimensionError") )
  }
  X1 = as.numeric(X1)
  if (NCOL(X2) > 1){
    stop(errorCondition(
      message = paste0("X2 should be univariate. Here it has ",
                       NCOL(X2), " columns"),
      class = "WrongDimensionError") )
  }
  X2 = as.numeric(X2)
  if (NCOL(Z) > 1){
    stop("'simpA.kendallReg' is currently only implemented for univariate Z")
  }
  Z = as.numeric(Z)

  n = length(X1)
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
    typeEstCKT = typeEstCKT,
    progressBar = 0)

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

  designMatrix_withIntercept = cbind(1 , designMatrixZ)
  if (is.null(names(listPhi))){
    coefNames = c("Intercept", paste0("phi", 1:length(listPhi)))
  } else {
    coefNames = c("Intercept", names(listPhi))
  }
  colnames(designMatrix_withIntercept) <- coefNames

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
    reg = stats::lm.fit(x = designMatrix_withIntercept[whichFinite, ],
                        y = LambdaCKT[whichFinite])
    vector_hat_beta = reg$coefficient
    fitted.values = reg$fitted.values
  } else {
    # Penalized estimation
    reg = glmnet::glmnet(x = designMatrix_withIntercept[whichFinite, ],
                         y = LambdaCKT[whichFinite],
                         family = "gaussian",
                         intercept = FALSE)

    vector_hat_beta = stats::coef(reg, s = lambda)[-1]
    names(vector_hat_beta) <- coefNames
    fitted.values = stats::predict(reg, s = lambda,
                                   newx = designMatrix_withIntercept[whichFinite, ])
  }

  # non_zeroCoeffs = which(vector_hat_beta != 0)
  # if (length(non_zeroCoeffs) == 0){
  #   warning("All coefficients are found to be equal to 0. ",
  #           "This means that lambda is probably too high. Currently lambda = ",
  #           lambda)
  #
  #   return (NULL)
  # }

  # Computation of the asymptotic variance matrix

  resultWn = compute_for_Wn(
    vectorZ = Z,
    vectorZToEstimate = vectorZToEstimate[whichFinite],
    vector_hat_CKT_NP = vectorEstimate_1step[whichFinite],
    vector_hat_beta = vector_hat_beta,
    matrixSignsPairs = matrixSignsPairs,
    inputMatrix = designMatrix_withIntercept[whichFinite, ],
    h = h_kernel,
    kernel.name = "Epa",
    intK2 = 3/5,
    Lambda_deriv = Lambda_deriv,
    progressBar = TRUE #### TODO : manage the progressbars
  )

  # computation of the variance-covariance matrix
  varCov = resultWn$matrix_Vn_completed / (n * h_kernel)

  df = length(listPhi)

  # Computation of the test statistic W_n
  statWn = n * h_kernel * t(vector_hat_beta[-1]) %*%
    solve(resultWn$matrix_Vn_completed[-1, -1]) %*% t( t(vector_hat_beta[-1]) )

  # Conversion to numeric type and computation of the p-value
  statWn = as.numeric(statWn)
  pval_Wn = 1 - stats::pchisq(statWn, df = df)

  result = list(p_val = pval_Wn,
                statWn = statWn, df = df,
                coef = vector_hat_beta,
                resultWn = resultWn,
                varCov = varCov,
                designMatrix_withIntercept = designMatrix_withIntercept,
                vectorZToEstimate = vectorZToEstimate,
                vector_hat_CKT_NP = vectorEstimate_1step,
                n = n,
                h_kernel = h_kernel,
                fitted.values = fitted.values,
                Lambda = Lambda,
                Lambda_inv = Lambda_inv,
                Lambda_deriv = Lambda_deriv,
                lambda = lambda)

  class(result) <- "simpA_kendallReg_test"
  return (result)
}


#' Print method
#'
#' @param x object of class \code{simpA_kendallReg_test}
#'
#' @param ... other arguments, unused
#'
#' @export
#' @rdname simpA.kendallReg
print.simpA_kendallReg_test <- function(x, ...)
{
  # Use stringi::stri_escape_unicode to get the unicode characters escaped
  cat("Kendall regression: \u039b(\U0001d70f) = \u03b20 + \u03b2\' phi(z) \n")
  cat("where \U0001d70f is conditional Kendall's tau between X1 and X2 given Z = z \n \n")
  cat("Coefficients: \n")
  std_errors = diag(x$varCov)
  z_values = x$coef / diag(x$varCov)
  coef = cbind(Estimate = x$coef,
               `Std. Error` = std_errors,
               `z value` = z_values,
               `p-value` = 2 * stats::pnorm(abs(z_values), lower.tail = FALSE))

  stats::printCoefmat(coef)
  cat("\n")
  cat("Wald statistics \u03b2\' V^{-1} \u03b2\ :", x$statWn, "\n")
  cat("Test of the simplifying assumption with", x$df, "DF,  ",
      "pvalue: ", x$p_val)
  cat("\n")
}


#' Plot method
#'
#' @param x object of class \code{simpA_kendallReg_test}
#'
#' @param ... other arguments, unused
#'
#' @export
#' @rdname simpA.kendallReg
plot.simpA_kendallReg_test <- function(x, ylim = c(-1.5, 1.5), ...)
{
  z = x$vectorZToEstimate
  est_CKT_NP = x$vector_hat_CKT_NP
  asympt_se_np = sqrt(x$resultWn$vect_H_ii) / sqrt(x$n * x$h_kernel)

  plot(z, est_CKT_NP, type = "l", ylim = ylim, col = "red",
       ylab = "Conditional Kendall's tau given Z = z")
  graphics::lines(z, est_CKT_NP + 1.96 * asympt_se_np, type = "l", col = "red")
  graphics::lines(z, est_CKT_NP - 1.96 * asympt_se_np, type = "l", col = "red")

  est_CKT_reg = vapply(X = x$fitted.values, FUN = x$Lambda_inv, FUN.VALUE = 1)

  # Remember the chain rule (lambda^{-1})' = 1 / (lambda' o lambda^{-1}))
  lambdainvprime_estCKT = 1 / vapply(X = x$fitted.values, FUN = x$Lambda_deriv, FUN.VALUE = 1)

  # vcov matrix for beta' phi(z)
  designMatrix = x$designMatrix_withIntercept
  varCov = x$varCov
  vcova = designMatrix %*% varCov %*% t(designMatrix)

  # vcov matrix for lambdainv(beta' phi(z))
  # by the delta method, this is given by the sandwich rule
  vcova = diag(lambdainvprime_estCKT) %*% vcova %*% diag(lambdainvprime_estCKT)

  asympt_se_reg = sqrt(diag(vcova))

  graphics::lines(z, est_CKT_reg, type = "l", ylim = ylim, col = "blue")
  graphics::lines(z, est_CKT_reg + 1.96 * asympt_se_reg, type = "l", col = "blue")
  graphics::lines(z, est_CKT_reg - 1.96 * asympt_se_reg, type = "l", col = "blue")

  df = cbind(
    z,
    est_CKT_NP,
    asympt_se_np,
    est_CKT_NP_q025 = est_CKT_NP - 1.96 * asympt_se_np,
    est_CKT_NP_q975 = est_CKT_NP + 1.96 * asympt_se_np,
    est_CKT_reg,
    asympt_se_reg,
    est_CKT_reg_q025 = est_CKT_reg - 1.96 * asympt_se_reg,
    est_CKT_reg_q975 = est_CKT_reg + 1.96 * asympt_se_reg
  )

  return (invisible(df))
}


# Compute the elements for Wald-type test statistic related to
# the hypothesis $\beta = 0$ against $\beta != 0$
# vectorZ is the vector of Z in the database of length n
# vectorZToEstimate is the vector of z'_i
# vector_hat_CKT_NP is the vector of the nonparametric estimators of the CKT at every z'_i
# vector_hat_beta is the vector of coefficients, of size p'
# matrixSignsPairs is a matrix of size n * n
# inputMatrix is the inputMatrix of size n' * p'
compute_for_Wn = function(
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
      (Lambda_deriv(pointZ))^2 * abs(Gn_zipr[iprime] - (vector_hat_CKT_NP[iprime])^2)
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

  return (list(Gn_zipr = Gn_zipr,
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
# We want to compute sum_{i,j,k} w_i w_j w_k g_{i,k} g_{j,k}
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


