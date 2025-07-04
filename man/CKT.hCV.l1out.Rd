% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimationCKT.kernelCV.R
\name{CKT.hCV.l1out}
\alias{CKT.hCV.l1out}
\alias{CKT.hCV.Kfolds}
\title{Choose the bandwidth for kernel estimation of
conditional Kendall's tau using cross-validation}
\usage{
CKT.hCV.l1out(
  X1 = NULL,
  X2 = NULL,
  Z = NULL,
  range_h,
  matrixSignsPairs = NULL,
  nPairs = NULL,
  typeEstCKT = "wdm",
  kernel.name = "Epa",
  progressBar = TRUE,
  verbose = FALSE,
  observedX1 = NULL,
  observedX2 = NULL,
  observedZ = NULL
)

CKT.hCV.Kfolds(
  X1,
  X2,
  Z,
  ZToEstimate,
  range_h,
  matrixSignsPairs = NULL,
  typeEstCKT = "wdm",
  kernel.name = "Epa",
  Kfolds = 5,
  progressBar = TRUE,
  verbose = FALSE,
  observedX1 = NULL,
  observedX2 = NULL,
  observedZ = NULL
)
}
\arguments{
\item{X1}{a vector of \code{n} observations of the first variable}

\item{X2}{a vector of \code{n} observations of the second variable}

\item{Z}{vector of observed values of Z.
If Z is multivariate, then this is a matrix whose rows correspond
to the observations of Z}

\item{range_h}{vector containing possible values for the bandwidth.}

\item{matrixSignsPairs}{square matrix of signs of all pairs,
produced by \code{\link{computeMatrixSignPairs}(observedX1, observedX2)}.
Only needed if \code{typeEstCKT} is not the default 'wdm'.}

\item{nPairs}{number of pairs used in the cross-validation criteria.
Use \code{nPairs = "all"} to choose all pairs.
The default is \code{nPairs = 10 * n}, where \code{n} is the sample size.}

\item{typeEstCKT}{type of estimation of the conditional Kendall's tau.}

\item{kernel.name}{name of the kernel used for smoothing.
Possible choices are \code{"Gaussian"} (Gaussian kernel)
and \code{"Epa"} (Epanechnikov kernel).}

\item{progressBar}{if \code{TRUE}, a progressbar for each h is displayed
to show the progress of the computation.}

\item{verbose}{if \code{TRUE}, print the score of each h during the procedure.}

\item{observedX1, observedX2, observedZ}{old parameter names for \code{X1},
\code{X2}, \code{Z}. Support for this will be removed at a later version.}

\item{ZToEstimate}{vector of fixed conditioning values at which
the difference between the two conditional Kendall's tau should be computed.
Can also be a matrix whose lines are the conditioning vectors at which
the difference between the two conditional Kendall's tau should be computed.}

\item{Kfolds}{number of subsamples used.}
}
\value{
Both functions return a list with two components:
\itemize{
    \item \code{hCV}: the chosen bandwidth
    \item \code{scores}: vector of the same length as range_h giving the
    value of the CV criteria for each of the h tested.
    Lower score indicates a better fit.
}
}
\description{
Let \eqn{X_1} and \eqn{X_2} be two random variables.
The goal here is to estimate the conditional Kendall's tau
(a dependence measure) between \eqn{X_1} and \eqn{X_2} given \eqn{Z=z}
for a conditioning variable \eqn{Z}.
Conditional Kendall's tau between \eqn{X_1} and \eqn{X_2} given \eqn{Z=z}
is defined as:
\deqn{P( (X_{1,1} - X_{2,1})(X_{1,2} - X_{2,2}) > 0 | Z_1 = Z_2 = z)}
\deqn{- P( (X_{1,1} - X_{2,1})(X_{1,2} - X_{2,2}) < 0 | Z_1 = Z_2 = z),}
where \eqn{(X_{1,1}, X_{1,2}, Z_1)} and \eqn{(X_{2,1}, X_{2,2}, Z_2)}
are two independent and identically distributed copies of \eqn{(X_1, X_2, Z)}.
For this, a kernel-based estimator is used, as described in
(Derumigny & Fermanian (2019)).
These functions aims at finding the best bandwidth \code{h} among a given
\code{range_h} by cross-validation. They use either:
\itemize{
   \item \strong{leave-one-out} cross-validation:
   function \code{CKT.hCV.l1out}

   \item or \strong{K-folds} cross-validation:
   function \code{CKT.hCV.Kfolds}
}
}
\examples{
# We simulate from a conditional copula
set.seed(1)
N = 200
Z = rnorm(n = N, mean = 5, sd = 2)
conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
simCopula = VineCopula::BiCopSim(N=N , family = 1,
    par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
X1 = qnorm(simCopula[,1])
X2 = qnorm(simCopula[,2])

newZ = seq(2,10,by = 0.1)
range_h = 3:10

resultCV <- CKT.hCV.l1out(X1 = X1, X2 = X2, Z = Z,
                          range_h = range_h, nPairs = 100)

resultCV <- CKT.hCV.Kfolds(X1 = X1, X2 = X2, Z = Z,
                           range_h = range_h, ZToEstimate = newZ)

plot(range_h, resultCV$scores, type = "b")

}
\references{
Derumigny, A., & Fermanian, J. D. (2019).
On kernel-based estimation of conditional Kendall’s tau:
finite-distance bounds and asymptotic behavior.
Dependence Modeling, 7(1), 292-321.
Page 296, Equation (4).
\doi{10.1515/demo-2019-0016}
}
\seealso{
\code{\link{CKT.kernel}} for the corresponding
estimator of conditional Kendall's tau by kernel smoothing.
}
