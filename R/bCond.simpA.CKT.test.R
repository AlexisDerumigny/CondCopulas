

##
##  WE USE HERE hat.tau.1 instead of ordinary hat.tau
##

# Computes asymptotic covariance matrix from Theorem 1
# for a three dimensional random sample
# Notations from the version "October 2019"
#
# Input:
#   x:          Numeric matrix or data frame on ORIGINAL scale
#               with two columns.
#   ind:        Matrix of size n * m where n is the sample size and m the index of the box
#               where w[i,j] = 1(X_{i,J} \in A_{j,J})
#
#   hat.tau:    vector of size N.boxes of estimated CKT
#   tau.cond.true : true unconditional CKT under H0
#
#
# Output:
#   List containing the value of the test statistic.
#   $statistic: Numeric vector of length one giving the value of the test statistic.
#   $cov:       Estimates asymptotic covariance matrix
fun.mult.norm.prop.2 <- function(x, ind, hat.tau, tau.cond.true)
{

  # Denote sample size and number of threshold points
  n        <- nrow(x)
  N.boxes  <- ncol(ind)


  ##############################################
  # Determine weights "w" from the indicators "ind"

  w <- ind

  for(k in 1:N.boxes)
  {
    # Weights "w"
    w[,k] <- w[,k]/sum(w[,k])
  }


  ##############################################
  # Compute hat.D.k and g_0^* the v_i
  hat.D.k <-  array(0, N.boxes)
  g.star <-  array(0, c(n, N.boxes))
  bar.v.i <-  array(0, c(n, N.boxes))

  for(k in 1:N.boxes)
  {
    for(i in 1:n)
    {
      # Compute hat.D.k

      # This part was vectorized
      ARG        <- ind[i,k] * ind[-i,k] * (x[i,1]<x[-i,1]) * (x[i,2]<x[-i,2])
      hat.D.k[k] <- hat.D.k[k] + sum(ARG)


      # Compute g_0^*

      temp.integral <- 0
      # This part was vectorized
      pi.k           <- ( (x[i,1]<x[ ,1]) * (x[i,2]<x[ ,2]) + (x[i,1]>x[ ,1]) * (x[i,2]>x[ ,2]) )/2
      temp.integral  <- ind[ ,k]*pi.k

      temp  <- sum(temp.integral)/n
      g.star[i,k]   <- ind[i,k]*temp
    }
    bar.v.i[ ,k] <- 2*(g.star[ ,k] - hat.D.k[k])
  }
  hat.D.k <-  hat.D.k/(n*(n-1))


  ##############################################
  # Estimated probs
  est.p.k <- apply(ind[ , 1:N.boxes], 2, mean)


  # covariance of bar.v.i
  COVAR1 <- stats::cov(bar.v.i)
  TERM1  <- diag( (diag(COVAR1) + 4*hat.D.k^2)/(est.p.k^4))
  TERM2  <- diag( (1+hat.tau)^2/(4*est.p.k) )
  COVAR  <- 16*(TERM1 - TERM2)

  # stat   <- sqrt(n) * (hat.tau - tau.cond.true)/sqrt(diag(COVAR))

  out = list(
    # "stat" = stat ,
    "COVAR" = COVAR ,
    "hat.D.k" = hat.D.k,
    "hat.tau" = hat.tau,
    "bar.v.i" = bar.v.i,
    "est.p.k" =  est.p.k)
  return(out)
}


# Computes asymptotic covariance matrix from Proposition 4
# for a d dimensional random sample
# Notations from the version "October 2019"
# This function uses the R-function "fun.mult.norm.prop.2" for d=2
#
# Input:
#   xI:         Numeric matrix or data frame on ORIGINAL scale
#               with d columns. Last vector is conditioning variable
#
#   ind:        Matrix of indicators of size n * p(p-1) / 2
#
#   hat_CKT:    Matrix of estimated conditional Kendall's taus
#               of size p(p-1) / 2 * N.boxes
#               such as given by the function Get(Traverse(estimatedTree, filterFun = isLeaf), "CKT")
#
#   tau.cond.true: True conditional taus as vector of length p(p-1)/2
#                  First row, second row etc of the upper diagonal matrix
#
#
# Output:
#   List containing the value of the test statistic.
#   $statistic: Numeric vector of length one giving the value of the test statistic.
#   $cov:       Estimates asymptotic covariance matrix
fun.mult.norm.prop.4 <- function(xI, ind,  hat_CKT, tau.cond.true)
{

  # Denote sample size and number of threshold points
  n        <- nrow(xI)
  p        <- ncol(xI)

  N.boxes  <- ncol(ind)  # in the paper, it is $m$

  tilde.v   <- c()
  hat.D.k   <- c()

  # Test statistics
  dd <- p*(p-1)/2
  TT <- numeric(dd)

  # All possible pairs
  iPair = 1
  for(i in 1:(p-1))
  {
    for(j in (i+1):p )
    {
      # Two dimenional data
      xxx <- cbind(xI[, i], xI[, j])

      # The test statistic
      out.temp <- fun.mult.norm.prop.2(x = xxx,  ind = ind, hat.tau = hat_CKT[iPair, ],
                                       tau.cond.true[iPair])

      # Dimension of tilde.v = N.boxes * p*(p-1)/2
      tilde.v   <- cbind(tilde.v , out.temp$bar.v.i)

      hat.D.k   <- c(hat.D.k , out.temp$hat.D.k)

      # Updating of the pair counter
      iPair = iPair + 1
    }
  }
  est.p.k   <- out.temp$est.p.k

  # Covariance matrix
  DELTA <- array(0, c(N.boxes*dd, N.boxes*dd))

  for(i in 1:dd)
  {
    for(j in 1:dd)
    {
      index.row    <- c( ((i-1)*N.boxes+1) : (i*N.boxes) )
      index.column <- c( ((j-1)*N.boxes+1) : (j*N.boxes) )

      # covariance of bar.v.i
      COVAR1 <- stats::cov(tilde.v[ ,index.row], tilde.v[ ,index.column])
      TERM1  <- diag( (diag(COVAR1) +
                         4 * hat.D.k[index.row] * hat.D.k[index.column])/(est.p.k^4))
      TERM2  <- diag( (1+hat_CKT[i,]) * (1+hat_CKT[j,]) / (4*est.p.k) )
      COVAR  <- 16*(TERM1 - TERM2)

      DELTA[index.row, index.column] <- COVAR
    }
  }

  round(DELTA, digits=3)

  # stat   <- sqrt(n) * (as.numeric(hat_CKT) - tau.cond.true)/sqrt(diag(DELTA))

  out = list(
    # "stat" = stat ,
    "DELTA" = DELTA ,
    "hat.D.k" = hat.D.k,
    "hat_CKT" = hat_CKT,
    "tilde.v"= tilde.v,
    "n"= n,
    "p"= p,
    "N.boxes" = N.boxes)
  return(out)
}


fun.Maximum.Stat.gen.dim <- function(hat_CKT, tau.cond.true) {

  # Computes Maximum Statistic based on Proposition 4
  # Notations from the version "October 2019"
  # Contrast Matrix is fixed (see paper)
  #
  # Input:
  #   R.Object: Output
  #             of the R-function "fun.mult.norm.prop.4"
  #
  #
  # Output:
  #   List containing the value of the test statistic.
  #   $statistic: Numeric vector of length one giving the value of the test statistic.


  # Contrast Matrix
  dd <- nrow(hat_CKT)
  N.boxes <- ncol(hat_CKT)

  TT   <- array(0, c( dd*(N.boxes-1) , dd*N.boxes) )

  temp <- -diag(N.boxes-1)
  TT.temp   <- cbind(rep(1, N.boxes-1), temp)

  for(i in 1:dd)
  {
    index.row    <- c( ((i-1)*(N.boxes-1)+1) : (i*(N.boxes-1)) )
    index.column <- c( ((i-1)*N.boxes+1) : (i*N.boxes) )

    TT[index.row,index.column] <- TT.temp
  }

  stat <- max( abs( TT %*% (as.numeric(hat_CKT) - as.numeric(tau.cond.true)) ) )

  return(stat)
}



# Computes Wald Statistic based on Proposition 4
# for a d-dimensional
# Notations from the version "October 2019"
# Contrast Matrix is fixed (see paper)
#
# Input:
#   R.Object: Output
#             of the R-function "fun.mult.norm.prop.4"
#
#
# Output:
#   List containing the value of the test statistic.
#   $statistic: Numeric vector of length one giving the value of the test statistic.
#   $p.value: Numeric vector of length one giving the p-value.
fun.Wald.Stat.gen.dim <- function(R.Object)
{
  # Finding dimension values
  p  <- R.Object$p
  n <- R.Object$n
  N.boxes <- R.Object$N.boxes

  # Contrast Matrix
  dd <- (p-1)*p/2

  TT   <- array(0, c( dd*(N.boxes-1) , dd*N.boxes) )


  # The first "contrast" matrix. The first column will be compared with the other ones
  temp      <- -diag(N.boxes-1)
  TT.temp   <- cbind(rep(1, N.boxes-1), temp)
  dim.temp  <- dim(TT.temp)[2] # number of columns

  # Permute columns randomly. The randomly chosen "main box" will be compared with other ones
  # It is from the conversation from November 4th 2019
  ind.temp <- sample(1:dim.temp, replace = FALSE)
  TT.temp  <- TT.temp[, ind.temp]
  # rm(ind.temp, dim.temp)

  for(i in 1:dd)
  {
    index.row    <- c( ((i-1)*(N.boxes-1)+1) : (i*(N.boxes-1)) )
    index.column <- c( ((i-1)*N.boxes+1) : (i*N.boxes) )

    TT[index.row,index.column] <- TT.temp
  }


  WW    <- as.numeric(R.Object$hat_CKT)
  Delta <- R.Object$DELTA

  covMatrix = n * t(TT) %*% solve(TT %*% Delta %*% t(TT)) %*% TT

  stat <- t(WW) %*% covMatrix %*%WW

  # In practice this matrix may not be positive definite
  # so we take the absolute value in these extreme cases
  p.value <- 1 - stats::pchisq(q = as.numeric(abs(stat)), df = (N.boxes)*dd)

  out = list("stat" = stat, "p.value" = p.value, "covMatrix" = covMatrix)
  return(out)
}



