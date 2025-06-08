

#' Get Gaussian quadrature nodes and weights for integration on a specific interval
#'
#' @param nGrid number of points in the grid
#' @param kind of Gaussian quadrature method, see \code{?statmod::gauss.quad}
#' @param center center of the interval
#' @param half_length half-length of the interval
#'
#' @returns a list with components \code{nodes} and \code{weights}
#' for integration on the interval \code{[center +- half_length]}.
#'
#' @examples
#' # Interval [0, 1] = [1/2 +- 1/2]
#' grid = get.gauss.quad(nGrid = 10, kind = "legendre",
#'                       center = 1/2, half_length = 1/2)
#'
#' # Integral of \eqn{x \mapsto x} on [0, 1]
#' stopifnot(isTRUE(all.equal(sum(grid$weights) , 1) ) )
#'
#' # Integral of \eqn{x \mapsto x} on [1/4, 3/4] = [1/2 +- 1/4]
#' grid = get.gauss.quad(nGrid = 10, kind = "legendre",
#'                       center = 1/2, half_length = 1/4)
#' stopifnot(isTRUE(all.equal(sum(grid$weights) , 1/2) ) )
#'
#' @noRd
get.gauss.quad <- function(nGrid = nGrid, kind = "legendre",
                           center = 1/2, half_length)
{
  # Original grid on the interval [-1 , 1]
  grid = statmod::gauss.quad(n = nGrid, kind = kind)

  # Change of range to be on [center +- half_length]
  grid$nodes <- grid$nodes * half_length + center

  # Adapting the weights accordingly
  grid$weights <- grid$weights * half_length

  return (grid)
}



