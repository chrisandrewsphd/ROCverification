#' Compute area by trapezoid rule
#'
#' @param x horizontal coordinates
#' @param y vertical coordinates
#'
#' @return area under curve
#' @export
#'
#' @examples trapezoid(seq(0, 10), rnorm(11, 5, 1))
trapezoid <- function(x, y) {
  if ((n <- length(x)) != length(y)) stop("x and y must have equal length")
  y <- y[order(x)]
  x <- x[order(x)]
  
  midy <- (y[-1] + y[-n])/2
  deltax <- diff(x)
  area <- sum(midy * deltax)
  
  return(area)
}
