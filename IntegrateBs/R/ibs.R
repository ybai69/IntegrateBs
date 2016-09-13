#' Integration for B-splines
#'
#' Calculate the integral of a B-spline function. 
#' @param x Numerical value or vector. The value(s) at which to evaluate the integral of the B-spline; must be in the bewteen min(knots) and max(knots).
#' @param knots A numeric vector of knot positions.
#' @param ord An integer >=1. The order of the B-spline integrand function to be integrated. Equals degree plus 1.
#' @param coef A numerical vector. The coefficients (de Boor points) defining the B-spline integrand function.
#' @return A numerical equal to the integral(s).
#' @details The function returns the integral(s) of the B-spline function specified by knots knots, order ord, and coefficients coef, from the minimum knot position to each x value. The evaluation is based on a closed form expression of the integral in terms of higher order B-splines, given on page 128 of de Boor (2001).
#' @export
#' @references de Boor, C (2001)  \emph{A Practical Guide to Splines}. New York: Springer.
#' @examples
#' library(splines)
#' f <- function(x) x + 2 * x^2 - 3 * x^3 
#' n <- 200
#' set.seed(123)
#' x <- runif(n)
#' y <- f(x) + rnorm(n, sd = 0.1)
#' kns <- c(rep(0, 4), 1:4 * 0.2, rep(1, 4))
#' bs.c <- splineDesign(kns, x, 4)
#' coeff <- as.matrix(lm(y ~ bs.c-1)$coefficients)
#' f.b <- function(x, coeff) splineDesign(kns, x, 4) %*% coeff
#' integrate(f.b, 0, 1, coeff)
#' ibs(1,kns,4,coeff)
#' integrate(f, 0, 1)
#' plot(x,y)
#' curve(f(x), add = TRUE)
#' points(x,fitted(lm(y~bs.c-1)),col="blue",lty=1)
ibs<-function(x,knots=NULL,ord = 4,coef=rep(1,length(knots)-ord)){
  if (requireNamespace("splines", quietly = TRUE)) {
    if(ord<1)stop("ord should be greater or equal to 1!")
    if(ord!=as.integer(ord))stop("ord should be integer!") 
    if(length(coef)!=length(knots)-ord)stop("length(knots)-ord!=length(coef)!")
    nk<-length(knots)
    knots<-c(knots,knots[nk])
    int<-0
    for(i in 1:(nk-ord)){
      s<-0
      for(j in 1:i){
        s<-s+coef[j]*(knots[j+ord]-knots[j])/ord
      }
    int<-int+s*splines::splineDesign(knots,x,ord+1)[i]
    }
  }
  return(int)
}
