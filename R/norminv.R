#' Inverse of the normal cumulative distribution function (cdf)
#'
#'Returns the inverse cdf for the normal distribution with mean MU and standard deviation SIGMA at P value
#' Reference: https://rdrr.io/github/maxto/qapi/src/R/stats.R
#' @param p probability value in range 0-1
#' @param mu mean value
#' @param sigma standard deviation
#'
#' @return numeric
#' @export
#'
#' @examples
#' xx <- c(0.003,0.026,0.015,-0.009,-0.014,-0.024,0.015,0.066,-0.014,0.039)
#' norminv(0.01,mean(xx),sd(xx))
norminv <- function(p,mu=0,sigma=1) {
  x0 <- -sqrt(2)*pracma::erfcinv(2*p)
  x0*sigma+mu
}
