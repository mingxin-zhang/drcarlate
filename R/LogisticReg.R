#' @title  Logistic Regression Function
#'
#' @description  Regression-adjusted Estimation of Quantile Treatment Effects under
#'   Covariate-Adaptive Randomizations.
#'
#' @usage LogisticReg(x)
#'
#' @param x A nx1 matrix.
#'
#' @return y A nx1 matrix.
#'   y equals to exp(x)/(1+exp(x)) if y is not NA and 0 else.
#'
#' @export
#'
#' @examples
#' x <- pracma::rand(5,1)
#' y <- LogisticReg(x = x)
#'
LogisticReg <- function(x) {
  y <- exp(x)/(1+exp(x))

  y[is.na(y)] <- 0 # If y=NaN, then set y=0. This corresponds to no adjustment.

  return(y)
}
