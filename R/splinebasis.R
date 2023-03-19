#' @title Spline The Matrix
#' @description
#'    Splinebasis is based on Yubo's code.
#'    It generates polynomial spline basis in Chen (2007) p 5571.
#' @param X The extra covariate, a nxK matrix. No constant included.
#' @return H A nxK matrix whose each element is bigger than its column median.
#' @export
#'
#' @examples
#' library(pracma)
#' X <- rand(4,4)
#' H <- splinebasis(X = X)
#'
splinebasis <- function(X) {
  n <- size(X,1)
  K <- size(X,2)

  H <- zeros(n, K)

  for (k in 1:K) {
    vTmp1 = X[,k] - median(X[,k])
    vTmp2 = vTmp1 * (vTmp1 > 0)

    H[,k] <- vTmp2
  }
  return(H)
}
