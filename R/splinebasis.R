#' @title For each column of an input matrix, elements which are less than the median of that column are set to 0, leaving the rest of the elements unchanged
#' @description
#'    For each column of an input matrix, elements which are less than the median of that column are set to 0, leaving the rest of the elements unchanged.
#' @param X The extra covariates, a n x K matrix. No constant included.
#' @return H A n x K matrix. All elements of the X that are less than
#' the median of their corresponding columns are set to 0, leaving the rest unchanged.
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
