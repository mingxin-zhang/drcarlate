#' @title Compute Estimated Tau
#' @description `tau` computes the estimated tau in Jiang et al.(2022); scalar.
#' @param muY1 A nx1 vector of hat\{mu\}^Y(A=1)s
#' @param muY0 A nx1 vector of hat\{mu\}^Y(A=0)s
#' @param muD1 A nx1 vector of hat\{mu\}^D(A=1)s
#' @param muD0 A nx1 vector of hat\{mu\}^D(A=0)s
#' @param A A nx1 vector.
#' @param S A nx1 vector.
#' @param Y A nx1 vector.
#' @param D A nx1 vector.
#'
#' @return Tauhat scalar; LATE estimate.
#'
#' @export
#' @references Jiang L, Linton O B, Tang H, et al. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.

#' @examples
#' DGP <- FuncDGP(dgptype = 1, rndflag = 1, n = 200, g = 4, pi = 0.5)
#' muY1 <- DGP[["Y1"]]
#' muY0 <- DGP[["Y0"]]
#' muD1 <- DGP[["D1"]]
#' muD0 <- DGP[["D0"]]
#' A <- DGP[["A"]]
#' S <- DGP[["S"]]
#' Y <- DGP[["Y"]]
#' D <- DGP[["D"]]
#' tau(muY1, muY0, muD1, muD0, A, S, Y, D)
#'
#'
tau <- function(muY1, muY0, muD1, muD0, A, S, Y, D) {
  vPihat <-  pihat(A = A, S = S)
  iDe <- mean(A*(D-muD1)/vPihat-(1-A)*(D-muD0)/(1-vPihat)+muD1-muD0)
  iNu <- mean(A*(Y-muY1)/vPihat-(1-A)*(Y-muY0)/(1-vPihat)+muY1-muY0)
  tauhat <- iNu/iDe
  if (is.nan(tauhat)) {
#    print("iNu:")
#    print(iNu)
#    print("iDe:")
#    print(iDe)
    stop("Error: tauhat==NaN")
#    return(list(muY1 = muY1,
#                muY0 = muY0,
#                muD1 = muD1,
#                muD0 = muD0,
#                iDe = iDe,
#                iNu = iNu))
  } else {
    return(tauhat)
  }
}
