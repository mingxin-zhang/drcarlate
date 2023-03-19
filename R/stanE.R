#' @title Compute Estimated Sigma
#' @description stanE_v1 computes hat\{sigma\} in Jiang et al.(2022).
#' @param muY1 A nx1 vector of hat\{mu\}^Y(A=1)s
#' @param muY0 A nx1 vector of hat\{mu\}^Y(A=0)s
#' @param muD1 A nx1 vector of hat\{mu\}^D(A=1)s
#' @param muD0 A nx1 vector of hat\{mu\}^D(A=0)s
#' @param A A nx1 vector.
#' @param S A nx1 vector.
#' @param Y A nx1 vector.
#' @param D A nx1 vector.
#' @param tauhat A scalar; LATE estimate
#'
#' @return The estimated sigma in Jiang et al.(2022); scalar.
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
#' tauhat <- tau(muY1, muY0, muD1, muD0, A, S, Y, D)
#' stanE(muY1, muY0, muD1, muD0, A, S, Y, D, tauhat)
#'
stanE <- function(muY1, muY0, muD1, muD0, A, S, Y, D, tauhat) {
  n <-  length(S)
  vPihat <-  pihat(A = A, S = S)
  vXi_til_1 <-  ((1-1/vPihat)*muY1 - muY0 + Y/vPihat) -
    tauhat*((1-1/vPihat)*muD1 - muD0 + D/vPihat)
  vXi_til_0 <- ((1/(1-vPihat) - 1)*muY0 + muY1 - Y/(1-vPihat)) -
    tauhat*((1/(1-vPihat)-1)*muD0 + muD1 - D/(1-vPihat))
  vXi_2 <- NaN*ones(n,1)
  vXi_hat_1 <- NaN*ones(n,1)
  vXi_hat_0 <- NaN*ones(n,1)

  for (s in 1:max(S)) {
    vXi_hat_1[S==s] <- vXi_til_1[S==s] - mean(vXi_til_1[S==s & A==1])
    vXi_hat_0[S==s] <- vXi_til_0[S==s] - mean(vXi_til_0[S==s & A==0])
    vXi_2[S==s] <- mean(Y[S==s & A==1] - tauhat*D[S==s & A==1]) -
      mean(Y[S==s & A==0] - tauhat*D[S==s & A==0])
  }

  iNu <- mean(A*vXi_hat_1^2 + (1-A)*vXi_hat_0^2)
  iDe <- (mean(A*(D-muD1)/vPihat - (1-A)*(D-muD0)/(1-vPihat) +
                 muD1 - muD0))^2
  sigmahat <- sqrt(iNu/iDe)

  return(sigmahat)
}

