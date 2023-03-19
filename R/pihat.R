#' @title Compute Estimated Pi
#' @description
#'  Pihat computes the estimated pi in Jiang et al.(2022) and stacks them in an nx1 vector.
#' @param A A nx1 vector.
#' @param S A nx1 vector.
#'
#' @return A nx1 cector, each element corresponds to  the estimated pi in Jiang et al.(2022).
#' @export
#' @references Jiang L, Linton O B, Tang H, et al. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' DGP <-FuncDGP(dgptype = 1,rndflag = 2,n = 100,g = 4,pi = 0.5)
#' A <- DGP[["A"]]
#' S <- DGP[["S"]]
#' pihat(A = A, S = S)
#'
pihat <- function(A, S) {
  vPihat  <- NaN*ones(length(A),1)

  for (s in 1:max(S)) {
    ns <-  sum(S == s)
    n1s <- sum(S == s & A == 1)
    pi_s <- n1s/ns
    vPihat[S == s] <- pi_s
  }

  return(vPihat)
}
