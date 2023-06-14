#' @title Compute Estimated Treatment Assignment Probabilities
#' @description
#'  Pihat computes the targeted treatment assignment probabilities across all strata in Jiang et al. (2022) and stacks them in an nx1 vector.
#' @param A A nx1 vector.
#' @param S A nx1 vector.
#' @param stratum A vector containing the stratum identifiers; no repetition allowed. For example, if S=c(2, 3, 2, 4), then we would set stratum=c(2, 3, 4). Note that if S contains all the integers from 1 to max(S), such as S = c(4, 1, 2, 4, 5, 3), there is no need specifying the stratum parameter. Its default value is NULL.

#' @return A nx1 cector, each element corresponds to the targeted treatment assignment probabilities across all strata in Jiang et al. (2022).
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' DGP <-FuncDGP(dgptype = 1,rndflag = 2,n = 100,g = 4,pi = c(0.5, 0.5, 0.5, 0.5))
#' A <- DGP[["A"]]
#' S <- DGP[["S"]]
#' pihat(A = A, S = S)
#'
pihat <- function(A, S, stratum = NULL) {
  vPihat  <- NaN*ones(length(A),1)

  if (is.null(stratum)) {
    for (s in 1:max(S)) {
      ns <-  sum(S == s)
      n1s <- sum(S == s & A == 1)
      pi_s <- n1s/ns
      vPihat[S == s] <- pi_s
    }
  } else if (!is.null(stratum)) {
    for (j in 1:length(stratum)) {
      s <- stratum[j]
      ns <-  sum(S==s)
      nls <-  sum(S == s & A == 1)
      pi_s <- nls/ns
      vPihat[S == s] <- pi_s
    }
  }

  return(vPihat)
}
