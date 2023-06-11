#' Calculates the true ATE effect.
#' @description
#' ATETrueValue is the version of TrueValue under full compliance.
#' @param dgptype A scalar. The value can be string 1, 2, or 3,
#'    respectively corresponding to the three DGP schemes in the paper (See Jiang et al. (2022) for DGP details).
#' @param vIdx A 1xR vector. The authors set vIdx=[1 2 3 4] in Jiang et al. (2022). Every number declares the method of covariate-adaptive randomization.
#'    1-SRS; 2-WEI; 3-BCD; 4-SBR.
#' @param n Sample size.
#' @param g Number of strata. The authors set g=4 in Jiang et al. (2022).
#' @param pi Targeted assignment probability across strata.
#'
#' @importFrom purrr map2
#' @return A 1xR vector. Simulated true ATE effect.
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' \donttest{
#'  ATETrueValue(dgptype = 1, vIdx = c(1,2,3,4), n = 100, g = 4, pi = c(0.5,0.5,0.5,0.5))
#'  ATETrueValue(dgptype = 2, vIdx = c(1,2,3,4), n = 100, g = 4, pi = c(0.5,0.5,0.5,0.5))
#'  ATETrueValue(dgptype = 3, vIdx = c(1,2,3,4), n = 100, g = 4, pi = c(0.5,0.5,0.5,0.5))
#'  }


ATETrueValue <- function(dgptype, vIdx, n, g, pi) {
  R <- length(vIdx)
  Nsim <- 1000
  Tau <- NaN*ones(Nsim,R)

  ATETrueValue_component <- function(sim, R) {
    for (r in 1:R) {
      DGP <- ATEDGP(dgptype = dgptype, rndflag = vIdx[r], n = n, g = g, pi = pi)
      Y1 <- DGP[["Y1"]]
      Y0 <- DGP[["Y0"]]

      Tau[sim, r] <<- mean(Y1 - Y0)

          }
    print(str_c('currently at ' ,sim, 'th sample for simulating true value'))
  }

  map2(.x = 1:Nsim, .y = R, .f = ATETrueValue_component)

  tau <- colMeans(Tau)
  return(tau)
}
