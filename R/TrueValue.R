#' Calculate the True LATE tau.
#' @description  Calculate the true LATE tau in Jiang et al. (2022).
#' @param dgptype A scalar. The value can be string 1, 2, or 3,
#'    respectively corresponding to the three random data generation methods in the paper (See Jiang et al. (2022)for DGP details)
#' @param vIdx A 1xR vector. The authors set vIdx=[1 2 3 4]. Every number declares the method of covariate-adaptive randomization
#' which simulates the LATE across different CAR schemes: 1-SRS; 2-WEI; 3-BCD; 4-SBR.
#' @param n Sample size.
#' @param g Number of strata. The authors set g=4 in Jiang et al. (2022).
#' @param pi Targeted assignment probability across strata.
#' @importFrom purrr map2
#'
#' @return A list containing two vectors named tau and mPort.
#'    tau is a 1xR vector which Simulated true LATE effect, mPort is a 3xR vector.
#'    The 1st row of mPort: the LATE of never takers across varies CAR schemes,
#'    the 2nd row of mPort: the LATE of compilers across varies CAR schemes,
#'    the 3rd row of mPort: the LATE of always takers across varies CAR schemes.
#'
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' \donttest{
#'  TrueValue(dgptype = 1, vIdx = c(1,2,3,4), n=100, g = 4, pi = c(0.5,0.5,0.5,0.5))
#'  TrueValue(dgptype = 2, vIdx = c(1,2,3,4), n=100, g = 4, pi = c(0.5,0.5,0.5,0.5))
#'  TrueValue(dgptype = 3, vIdx = c(1,2,3,4), n=100, g = 4, pi = c(0.5,0.5,0.5,0.5))
#'  }

TrueValue <- function(dgptype, vIdx, n, g, pi) {
  R <- length(vIdx)
  Nsim <- 1000
  Tau <- NaN*ones(Nsim,R)
  mCom <- NaN*ones(Nsim,R)
  mAT <- NaN*ones(Nsim,R)
  mNT <- NaN*ones(Nsim,R)

  TrueValue_component <- function(sim, R) {
    for (r in 1:R) {
      DGP <- FuncDGP(dgptype = dgptype, rndflag = vIdx[r], n = n, g = g, pi = pi)
      Y1 <- DGP[["Y1"]]
      Y0 <- DGP[["Y0"]]
      D1 <- DGP[["D1"]]
      D0 <- DGP[["D0"]]

      vTmp1 <- Y1 - Y0
      Tau[sim, r] <<- mean(vTmp1[D1==1 & D0==0])
      vTmp2 <- ones(n,1)
      mAT[sim,r] <<- length(vTmp2[D1==1 & D0==1])/n
      mNT[sim,r] <<- length(vTmp2[D1==0 & D0==0])/n
      mCom[sim,r] <<- length(vTmp2[D1==1 & D0==0])/n
    }
    print(str_c('currently at ' ,sim, 'th sample for simulating true value'))
  }

  map2(.x = 1:Nsim, .y = R, .f = TrueValue_component)

  tau <- colMeans(Tau)
  mPort <- rbind(colMeans(mNT), colMeans(mCom,1), colMeans(mAT,1))

  result_list <- list(tau = tau, mPort = mPort)
  return(result_list)

}




