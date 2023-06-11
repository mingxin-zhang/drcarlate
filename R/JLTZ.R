#' Reproduce the results of the Jiang et al. (2022)
#' @description  Helps the user reproduce the results of the data simulation section of Jiang et al. (2022).
#' @param iMonte A scalar. Monte Carlo sizes.
#' @param dgptype A scalar. The value can be string 1, 2, or 3,
#'    respectively corresponding to the three random data generation methods in the paper (See Jiang et al. (2022) for DGP details).
#' @param n Sample size.
#' @param g Number of strata. We set g=4 in Jiang et al. (2022).
#' @param pi Targeted assignment probability across strata.
#' @param iPert A scalar. iPert = 0 means size. Otherwise means power: iPert is the perturbation of false null.
#' @param iq A scalar. Size of hypothesis testing. The authors set iq = 0.05.
#' @param iridge A scalar. The penalization parameter in ridge regression.
#' @param seed A scalar. The random seed, the authors set seed = 1 in Jiang et al. (2022).

#' @return A table summarizing the estimated results, mProd.
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' # size, iPert = 0
#' ATEJLTZ(iMonte = 10, dgptype = 1, n = 20, g = 4,
#'     pi = c(0.5, 0.5, 0.5, 0.5), iPert = 0, iq = 0.05, iridge = 0.001, seed = 1)
#'
#' # power, iPert = 1
#' JLTZ(iMonte = 10, dgptype = 1, n = 20, g = 4,
#'     pi = c(0.5, 0.5, 0.5, 0.5), iPert = 1, iq = 0.05, iridge = 0.001, seed = 1)

JLTZ <- function(iMonte, dgptype, n, g, pi, iPert, iq = 0.05, iridge = 0.001, seed = 1) {

# Start
tic()

# Set random seed
set.seed(seed)

# test the length of pi
if (g != length(pi)){
  stop("g not equal to length(pi)!")
}

if (seed == 1) {
  # Simulate the LATE tau when set.seed(1)
  if (dgptype == 1) {
    vtau <- c(0.9203193, 0.9211907, 0.9224011, 0.9216388) # DGP 1
  } else if (dgptype == 2) {
    vtau <- c(0.9220499, 0.9216538, 0.9214587, 0.9214584) # DGP 2
  } else if (dgptype == 3) {
    vtau <- c(0.9225413,  0.9215747, 0.9215092, 0.9221616) #DGP 3
  }
} else {
  truevalue_result <- TrueValue(dgptype = dgptype, vIdx = 1:4, n = 10000, g = g, pi = pi)
  tau <- truevalue_result[["tau"]]
  print("Finish simulating true tau.")

}

# Monte Carlo Experiment
# If a value exists, cols are (1)NA (2)LP (3)LG (4)F (5)NP (6)lasso (7)2SLS (8)R
# 1-SRS
mtauhat_d1 <- NaN*ones(iMonte, 12)
msighat_d1 <- NaN*ones(iMonte, 12)
mstat_d1   <- NaN*ones(iMonte, 12)
mdeci_d1   <- NaN*ones(iMonte, 12)

# 2-WEI
mtauhat_d2 <- NaN*ones(iMonte, 12)
msighat_d2 <- NaN*ones(iMonte, 12)
mstat_d2   <- NaN*ones(iMonte, 12)
mdeci_d2   <- NaN*ones(iMonte, 12)

# 3-BCD
mtauhat_d3 <- NaN*ones(iMonte, 12)
msighat_d3 <- NaN*ones(iMonte, 12)
mstat_d3   <- NaN*ones(iMonte, 12)
mdeci_d3   <- NaN*ones(iMonte, 12)

# 4-SBR
mtauhat_d4 <- NaN*ones(iMonte, 12)
msighat_d4 <- NaN*ones(iMonte, 12)
mstat_d4   <- NaN*ones(iMonte, 12)
mdeci_d4   <- NaN*ones(iMonte, 12)


for (i in 1:iMonte) {
  print("1-SRS")
  # 1-SRS
  result_srs <- tryCatch({Output(ii = i, tau = vtau[1], dgptype = dgptype, rndflag = 1,
                  n = n, g = g, pi = pi, iPert = iPert, iq = iq, iridge = iridge)},
           error = function(e) {list(vtauhat = NA, vsighat = NA,
                                                    vstat = NA, vdeci = NA)})

  vtauhat_d1 <- result_srs[["vtauhat"]]
  vsighat_d1 <- result_srs[["vsighat"]]
  vstat_d1   <- result_srs[["vstat"]]
  vdeci_d1   <- result_srs[["vdeci"]]

  mtauhat_d1[i,] <- vtauhat_d1
  msighat_d1[i,] <- vsighat_d1
  mstat_d1[i,]   <- vstat_d1
  mdeci_d1[i,]   <- vdeci_d1

  # 2-WEI
  print("2-WEI")
  result_wei <- tryCatch({Output(ii = i, tau = vtau[2], dgptype = dgptype, rndflag = 2,
                  n = n, g = g, pi = pi, iPert = iPert, iq = iq, iridge = iridge)},
           error = function(e) {list(vtauhat = NA, vsighat = NA,
                                            vstat = NA, vdeci = NA)})
  vtauhat_d2 <- result_wei[["vtauhat"]]
  vsighat_d2 <- result_wei[["vsighat"]]
  vstat_d2   <- result_wei[["vstat"]]
  vdeci_d2   <- result_wei[["vdeci"]]

  mtauhat_d2[i,] <- vtauhat_d2
  msighat_d2[i,] <- vsighat_d2
  mstat_d2[i,]   <- vstat_d2
  mdeci_d2[i,]   <- vdeci_d2

  # 3-BCD
  print("3-BCD")
  result_bcd <- tryCatch({Output(ii = i, tau = vtau[3], dgptype = dgptype, rndflag = 3,
                  n = n, g = g, pi = pi, iPert = iPert, iq = iq, iridge = iridge)},
           error = function(e) {list(vtauhat = NA, vsighat = NA,
                                                    vstat = NA, vdeci = NA)})


  vtauhat_d3 <- result_bcd[["vtauhat"]]
  vsighat_d3 <- result_bcd[["vsighat"]]
  vstat_d3   <- result_bcd[["vstat"]]
  vdeci_d3   <- result_bcd[["vdeci"]]

  mtauhat_d3[i,] <- vtauhat_d3
  msighat_d3[i,] <- vsighat_d3
  mstat_d3[i,]   <- vstat_d3
  mdeci_d3[i,]   <- vdeci_d3

  # 4-SBR
  print("4-SBR")
  result_sbr <- tryCatch({Output(ii = i, tau = vtau[4], dgptype = dgptype, rndflag = 4,
                  n = n, g = g, pi = pi, iPert = iPert, iq = iq, iridge = iridge)},
           error = function(e) {list(vtauhat = NA, vsighat = NA,
                                                    vstat = NA, vdeci = NA)})

  vtauhat_d4 <- result_srs[["vtauhat"]]
  vsighat_d4 <- result_srs[["vsighat"]]
  vstat_d4   <- result_srs[["vstat"]]
  vdeci_d4   <- result_srs[["vdeci"]]

  mtauhat_d4[i,] <- vtauhat_d4
  msighat_d4[i,] <- vsighat_d4
  mstat_d4[i,]   <- vstat_d4
  mdeci_d4[i,]   <- vdeci_d4
}

# evaluation
vProb_d1 <- colMeans(mdeci_d1, na.rm = TRUE) # 1x12. size or power depending on the context.
vProb_d2 <- colMeans(mdeci_d2, na.rm = TRUE)
vProb_d3 <- colMeans(mdeci_d3, na.rm = TRUE)
vProb_d4 <- colMeans(mdeci_d4, na.rm = TRUE)
mProd <- t(rbind(vProb_d1,vProb_d2,vProb_d3,vProb_d4))

return(mProd)

toc()
}



