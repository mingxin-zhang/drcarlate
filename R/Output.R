#' @title Computes All the Estimators
#' @description
#'    Output is an integrated function that computes all the estimates used in Jiang et al.(2022).
#'    See the paper for more details.
#'
#' @param ii Monte carlo index.
#' @param tau A scalar. The simulated true LATE effect.
#' @param dgptype A Scalar. 1, 2, 3 (See Jiang et al.(2022) for DGP details).
#' @param rndflag Method of CAR(covariate-adaptive randomizations).
#'    Its value can be 1, 2, 3 or4. 1-SRS; 2-WEI; 3-BCD; 4-SBR.
#'    See Jiang et al.(2022) for more details about CAR.
#' @param n Sample size.
#' @param g Number of strata. We set g=4 in Jiang et al.(2022).
#' @param pi Targeted assignment probability.
#' @param iPert A scalar. iPert =0 means size. Otherwise means power: iPert is the perturbation of false null.
#' @param iq  Size of hypothesis testing. We set iq = 0.05.
#' @param iridge  A scalar. The penalization parameter in ridge regression.
#'
#' @return A list containing four matrices named vtauhat, vsighat, vstat and vdeci respectively.
#'    vtauhat is a 1x8 vector: (1) NA (2) LP (3) LG (4) F (5) NP (6) lasso (7) 2SLS (8) NP+lasso.
#'    vsighat is a 1x8 vector: unscaled standard errors for vtauhat.
#'    vstat is a 1x8 vector: test statistic.
#'    vdeci is a 1x8 logical vector: if applicable, 1 means rejecting the null. 0 means not rejecting the null.
#'
#' @export
#' @references Jiang L, Linton O B, Tang H, et al. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.#'
#' @examples
#' Output(ii = 1, tau = 0.9122762, dgptype = 1,
#'        rndflag = 4, n = 2000, g = 4, pi = 0.5,
#'        iPert = 1, iq = 0.05, iridge = 0.001)
Output <- function(ii, tau, dgptype, rndflag, n, g, pi, iPert, iq, iridge) {

  # Generate data
  DGP <- FuncDGP(dgptype = dgptype, rndflag = rndflag, n = n, g = g, pi = pi)

  Y <- DGP[["Y"]]
  X <- DGP[["X"]]
  S <- DGP[["S"]]
  A <- DGP[["A"]]
  D <- DGP[["D"]]
  vtauhat <- NaN*ones(1,8)
  vsighat <- NaN*ones(1,8)
  vstat <- NaN*ones(1,8)
  vdeci   <- NaN*ones(1,8)

  # No adjustment (z)
  muY0_z <- zeros(n,1)
  muY1_z <- zeros(n,1)
  muD0_z <- zeros(n,1)
  muD1_z <- zeros(n,1)
  vtauhat[1] <- tau(muY1 = muY1_z, muY0 = muY0_z, muD1 = muD1_z, muD0 = muD0_z,
                    A = A, S = S, Y = Y, D = D)
  vsighat[1] <- stanE(muY1 = muY1_z, muY0 = muY0_z, muD1 = muD1_z, muD0 = muD0_z,
                      A = A, S = S, Y = Y, D = D, tauhat = vtauhat[1])
  vstat[1] <- abs(sqrt(n)*(vtauhat[1]-tau-iPert))/vsighat[1]
  vdeci[1] <- vstat[1] > norminv(1-iq/2)

  if (dgptype != 3) {

    # 2SLS
    mS_iv <- zeros(n, max(S))
    for (s in 1:max(S)) {
      mS_iv[S==s, s] <- 1
    }

    mX_iv <-cbind(D, X, mS_iv)
    mZ_iv <- cbind(A, X, mS_iv)
    vPara_iv <- mldivide(t(mZ_iv) %*% mZ_iv, t(mZ_iv) %*% Y)
    mE_iv2 <- diag(((Y - mX_iv %*% vPara_iv)^2)[,])
    m0_iv <- mldivide(t(mZ_iv) %*% mX_iv/n, t(mZ_iv) %*% mE_iv2 %*% mZ_iv/n) %*% inv(t(mX_iv) %*% mZ_iv/n)
    vtauhat[7] <- vPara_iv[1]
    vsighat[7] <- sqrt(m0_iv[1,1])
    vstat[7] <- abs(sqrt(n)*(vtauhat[7]-tau-iPert))/vsighat[7]
    vdeci[7] <- vstat[7]>norminv(1-iq/2)

    # Linear + Linear model (a)
    muY0_a <- NaN*ones(n,1)
    muY1_a <- NaN*ones(n,1)
    muD0_a <- NaN*ones(n,1)
    muD1_a <- NaN*ones(n,1)

    for (s in 1:max(S)) {
      result <- LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = s, modelflag = 1, iridge = iridge)

      theta_0s_a <- result[["theta_0s"]]
      theta_1s_a <- result[["theta_1s"]]
      beta_0s_a <- result[["beta_0s"]]
      beta_1s_a <- result[["beta_1s"]]

      muY0_a[S==s] <- X[S==s,] %*% theta_0s_a
      muY1_a[S==s] <- X[S==s,] %*% theta_1s_a
      muD0_a[S==s] <- X[S==s,] %*% beta_0s_a
      muD1_a[S==s] <- X[S==s,] %*% beta_1s_a
    }

    vtauhat[2] <- tau(muY1 = muY1_a, muY0 = muY0_a, muD1 = muD1_a, muD0 = muD0_a,
                      A = A, S = S, Y = Y, D = D)
    vsighat[2] <- stanE(muY1 = muY1_a, muY0 = muY0_a, muD1 = muD1_a, muD0 = muD0_a,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[2])
    vstat[2] <- abs(sqrt(n)*(vtauhat[2] - tau - iPert))/vsighat[2]
    vdeci[2] <- vstat[2] > norminv(1-iq/2)

    # Linear + Logistic model (b)
    muY0_b <- NaN*ones(n,1)
    muY1_b <- NaN*ones(n,1)
    muD0_b <- NaN*ones(n,1)
    muD1_b <- NaN*ones(n,1)

    for (s in 1:max(S)) {
      result <- LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = s, modelflag = 2, iridge = iridge)

      theta_0s_b <- result[["theta_0s"]]
      theta_1s_b <- result[["theta_1s"]]
      beta_0s_b <- result[["beta_0s"]]
      beta_1s_b <- result[["beta_1s"]]

      muY0_b[S==s] <- X[S==s,] %*% matrix(theta_0s_b[2:length(theta_0s_b)]) + theta_0s_b[1]
      muY1_b[S==s] <- X[S==s,] %*% matrix(theta_1s_b[2:length(theta_0s_b)]) + theta_1s_b[1]
      muD0_b[S==s] <- LogisticReg(x = (X[S==s,] %*% matrix(beta_0s_b[-1]) + beta_0s_b[1]))
      muD1_b[S==s] <- LogisticReg(x = (X[S==s,] %*% matrix(beta_1s_b[-1]) + beta_1s_b[1]))
    }

    vtauhat[3] <- tau(muY1 = muY1_b, muY0 = muY0_b, muD1 = muD1_b, muD0 = muD0_b,
                      A = A, S = S, Y = Y, D = D)
    vsighat[3] <- stanE(muY1 = muY1_b, muY0 = muY0_b, muD1 = muD1_b, muD0 = muD0_b,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[3])
    vstat[3] <- abs(sqrt(n)*(vtauhat[3] - tau - iPert))/vsighat[3]
    vdeci[3] <- vstat[3] > norminv(1-iq/2)

    # Further Efficiency Improvement (c)
    X_c <- cbind(X, muD1_b, muD0_b)
    muY0_c <- NaN*ones(n,1)
    muY1_c <- NaN*ones(n,1)
    muD0_c <- NaN*ones(n,1)
    muD1_c <- NaN*ones(n,1)

    for (s in 1:max(S)) {
      result <- LinearLogit(Y = Y, D = D, A = A, X = X_c, S = S, s = s, modelflag = 1, iridge = iridge)

      theta_0s_c <- result[["theta_0s"]]
      theta_1s_c <- result[["theta_1s"]]
      beta_0s_c <- result[["beta_0s"]]
      beta_1s_c <- result[["beta_1s"]]

      muY0_c[S==s] <- X_c[S==s,] %*% theta_0s_c
      muY1_c[S==s] <- X_c[S==s,] %*% theta_1s_c
      muD0_c[S==s] <- X_c[S==s,] %*% beta_0s_c
      muD1_c[S==s] <- X_c[S==s,] %*%beta_1s_c
    }

    vtauhat[4] <- tau(muY1 = muY1_c, muY0 = muY0_c, muD1 = muD1_c, muD0 = muD0_c,
                      A = A, S = S, Y = Y, D = D)
    vsighat[4] <- stanE(muY1 = muY1_c, muY0 = muY0_c, muD1 = muD1_c, muD0 = muD0_c,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[4])
    vstat[4] <- abs(sqrt(n)*(vtauhat[4] - tau - iPert))/vsighat[4]
    vdeci[4] <- vstat[4] > norminv(1-iq/2)

    # Nonparametric (d)
    mH <-splinebasis(X = X)
    X_d <- cbind(X, X^2,mH, mH[,1]*mH[,2], X[,1]*X[,2])
    muY0_d <- NaN*ones(n,1)
    muY1_d <- NaN*ones(n,1)
    muD0_d <- NaN*ones(n,1)
    muD1_d <- NaN*ones(n,1)

    for (s in 1:max(S)) {
      #for the case A==0
      X_d_tmp_0 <- X_d
      X_d_tmp_0_a <- X_d_tmp_0[S==s & A==0,]
      X_d_tmp_0 <- X_d_tmp_0[, colSums(X_d_tmp_0_a^2) != 0]

      result <- LinearLogit(Y = Y, D = D, A = A, X = X_d_tmp_0, S = S, s = s,
                            modelflag = 2, iridge = iridge)

      theta_0s_d <- result[["theta_0s"]]
      beta_0s_d  <- result[["beta_0s"]]

      muY0_d[S==s] <- X_d_tmp_0[S==s,] %*% matrix(theta_0s_d[-1]) + theta_0s_d[1]
      muD0_d[S==s] <- LogisticReg(X_d_tmp_0[S==s,] %*% matrix(beta_0s_d[-1]) + beta_0s_d[1])

      # for the case A==1
      X_d_tmp_1 <- X_d
      X_d_tmp_1_a <- X_d_tmp_1[S==s & A==1,]
      X_d_tmp_1 <- X_d_tmp_1[, colSums(X_d_tmp_1_a^2) != 0]

      result <- LinearLogit(Y = Y, D = D, A = A, X = X_d_tmp_1, S = S, s = s,
                            modelflag = 2, iridge = iridge)
      theta_1s_d <- result[["theta_1s"]]
      beta_1s_d  <- result[["beta_1s"]]

      muY1_d[S==s] <- X_d_tmp_1[S==s,] %*% matrix(theta_1s_d[-1]) + theta_1s_d[1]
      muD1_d[S==s] <- LogisticReg(X_d_tmp_1[S==s,] %*% matrix(beta_1s_d[-1]) + beta_1s_d[1])
    }

    vtauhat[5] <- tau(muY1 = muY1_d, muY0 = muY0_d, muD1 = muD1_d, muD0 = muD0_d,
                      A = A, S = S, Y = Y, D = D)
    vsighat[5] <- stanE(muY1 = muY1_d, muY0 = muY0_d, muD1 = muD1_d, muD0 = muD0_d,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[5])
    vstat[5] <- abs(sqrt(n)*(vtauhat[5] - tau - iPert))/vsighat[5]
    vdeci[5] <- vstat[5] > norminv(1-iq/2)

    # Nonparametric + Lasso (f)
    X_f <- cbind(X, X^2, X[,1]*X[,2], mH, mH^2, mH[,1]*mH[,2])
    muY0_f <- NaN*ones(n,1)
    muY1_f <- NaN*ones(n,1)
    muD0_f <- NaN*ones(n,1)
    muD1_f <- NaN*ones(n,1)

    for (s in 1:max(S)) {
      #for the case A==0
      X_f_tmp_0 <- X_f
      X_f_tmp_0_a <- X_f_tmp_0[S==s & A==0,]
      X_f_tmp_0 <- X_f_tmp_0[, colSums(X_f_tmp_0_a^2) != 0]

      result <- LinearLogit(Y = Y, D = D, A = A, X = X_f_tmp_0, S = S, s = s,
                            modelflag = 3, iridge = iridge)

      theta_0s_f <- result[["theta_0s"]]
      beta_0s_f  <- result[["beta_0s"]]

      muY0_f[S==s] <- X_f_tmp_0[S==s,] %*% matrix(theta_0s_f[-1]) + theta_0s_f[1]
      muD0_f[S==s] <- LogisticReg(X_f_tmp_0[S==s,] %*% matrix(beta_0s_f[-1]) + beta_0s_f[1])

      # for the case A==1
      X_f_tmp_1 <- X_f
      X_f_tmp_1_a <- X_f_tmp_1[S==s & A==1,]
      X_f_tmp_1 <- X_f_tmp_1[, colSums(X_f_tmp_1_a^2) != 0]

      result <- LinearLogit(Y = Y, D = D, A = A, X = X_f_tmp_1, S = S, s = s,
                            modelflag = 3, iridge = iridge)
      theta_1s_f <- result[["theta_1s"]]
      beta_1s_f  <- result[["beta_1s"]]

      muY1_f[S==s] <- X_f_tmp_1[S==s,] %*% matrix(theta_1s_f[-1]) + theta_1s_f[1]
      muD1_f[S==s] <- LogisticReg(X_f_tmp_1[S==s,] %*% matrix(beta_1s_f[-1]) + beta_1s_f[1])
    }

    vtauhat[8] <- tau(muY1 = muY1_f, muY0 = muY0_f, muD1 = muD1_f, muD0 = muD0_f,
                      A = A, S = S, Y = Y, D = D)
    vsighat[8] <- stanE(muY1 = muY1_f, muY0 = muY0_f, muD1 = muD1_f, muD0 = muD0_f,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[8])
    vstat[8] <- abs(sqrt(n)*(vtauhat[8] - tau - iPert))/vsighat[8]
    vdeci[8] <- vstat[8] > norminv(1-iq/2)

  } else if (dgptype == 3) {

    #High-dimensional (e)
    muY0_e <- NaN*ones(n,1)
    muY1_e <- NaN*ones(n,1)
    muD0_e <- NaN*ones(n,1)
    muD1_e <- NaN*ones(n,1)

    for (s in 1:max(S)) {

      result <-LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = s, modelflag = 3, iridge = iridge)

      theta_0s_e <- result[["theta_0s"]]
      theta_1s_e <- result[["theta_1s"]]
      beta_0s_e <- result[["beta_0s"]]
      beta_1s_e <- result[["beta_1s"]]

      muY0_e[S==s] <- X[S==s,] %*% matrix(theta_0s_e[-1]) + theta_0s_e[1]
      muY1_e[S==s] <- X[S==s,] %*% matrix(theta_1s_e[-1]) + theta_1s_e[1]
      muD0_e[S==s] <- LogisticReg(X[S==s,] %*% matrix(beta_0s_e[-1]) + beta_0s_e[1])
      muD1_e[S==s] <- LogisticReg(X[S==s,] %*% matrix(beta_1s_e[-1]) + beta_1s_e[1])
    }

    #return(tau(muY1 = muY1_e, muY0 = muY0_e, muD1 = muD1_e, muD0 = muD0_e,
    #           A = A, S = S, Y = Y, D = D))

    vtauhat[6] <- tau(muY1 = muY1_e, muY0 = muY0_e, muD1 = muD1_e, muD0 = muD0_e,
                      A = A, S = S, Y = Y, D = D)

    vsighat[6] <- stanE(muY1 = muY1_e, muY0 = muY0_e, muD1 = muD1_e, muD0 = muD0_e,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[6])
    vstat[6] <- abs(sqrt(n)*(vtauhat[6] - tau - iPert))/vsighat[6]
    vdeci[6] <- vstat[6] > norminv(1-iq/2)
  }

  # Monte Carlo Tracking
  print(str_c("Currently at ", ii, " th sample!"))

  #return results
  output_result <- list(vtauhat = vtauhat, vsighat = vsighat,
                        vstat = vstat, vdeci = vdeci)

  return(output_result)
}
