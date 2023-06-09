#' Simulates the data for ATE estimators
#' @description
#' ATEDGP is the version of FuncDGP under full compliance.
#' @param dgptype A scalar. 1, 2, 3 (Almost the same as 1-3 in the paper except that it does not have the DGP for D(1) or D(0)).
#' @param rndflag A scalar. method of covariate-adaptive randomization. 1-SRS; 2-WEI; 3-BCD; 4-SBR.
#' @param n Sample size.
#' @param g Number of strata. The authors set g = 4 in the Jiang et al. (2022).
#' @param pi A g x 1 vector. Targeted assignment probabilities across strata.
#'
#' @return ATEDGP returns a list containing 7 nx1 vectors named Y, X, S, A, Y1, Y0 and D.
#'    These seven vectors are the same as defined in Jiang et al. (2022).
#'    Note that vector X does not contain the constant term.
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' ATEDGP(dgptype = 1, rndflag = 1, n = 200, g = 4, pi = c(0.5, 0.5, 0.5, 0.5))

ATEDGP <- function(dgptype, rndflag, n, g, pi) {
  if (dgptype == 1) {
    Z <- (matrix(rbeta(n = n,shape1 = 2,shape2 = 2))-0.5)*sqrt(20)
    X <- cbind((rand(n,1)-0.5)*4, randn(n,1)+Z)
    grid <- seq(from = -0.5*sqrt(20), to = 0.5*sqrt(20), by = sqrt(20)/g)

    S <- matrix(rowSums(repmat(Z,1,g) < repmat(grid[2:length(grid)],n,1)))
    A <- CovAdptRnd(rndflag = rndflag, S = S, pi = pi)
    Sigma <- toeplitz((0.5)^(0:3))
    Err <- mvrnorm(n = n,mu = zeros(1,4), Sigma = Sigma)

    Y1 <- 2 + 0.7*X[,1]^2 + X[,2] + 4*Z + Err[,1]
    Y0 <- 1 + 0.7*X[,1]^2 + X[,2] + 4*Z + Err[,2]
    #D0 <- -1 + 0.5*X[,1]^2 - 0.5*X[,2]^2- 0.5*Z^2 > 3*Err[,3]
    #D1 <- 1.3 + 0.5*X[,1]^2 - 0.5*X[,2]^2 - 0.5*Z^2 >3*Err[,4]
    #D1[D0==1] <- 1
    #D <- D1*A + D0*(1-A)
    D <- A
    Y <- Y1*D + Y0*(1-D)

    result_list <- list(Y = Y, X = X, S = S, A = A, Y1 = Y1, Y0 = Y0, D = D)
    return(result_list)

  } else if (dgptype == 2) {

    Z <- 4*(rand(n,1) - 0.5)
    grid <- seq(from = -2, to = 2, by = 4/g)
    X <- cbind((rand(n,1)-0.5)*4, randn(n,1))

    S <- matrix(rowSums(repmat(Z,1,g)<repmat(grid[2:length(grid)],n,1)))
    A <- CovAdptRnd(rndflag = rndflag, S = S, pi = pi)
    Sigma <- toeplitz((0.5)^(0:3))
    Err <- mvrnorm(n = n,mu = zeros(1,4), Sigma = Sigma)

    Y1 <- 2 - 0.8*X[,2]*X[,1] + Z^2 + Z*X[,1] + Err[,1]
    Y0 <- 1 - 0.8*X[,2]*X[,1] + Z^2 + Z*X[,1] + Err[,2]
    #D0 <- -1 + 0.5*X[,1]^2 - 0.5*X[,2]^2 - 0.5*Z^2 > Err[,3]*3
    #D1 <- 1 + 0.5*X[,1]^2 - 0.5*X[,2]^2 - 0.5*Z^2 > Err[,4]*3
    #D1[D0==1] <- 1
    #D <- D1*A + D0*(1-A)
    D <- A
    Y <- Y1*D + Y0*(1-D)

    result_list <- list(Y = Y, X = X, S = S, A = A, Y1 = Y1, Y0 = Y0,D = D)
    return(result_list)

  } else if (dgptype == 3) {

    K <- 20
    beta <- sqrt(6)/matrix(1:K)^2
    #gamma <- -2/matrix(1:K)^2
    Omega <- toeplitz(0.5^(0:(K-1)))
    Sigma <- toeplitz(0.5^(0:3))

    X <- mvrnorm(n = n, mu = zeros(1,K), Omega)
    Z <- (matrix(rbeta(n = n, shape1 = 2, shape2 = 2)) - 0.5)*sqrt(20)
    grid <- seq(from = -0.5*sqrt(20), to = 0.5*sqrt(20), by = sqrt(20)/g)
    S <- matrix(rowSums(repmat(Z,1,g) < repmat(grid[2:length(grid)],n,1)))
    A <- CovAdptRnd(rndflag = rndflag, S, pi)
    Err <- mvrnorm(n = n,mu = zeros(1,4), Sigma = Sigma)

    Y1 <- 2 + Z + X %*% beta + Err[,1]
    Y0 <- 1 + Z + X %*% beta + Err[,2]
    #D0 <- -1 + X %*% gamma - Z > Err[,3]*sqrt(7)
    #D1 <- 2 + X %*% gamma - Z > Err[,4]*sqrt(7)
    #D1[D0==1] <- 1
    #D <- D1*A + D0*(1-A)
    D <- A
    Y <- Y1*D + Y0*(1-D)

    result_list <- list(Y = Y, X = X, S = S, A = A, Y1 = Y1, Y0 = Y0, D = D)
    return(result_list)

  } else {
    stop("The dgptype only takes values {1,2,3}!")
  }
}
