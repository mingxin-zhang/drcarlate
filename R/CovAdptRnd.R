#' @title Generate treatment assignment under various CARs
#' @description Generate treatment assignment under various CARs.
#' @param rndflag Index of the assignment rule. 1 for SRS; 2 for WEI; 3 for BCD; 4 for SBR
#' @param S A nx1 vector.
#' @param pi Targeted assignment probability across strata. It should be a vector with the length of max(S),
#' It should be noted that the treatment assignment process is independent of pi when rndflag == 2 or 3.
#' @return A nx1 treatment assignment vector generated according to the specified method.
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' CovAdptRnd(rndflag = 1, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
#' CovAdptRnd(rndflag = 2, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
#' CovAdptRnd(rndflag = 3, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
#' CovAdptRnd(rndflag = 4, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
#'
CovAdptRnd <- function(rndflag, S, pi) {
  n <- length(S)
  g <- max(S)
  A <- zeros(n,1)

  if (g != length(pi)) {
    stop("g not equal length(vpi)!")
  }

  if (rndflag == 1) {

    for (s in 1:g) {
      A[S==s] <- (rand(length(A[S==s]), 1) < pi[s])
    }

    return(A+0)

  } else if (rndflag == 2) {
    # WEI does not need pi. pi is endogenously determined (i.e., pi(s)=0.5 for all s)
    A[1,1] <- (rand(1,1) < 0.5)

    for (k in 2:n) {
      D <- 2*sum((A[1:(k-1),]-0.5) * (S[1:(k-1),1] == S[k,1]))

      if (is.nan(D/(sum(S[1:(k-1),1] == S[k,1])))) {
        A[k,1] <- 0
      } else {
      p <- (1-D/(sum(S[1:(k-1),1] == S[k,1])))/2
      A[k,1] <- (rand(1,1) < p)
      }
    }
    return(A)

  } else if (rndflag == 3) {
    # BCD does not need pi. pi is endogenously determined (i.e., pi(s)=0.5 for all s)
    A[1,1] <- (rand(1,1) < 0.5)

    for (k in 2:n) {
      D <- sum((A[1:(k-1),]-0.5) * (S[1:(k-1),1] == S[k,1]))
      p <- 0.5*(D==0) + 0.75*(D<0) + 0.25*(D>0)
      A[k,1] <- (rand(1,1) < p)
    }

    return(A)

  } else if (rndflag == 4) {
      A <- zeros(n,1)
      for (s in 1:g) {
        count1 <- (S==s)
        idx <- zeros(n,1)

        for (i in 1:n) {
          idx[i,1] <- sum(count1[1:i,1])
        }

        count2 <- cbind(matrix(1:n), count1*idx)
        ns <- sum(count1)
        ms <- floor(pi[s]*ns)
        idxa <- randsample(ns,ms,replacement = FALSE)
        idxd <- zeros(n,1)

        for (i in 1:ms) {
          idxd <- idxd+(count2[,2] == idxa[i])
        }

        A[idxd==1,1] <- 1
      }

      return(A)
    }
  }

