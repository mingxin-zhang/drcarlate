#' @title Linear Regression or Logit Regression
#' @description LinearLogit generates estimated pseudo true values for parametric models.
#'    Different estimation strategies are adopted according to different values of modelflag.
#'    See Jiang et al. (2022) for more details about different strategies.
#' @param Y The outcome vector. A nx1 vector.
#' @param D A nx1 vector.
#' @param A The treatment assignment. A nx1 vector.
#' @param X Extra covariate matrix, A nxK matrix without constant.
#' @param S The strata variable.
#' @param s A particular stratum.
#' @param modelflag Its value ranges from characters 1, 2, and 3, respectively declaring different estimation strategies.
#' 1-L; 2-NL; 3-R.
#' @param iridge A scalar. The penalization parameter in ridge regression.
#'
#' @return theta_0s, theta_1s, beta_0s, beta_1s are estimated coefficients vectors.
#'    The dimension is Kx1 if modelflag = 1; (K+1)x1 if modelflag = 2 or 3.
#' @export
#' @references Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
#' @examples
#' #' set.seed(1)
#' DGP <- FuncDGP(dgptype = 3, rndflag = 1, n = 10000, g = 4, pi = c(0.5, 0.5, 0.5, 0.5))
#' X <- DGP$X
#' Y <- DGP$Y
#' A <- DGP$A
#' S <- DGP$S
#' D <- DGP$D
#' LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = 1, modelflag = 1, iridge = 0.001)
#' LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = 2, modelflag = 2, iridge = 0.001)
#' LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = 3, modelflag = 3, iridge = 0.001)
#' LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = 4, modelflag = 3, iridge = 0.001)
#'
LinearLogit <- function(Y, D, A, X, S, s, modelflag, iridge) {
  n <- dim(X)[1]
  K <- dim(X)[2]

  if (modelflag == 1) {
    X_0s <-  X - repmat(a = colMeans(X[S==s & A ==0,]), n = n, m = 1)
    X_1s <-  X - repmat(a = colMeans(X[S==s & A ==1,]), n = n, m = 1)

    if (min(eig(t(X_0s[S==s & A==0,]) %*% X_0s[S==s & A==0,])) < 0.001) {
      theta_0s <- inv(t(X_0s[S==s & A==0,]) %*% X_0s[S==s & A==0,] + iridge*eye(K)) %*%
        (t(X_0s[S==s & A==0,]) %*% matrix(Y[S==s & A==0]))

      beta_0s <- inv(t((X_0s[S==s & A==0,])) %*% X_0s[S==s & A==0,] + iridge*eye(K)) %*%
        (t(X_0s[S==s & A==0,]) %*% matrix(D[S==s & A==0]))
    } else {
      theta_0s <- inv(t((X_0s[S==s & A==0,])) %*% X_0s[S==s & A==0,]) %*%
        (t(X_0s[S==s & A==0,]) %*% matrix(Y[S==s & A==0]))

      beta_0s <- inv(t((X_0s[S==s & A==0,])) %*% X_0s[S==s & A==0,]) %*%
        (t(X_0s[S==s & A==0,]) %*% matrix(D[S==s & A==0,]))
    }

    if (any(is.nan(theta_0s))) {
      theta_0s <- zeros(n = size(theta_0s,1), m = size(theta_0s,2))
    }
    if (any(is.nan(beta_0s))) {
      beta_0s <- zeros(n = size(beta_0s,1), m = size(beta_0s,2))
    }

    if (min(eig(t(X_1s[S==s & A==1,]) %*% X_1s[S==s & A==1,])) < 0.001) {
      theta_1s <- inv(t((X_1s[S==s & A==1,])) %*% X_1s[S==s & A==1,] + iridge*eye(K)) %*%
        (t(X_1s[S==s & A==1,]) %*% matrix(Y[S==s & A==1]))

      beta_1s <- inv(t((X_1s[S==s & A==1,])) %*% X_1s[S==s & A==1,] + iridge*eye(K)) %*%
        (t(X_1s[S==s & A==1,]) %*% matrix(D[S==s & A==1]))
    } else {
      theta_1s <- inv(t((X_1s[S==s & A==1,])) %*% X_1s[S==s & A==1,]) %*%
        (t(X_1s[S==s & A==1,]) %*% matrix(Y[S==s & A==1]))

      beta_1s <- inv(t((X_1s[S==s & A==1,])) %*% X_1s[S==s & A==1,]) %*%
        (t(X_1s[S==s & A==1,]) %*% matrix(D[S==s & A==1]))
    }

    if (any(is.nan(theta_1s))) {
      theta_1s <- zeros(n = size(theta_1s)[1], m = size(theta_1s)[2])
    }
    if (any(is.nan(beta_1s))) {
      beta_1s <- zeros(n = size(beta_1s)[1], m = size(beta_1s)[2])
    }

    result_list <- list(theta_0s = theta_0s,
                        theta_1s = theta_1s,
                        beta_0s = beta_0s,
                        beta_1s = beta_1s)

    return(result_list)

  } else if (modelflag == 2) {

    # test
    #print(min(eig(t(X[S==s & A==0,]) %*% X[S==s & A==0,])))
    #print(min(eig(t(X[S==s & A==0,]) %*% X[S==s & A==0,])))
    #print("==========================================================")
    #print(s)
    #print(size(X))
    #print(data.frame(S=S, A=A) %>%  group_by(S,A) %>% summarise(n=n()))
    #print(X[S==s & A==0,])
    #print(X[S==s & A==1,])
    #print("===========================================================")
    # end


    if (min(eig(t(X[S==s & A==0,]) %*% X[S==s & A==0,])) < 0.001) {
      mTmp0 <- cbind(ones(n,1), X)
      theta_0s <- inv(t(mTmp0[S==s & A==0,]) %*% mTmp0[S==s & A==0,] + iridge*eye(K+1)) %*%
        (t(mTmp0[S==s & A==0,]) %*% Y[S==s & A==0,])
      beta_0s <- matrix((glm(formula = D~.-1,
                     data = data.frame(D=D[S==s & A==0,], mTmp0[S==s & A==0,]),
                     family = binomial(link = "logit")))$coefficients)
    } else {
      mTmp0 <- cbind(ones(n,1), X)
      theta_0s <- inv(t(mTmp0[S==s & A==0,]) %*% mTmp0[S==s & A==0,]) %*%
        (t(mTmp0[S==s & A==0,]) %*% matrix(Y[S==s & A==0]))
      beta_0s <- matrix((glm(formula = D~.-1,
                             data = data.frame(D=D[S==s & A==0,], mTmp0[S==s & A==0,]),
                             family = binomial(link = "logit")))$coefficients)
    }

    if (any(is.nan(theta_0s))) {
      theta_0s <- zeros(size(theta_0s)[1], size(theta_0s)[2])
    }

    if (min(eig(t(X[S==s & A==1,]) %*% X[S==s & A==1,])) < 0.001) {
      mTmp1 <- cbind(ones(n,1), X)
      theta_1s <- inv(t(mTmp1[S==s & A==1,]) %*% mTmp1[S==s & A==1,] + iridge*eye(K+1)) %*%
        (t(mTmp1[S==s & A==1,]) %*% matrix(Y[S==s & A==1]))
      beta_1s <- matrix((glm(formula = D~.-1,
                             data = data.frame(D=D[S==s & A==1,], mTmp1[S==s & A==1,]),
                             family = binomial(link = "logit")))$coefficients)
    } else {
      mTmp1 <- cbind(ones(n,1), X)
      theta_1s <- inv(t(mTmp1[S==s & A==1,]) %*% mTmp1[S==s & A==1,]) %*%
        (t(mTmp1[S==s & A==1,]) %*% matrix(Y[S==s & A==1]))
      beta_1s <- matrix((glm(formula = D~.-1,
                             data = data.frame(D=D[S==s & A==1,], mTmp1[S==s & A==1,]),
                             family = binomial(link = "logit")))$coefficients)
    }

    if (any(is.nan(theta_1s))) {
      theta_1s <- zeros(size(theta_1s,1), size(theta_1s,2))
    }

    result_list <- list(theta_0s = theta_0s,
                        theta_1s = theta_1s,
                        beta_0s = beta_0s,
                        beta_1s = beta_1s)

    return(result_list)

  } else if(modelflag == 3) {
    X <- cbind(ones(n,1), X)

    if (size(X[S==s & A==0, ],1) < 6) {
      theta_0s <- zeros(K+1,1)
      beta_0s <- matrix(c(-100, zeros(K, 1)))
    } else {
      theta_0s <- feasiblePostLassoMatTool(y = matrix(Y[S==s & A==0]), x = X[S==s & A==0, ])
      if (identical(matrix(D[S==s & A==0,]),
                    zeros(size(matrix(D[S==s & A==0,]),1), size(matrix(D[S==s & A==0,]),2))) |
          identical(D[S==s & A==0,],
                              ones(size(matrix(D[S==s & A==0,]),1), size(matrix(D[S==s & A==0,]),2)))) {
        beta_0s <- matrix(glm(formula = D~.-1,
                       data = data.frame(D=D[S==s & A==0,],
                                         ones(size(D[S==s & A==0,],1),
                                              size(D[S==s & A==0,],2))),
                       family = binomial(link = "logit"))$coefficients)
        beta_0s <- rbind(beta_0s, zeros(K,1))
      } else {
        beta_0s <- feasiblePostLassoMatTool(y = matrix(D[S==s & A==0,]), x = X[S==s & A==0,],
                                            link = "logit")
      }
    }

    if (any(is.nan(theta_0s))) {
      theta_0s <- zeros(size(theta_0s,1), size(theta_0s,2))
    }

    if (size(X[S==s & A==1,],1) < 6) {
      theta_1s <- zeros(K+1,1)
      beta_1s <- matrix(c(-100, zeros(K, 1)))
    } else {
      theta_1s <- feasiblePostLassoMatTool(y = matrix(Y[S==s & A==1]), x = X[S==s & A==1, ])
      if (identical(matrix(D[S==s & A==1]),
                    zeros(size(matrix(D[S==s & A==1]),1), size(matrix(D[S==s & A==1]),2))) |
          identical(matrix(D[S==s & A==1]),
                              ones(size(matrix(D[S==s & A==1]),1), size(matrix(D[S==s & A==1]),2)))) {
        beta_1s <- matrix(glm(formula = D~.-1,
                       data = data.frame(D=matrix(D[S==s & A==1]),
                                         ones(size(matrix(D[S==s & A==1]),1),
                                              size(matrix(D[S==s & A==1]),2))),
                       family = binomial(link = "logit"))$coefficients)
        beta_1s <- rbind(beta_1s, zeros(K,1))
      } else {
        beta_1s <- feasiblePostLassoMatTool(y = matrix(D[S==s & A==1]), x = X[S==s & A==1,],
                                            link = "logit")
      }
    }
      if (any(is.nan(theta_1s))) {
        theta_1s <- zeros(size(theta_1s,1), size(theta_1s,2))
      }

    result_list <- list(theta_0s = theta_0s,
                        theta_1s = theta_1s,
                        beta_0s = beta_0s,
                        beta_1s = beta_1s)

    return(result_list)

  } else {
    print("Modelflag can take only 1, 2, 3!")
  }
}


