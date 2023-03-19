#' @title Feasible Post Lasso Mat Tool
#' @description Under the condition of high dimensional data,
#'    the function first selects covariables through lasso regression,
#'    then performs logit regression or linear regression according to the caller's requirements,
#'    and finally returns the adjusted Lasso regression coefficient vector.
#'    Note that HH tweaks a bit to accommodate D perfect separation.
#' @param x A nxk Matrix.
#' @param y A nx1 vector.
#' @param MaxIter Maximum iteration. The default value is 30.
#' @param UpsTol Upper limit of tolerance. The default value is 1e-6.
#' @param beta0 NULL.
#' @param clusterVar NULL.
#' @param Dist The default value is normal.
#' @param link Link can be identity or logit.
#'    This determines the method used for regression with the selected write variable after lasso.
#'    See Jiang et al.(2022) for more details.
#' @param glmTol  Maximum tolerance in GLM. The default value is 1e-8.
#' @param initScale Initail scale, the default value is 0.5.
#'
#' @return A kx1 cector, the coefficients b.
#' @export
#'
#' @examples
#' set.seed(1)
#' # Notice that when we set dgptype = 3, FuncDGP will generate a high dimensional data for us.
#' DGP <- FuncDGP(dgptype = 3, rndflag = 1, n = 10000, g = 4, pi = 0.5)
#' X <- DGP$X
#' Y <- DGP$Y
#' A <- DGP$A
#' S <- DGP$S
#' D <- DGP$D
#' feasiblePostLassoMatTool(x = X[S==1 & A==0,], y = Y[S==1 & A==0,])
#' feasiblePostLassoMatTool(x = X[S==1 & A==0,], y = D[S==1 & A==0,], link = "logit")
#'
#'
feasiblePostLassoMatTool <- function(x, y,
                                     MaxIter = 30,
                                     UpsTol = 1e-6,
                                     beta0 = c(),
                                     clusterVar = c(),
                                     Dist = "normal", link = "identity",
                                     glmTol = 1e-8, initScale = 0.5) {
  n <- size(x,1)
  p <- size(x,2)

  lambda = 1.1*sqrt(n)*norminv(1-(1/log(n)/(p)))
  lambda <- lambda/(4*n)

  if (isempty(clusterVar)) { # always goes into this loop
    if (isempty(beta0)) { # always goes into this loop
      if (str_detect(link, "logit")) {
        Syx <- x*((y-mean(y)) %*% ones(1,p)) #n x p
      } else {
        Syx <- x*((y-mean(y)) %*% ones(1,p))
      }

      Ups0 <- t(matrix(sqrt(colVars(Syx)))) #1xp

      # the following line prevents zero element of Ups0
      Ups0[Ups0 < 0.001] <- 1

      lasso_result <- glmnet(x = (x / (ones(n,1) %*% Ups0)),
                  y = y,
                  family = "gaussian",
                  lambda = initScale*lambda,
                  standardize = FALSE,
                  thresh = glmTol)

      b <- matrix(lasso_result$beta)
      use <- abs(b) > 0

      if (str_detect(link, "logit")) {
        e <- y- matrix(predict.glm(object = glm(formula = y ~ ., data = data.frame(y=y,x[,use]),
                                      family = binomial(link = "logit")),
                         newdata = data.frame(x[,use]),
                         type = "response"))
      } else {
        e <- y - (x[,use] %*% (mldivide(A = x[,use], B = y)))
      }
    }

    kk <- 1
    Syx <- x * (e %*% ones(1,p))

    Ups1 <- t(matrix(sqrt(colVars(Syx))))
    Ups1[Ups1<0.001] <- 1


    while (norm(Ups0 - Ups1,"2") > UpsTol & kk < MaxIter) {
      d0 <- norm(Ups0-Ups1, "2")

      lasso_result_a <-  glmnet(x = x/(ones(n,1) %*% Ups1), y = y, family = "gaussian",
                 lambda = lambda, standardize = FALSE, thresh = glmTol)

      b <- matrix(lasso_result_a$beta)
      use <- abs(b) > 0

      if (sum(use) < (size(x,2)-1)) {
        if (str_detect(link, "logit")) {
          e <-  y - matrix(predict.glm(object = glm(formula = y~., data = data.frame(y=y,x[,use]),
                                        family = binomial(link = "logit")),
                           newdata = data.frame(x[,use]),
                           type = "response"))
        } else {
          e <- y - x[,use] %*% mldivide(A = x[,use], B = y)
        }

        Ups0 <- Ups1

        Syx <- x * (e %*% ones(1,p))

        Ups1 <- t(matrix(sqrt(colVars(Syx))))
        Ups1[Ups1<0.001] <- 1

        kk <- kk+1
        d1 <- norm(Ups0 - Ups1, "2")

        if (d1 == d0) {
          Ups0 <- Ups1
        }
      } else {
        kk <- MaxIter + 1
      }
    }

    b <- b/t(Ups1)

    return(b)
  }
}
