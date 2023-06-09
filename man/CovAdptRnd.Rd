% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovAdptRnd.R
\name{CovAdptRnd}
\alias{CovAdptRnd}
\title{Generate treatment assignment under various CARs}
\usage{
CovAdptRnd(rndflag, S, pi)
}
\arguments{
\item{rndflag}{Index of the assignment rule. 1 for SRS; 2 for WEI; 3 for BCD; 4 for SBR}

\item{S}{A nx1 vector.}

\item{pi}{Targeted assignment probability across strata. It should be a vector with the length of max(S),
It should be noted that the treatment assignment process is independent of pi when rndflag == 2 or 3.}
}
\value{
A nx1 treatment assignment vector generated according to the specified method.
}
\description{
Generate treatment assignment under various CARs.
}
\examples{
CovAdptRnd(rndflag = 1, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
CovAdptRnd(rndflag = 2, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
CovAdptRnd(rndflag = 3, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))
CovAdptRnd(rndflag = 4, S = matrix(sample(1:4,100,TRUE)), pi = c(0.5, 0.5, 0.5, 0.5))

}
\references{
Jiang L, Linton O B, Tang H, Zhang Y. Improving estimation efficiency via regression-adjustment in covariate-adaptive randomizations with imperfect compliance [J]. 2022.
}
