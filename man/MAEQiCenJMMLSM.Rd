% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAEQiCenJMMLSM.R
\name{MAEQiCenJMMLSM}
\alias{MAEQiCenJMMLSM}
\title{A metric of prediction accuracy of joint model by comparing the predicted risk
with the empirical risks stratified on different predicted risk group.}
\usage{
MAEQiCenJMMLSM(
  seed = 100,
  object,
  landmark.time = NULL,
  horizon.time = NULL,
  obs.time = NULL,
  quadpoint = NULL,
  maxiter = 1000,
  n.cv = 3,
  quantile.width = 0.25,
  ...
)
}
\arguments{
\item{seed}{a numeric value of seed to be specified for cross validation.}

\item{object}{object of class 'iCenJMMLSM'.}

\item{landmark.time}{a numeric value of time for which dynamic prediction starts..}

\item{horizon.time}{a numeric vector of future times for which predicted probabilities are to be computed.}

\item{obs.time}{a character string of specifying a longitudinal time variable.}

\item{quadpoint}{the number of standard Gauss-Hermite quadrature points.}

\item{maxiter}{the maximum number of iterations of the EM algorithm that the 
function will perform. Default is 10000.}

\item{n.cv}{number of folds for cross validation. Default is 3.}

\item{quantile.width}{a numeric value of width of quantile to be specified. Default is 0.25.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
a list of matrices with conditional probabilities for subjects.
}
\description{
A metric of prediction accuracy of joint model by comparing the predicted risk
with the empirical risks stratified on different predicted risk group.
}
\seealso{
\code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
}
\author{
Shanpeng Li \email{lishanpeng0913@ucla.edu}
}
