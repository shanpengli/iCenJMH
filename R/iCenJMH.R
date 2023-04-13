#' @useDynLib iCenJMH, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats  as.formula  pnorm  pchisq complete.cases
#' @importFrom statmod  gauss.quad
#' @importFrom utils  read.table
#' @importFrom survival coxph
#' @importFrom parallel mclapply parLapply makeCluster stopCluster
#' @importFrom dplyr left_join summarise group_by_
#' @importFrom MASS mvrnorm
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom magrittr %>%
NULL
