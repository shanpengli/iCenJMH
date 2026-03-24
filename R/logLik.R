logLik.learn <- function(bw, data) {
  logLik_learn(
    bw = bw,
    Y = data$Y,
    X = data$X,
    Z = data$Z,
    W = data$W,
    X2 = data$X2,
    CH0 = data$CH0,
    HAZ0 = data$HAZ0,
    beta = data$beta,
    tau = data$tau,
    gamma = data$gamma,
    alpha = data$alpha,
    SigInv = data$SigInv,
    logdetSig = data$logdetSig,
    status = as.integer(data$status)
  )
}

logLik <- function(data, bw) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bw[1:p1a])
  w <- bw[p1a+1]
  sigma <- exp(data$W%*%data$tau + w)
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*sigma) + 0.5*log(sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%bw) +
    t(bw)%*%solve(data$Sig)%*%bw/2 + 0.5*log(det(data$Sig))
}