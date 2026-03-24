worker_postmode <- function(task_row, pSLR, pStol, Y, W, X1, Z, X2, H0Y, status,
                            beta, tau, gamma, alpha, SigInv, logdetSig, nsig) {
  
  i0  <- task_row$obs_start
  t0  <- task_row$surv_start
  nii <- task_row$ni
  t   <- task_row$t
  
  rowt <- t0 + t - 1
  
  if (pSLR[rowt, 4] <= pStol || is.nan(pSLR[rowt, 4])) {
    return(list(
      rowt = rowt,
      par = rep(0, nsig),
      hess_inv = matrix(0, nsig, nsig)
    ))
  }
  
  rows_i <- (i0 + (t - 1) * nii):(i0 + t * nii - 1)
  
  subY  <- Y[rows_i]
  subW  <- W[rows_i, , drop = FALSE]
  subX1 <- X1[rows_i, , drop = FALSE]
  subZ  <- Z[rows_i, , drop = FALSE]
  subX2 <- matrix(X2[rowt, ], nrow = 1)
  CH0   <- H0Y[rowt, 3]
  HAZ0  <- H0Y[rowt, 4]
  substatus <- status[rowt]
  
  data <- list(
    Y = subY,
    X = subX1,
    Z = subZ,
    W = subW,
    X2 = subX2,
    CH0 = CH0,
    HAZ0 = HAZ0,
    beta = beta,
    tau = tau,
    gamma = gamma,
    alpha = alpha,
    SigInv = SigInv,
    logdetSig = logdetSig,
    status = substatus
  )
  
  opt <- optim(
    rep(0, nsig),
    logLik.learn,
    data = data,
    method = "BFGS",
    hessian = TRUE
  )
  
  return(list(
    rowt = rowt,
    par = opt$par,
    hess_inv = solve(opt$hessian))
  )
}