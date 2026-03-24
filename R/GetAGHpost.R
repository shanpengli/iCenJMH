build_tasks <- function(nt, ni) {
  n <- length(nt)
  obs_start  <- cumsum(c(1, head(ni * nt, -1)))
  surv_start <- cumsum(c(1, head(nt, -1)))
  
  tasks <- vector("list", n)
  for (i in seq_len(n)) {
    tasks[[i]] <- data.frame(
      i = i,
      t = seq_len(nt[i]),
      obs_start = obs_start[i],
      surv_start = surv_start[i],
      ni = ni[i]
    )
  }
  do.call(rbind, tasks)
}

GetAGHpost <- function(nt, ni, pSLR, pStol, Y, W, X1, Z, X2, H0Y, status,
                       beta, tau, gamma, alpha, SigInv, logdetSig, nsig,
                       n.cores = 1) {
  
  tasks_df <- build_tasks(nt, ni)
  tasks <- split(tasks_df, seq_len(nrow(tasks_df)))
  tasks <- lapply(tasks, as.list)
  
  null_file <- if (.Platform$OS.type == "windows") "nul" else "/dev/null"
  cl <- parallel::makeCluster(n.cores, outfile = null_file)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "worker_postmode",
      "logLik.learn",
      "pSLR", "pStol", "Y", "W", "X1", "Z", "X2", "H0Y", "status",
      "beta", "tau", "gamma", "alpha", "SigInv", "logdetSig", "nsig"
    ),
    envir = environment()
  )
  
  parallel::clusterSetRNGStream(cl, iseed = 123)
  
  res <- parallel::parLapplyLB(
    cl,
    X = tasks,
    fun = function(task) {
      worker_postmode(
        task = task,
        pSLR = pSLR,
        pStol = pStol,
        Y = Y,
        W = W,
        X1 = X1,
        Z = Z,
        X2 = X2,
        H0Y = H0Y,
        status = status,
        beta = beta,
        tau = tau,
        gamma = gamma,
        alpha = alpha,
        SigInv = SigInv,
        logdetSig = logdetSig,
        nsig = nsig
      )
    }
  )
  
  posterior.mode <- matrix(0, nrow = sum(nt), ncol = nsig)
  posterior.var  <- matrix(0, nrow = sum(nt) * nsig, ncol = nsig)
  
  for (r in res) {
    rowt <- r$rowt
    posterior.mode[rowt, ] <- r$par
    posterior.var[((rowt - 1) * nsig + 1):(rowt * nsig), ] <- r$hess_inv
  }
  
  return(list(posterior.mode = posterior.mode, posterior.var = posterior.var))
}