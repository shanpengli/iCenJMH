##' @export
##'

bootsfitRI <- function(i, seed = 99, n = 100, phi = 0.04,
                       nc = 100,
                       covbw = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2),
                       lambda = 0.3, lambdaC = 0.05,
                       Cmin = 1,
                       Cmax = 5,
                       gamma = c(-0.05, 0.2, -0.1),
                       alpha = c(0.5, -0.5),
                       beta = c(5, 1, 2, -3, 3),
                       tau = c(-0.5, -0.1, -0.2, 0.8, 0.4),
                       increment = 2, maxiter = 1000,
                       quadpoint = 15) {
  
  data <- SimJMdata(seed = seed + i, n = n, phi = phi,
                    nc = nc,
                    covbw = covbw,
                    lambda = lambda, lambdaC = lambdaC,
                    Cmin = Cmin,
                    Cmax = Cmax,
                    gamma = gamma,
                    alpha = alpha,
                    beta = beta,
                    tau = tau,
                    increment = increment)
  
  Ydata <- data$Ydata
  Sdata <- data$Sdata
  Tdata <- data$Tdata
  
  iCen.info <- GetSupport(iCen.data = Sdata, iCen.tL = "Li", iCen.tR = "Ri",
                          ID = "ID", S = "Stime", weight = "Stime.w", weight.ID = "w.ID")

  a <- proc.time()
  fit <- try(iCenJMMLSM(Ydata = Ydata, Tdata = Tdata,
                        long.formula = Yij ~ X11 + X12,
                        surv.formula = Surv(Ei, status) ~ X21 + X22,
                        variance.formula = ~ X11 + X12,
                        random = ~ 1|ID,
                        timeVar = "Oij",
                        iCen.info = iCen.info,
                        maxiter = maxiter, epsilon = 1e-04,
                        epsilonH0 = 1e-03,
                        quadpoint = quadpoint, print.para = TRUE,
                        initial.para = TRUE), silent = TRUE)
  b <- proc.time()
  time <- (b - a)[3]
  
  coef <- vector()
  coefSE <- vector()
  count <- 1
  if ('try-error' %in% class(fit)) {
    coef <- rep(NA, 20)
    coefSE <- rep(NA, 18)
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
  } else if (is.null(fit$beta)) {
    coef <- rep(NA, 20)
    coefSE <- rep(NA, 18)
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
    return(coef)
  } else if (fit$iter == maxiter) {
    coef <- rep(NA, 20)
    coefSE <- rep(NA, 18)
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
    return(coef)
  } else {
    for (j in 1:length(fit$beta)) {
      coef[count] <- fit$beta[j]
      coefSE[count] <- fit$sebeta[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$tau)) {
      coef[count] <- fit$tau[j]
      coefSE[count] <- fit$setau[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$gamma)) {
      coef[count] <- fit$gamma[j]
      coefSE[count] <- fit$segamma[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$alpha)) {
      coef[count] <- fit$alpha[j]
      coefSE[count] <- fit$sealpha[j]
      count <- count + 1
    }
    
    for (j in 1:nrow(fit$Sig)) {
      coef[count] <- fit$Sig[j, j]
      coefSE[count] <- fit$seSig[j, j]
      count <- count + 1
    }
    
    coef[count] <- fit$Sig[1, 2]
    coefSE[count] <- fit$seSig[1, 2]
    count <- count + 1
    coef[count] <- time
    count <- count + 1
    coef[count] <- fit$iter
    
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
    
    return(coef)
  }

}