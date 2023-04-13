##' Print contents of iCenJMMLSM object.
##' @title Print iCenJMMLSM
##' @param x Object of class 'iCenJMMLSM'.
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{iCenJMMLSM}}
##' @author Shanpeng Li
##' @export
##'
print.iCenJMMLSM <- function(x, digits = 4, ...) {
  if (!inherits(x, "iCenJMMLSM"))
    stop("Not a legitimate \"iCenJMMLSM\" object.")

  cat("\nCall:\n", sprintf(format(paste(deparse(x$mycall, width.cutoff = 500), collapse = ""))), "\n\n")

  cat("Data Summary:\n")
  cat("Number of observations:", nrow(x$Ydata), "\n")
  cat("Number of groups:", nrow(x$Tdata), "\n\n")
  cat("Proportion of events:", round(x$PropComp[2, 2]/nrow(Tdata)*100, 2), "%\n")
  cat("\nNumerical intergration:\n")
  cat("Method: Standard Guass-Hermite quadrature\n")
  cat("Number of quadrature points: ", x$quadpoint, "\n")
  cat("Smoothing hazard kernel function: ", x$hazard.kernel, "\n")
  cat("Global bandwidth for kernel smoothing: ", x$c*x$H0[nrow(x$H0),1]/(8*nrow(x$H0)^0.2), "\n")
  cat("\nModel Type: joint modeling of longitudinal continuous and survival data with the presence of intra-individual variability", "\n\n")
  cat("Model summary:\n")
  cat("Longitudinal process: Mixed effects location scale model\n")
  cat("Event process: Cox proportional hazard model with non-parametric baseline hazard\n\n")
  ##cat("Loglikelihood: ", x$loglike, "\n\n")
  cat("Fixed effects in mean of longitudinal submodel: ",
      sprintf(format(paste(deparse(x$FunCall_long, width.cutoff = 500), collapse=""))), "\n")
  cat("\n")
  
  dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
  colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
  dat[, 1:3] <- round(dat[, 1:3], digits+1)
  dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
  print(dat)
  
  cat("\nFixed effects in variance of longitudinal submodel: \n",
      sprintf(format(paste(deparse(x$FunCall_longVar, width.cutoff = 500), collapse=""))), "\n")
  
  dat <- data.frame(x$tau, x$setau, x$tau/x$setau, 2 * pnorm(-abs(x$tau/x$setau)))
  colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
  dat[, 1:3] <- round(dat[, 1:3], digits+1)
  dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
  print(dat)
  
  cat("\nSurvival sub-model fixed effects: ",
      sprintf(format(paste(deparse(x$FunCall_survival, width.cutoff = 500), collapse=""))), "\n")
  cat("\n")
  dat <- data.frame(x$gamma, x$segamma, x$gamma/x$segamma, 2 * pnorm(-abs(x$gamma/x$segamma)))
  colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
  dat[, 1:3] <- round(dat[, 1:3], digits+1)
  dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
  print(dat)
  
  cat("\n Association parameters:                 \n")
  dat <- data.frame(x$alpha, x$sealpha, x$alpha/x$sealpha, 2 * pnorm(-abs(x$alpha/x$sealpha)))
  colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
  if (length(x$alpha) == 2) rownames(dat) <- c("alpha_b0", "alpha_w")
  if (length(x$alpha) == 3) rownames(dat) <- c("alpha_b0", "alpha_b1", "alpha_w")
  dat[, 1:3] <- round(dat[, 1:3], digits+1)
  dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
  print(dat)
  cat("\n")
  
  cat("\nRandom effects:                 \n")
  cat("  Formula:", format(as.formula(x$random)), "\n")
  
  
  if (nrow(x$Sig) == 2) {
    dat <- matrix(0, nrow = 3, ncol = 4)
    dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
    dat[2, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
    dat[3, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
    dat <- as.data.frame(dat)
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(dat) <- c("(Intercept)", "(Intercept):Residual", "Residual")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
  } else {
    dat <- matrix(0, nrow = 6, ncol = 4)
    dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
    dat[2, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
    dat[3, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
    for (i in 1:nrow(x$Sig)) dat[3+i, ] <- c(x$Sig[i,3], x$seSig[i,3], x$Sig[i,3]/x$seSig[i,3],
                                             2 * pnorm(-abs(x$Sig[i,3]/x$seSig[i,3])))
    dat <- as.data.frame(dat)
    slope <- all.vars(x$random)[1]
    interslope <- paste0("(Intercept):", slope)
    residslope <- paste0(slope, ":Residual")
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(dat) <- c("(Intercept)", interslope, slope, "(Intercept):Residual", residslope, "Residual")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    
  }
  
  
}
