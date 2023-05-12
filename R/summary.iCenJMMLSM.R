##' @title Summary of Fitted Joint Models
##' @name summary
##' @aliases summary.iCenJMMLSM
##' @description Produce result summaries of a joint model fit. 
##' @param object an object inheriting from class \code{iCenJMMLSM}.
##' @param process for which model (i.e., longitudinal model or survival model) to extract the estimated coefficients.
##' @param digits the number of significant digits to use when printing. Default is 4.
##' @param ... further arguments passed to or from other methods.
##' @return A table to summarize the model results.
##' @export
##'

summary.iCenJMMLSM <- function(object, process = c("longitudinal", "survival"), digits = 4, ...) {
  
  if (!inherits(object, "iCenJMMLSM"))
    stop("Use only with 'iCenJMMLSM' objects.\n")
  
  if (process == "longitudinal") {
    ##Estimates of betas
    Estimate <- object$beta
    SE <- object$sebeta
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
    out <- cbind(paste0("Mean_", rownames(out)), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    ##Estimates of tau
    Estimate <- object$tau
    SE <- object$setau
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out2 <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
    out2 <- cbind(paste0("Var_", rownames(out2)), out2)
    rownames(out2) <- NULL
    colnames(out2)[1] <- "Parameter"
    
    out3 <- rbind(out, out2)
    
    rownames(out3) <- NULL
    names(out3) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
    
    out3[, 2:ncol(out3)] <- round(out3[, 2:ncol(out3)], digits = digits)
    out3[, ncol(out3)] <- format(out3[, ncol(out3)], scientific = FALSE)
    
    return(out3)
    
  } else if (process == "survival") {
    ##gamma
    Estimate <- object$gamma
    SE <- object$segamma
    LowerLimit <- Estimate - 1.96 * SE
    expLL <- exp(LowerLimit)
    UpperLimit <- Estimate + 1.96 * SE
    expUL <- exp(UpperLimit)
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, expLL, expUL, pval)
    out <- cbind(rownames(out), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    outgamma <- out
    names(outgamma) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", 
                         "95%exp(Lower)", "95%exp(Upper)", "p-values")
    
    ##alpha
    Estimate <- object$alpha
    if (length(Estimate) == 3) names(Estimate) <- c("alpha_b0", "alpha_b1", "alpha_w0")
    if (length(Estimate) == 2) names(Estimate) <- c("alpha_b0", "alpha_w0")
    SE <- object$sealpha
    LowerLimit <- Estimate - 1.96 * SE
    expLL <- exp(LowerLimit)
    UpperLimit <- Estimate + 1.96 * SE
    expUL <- exp(UpperLimit)
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, expLL, expUL, pval)
    out <- cbind(rownames(out), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    outalpha <- out
    names(outalpha) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper",
                         "95%exp(Lower)", "95%exp(Upper)", "p-values")
    
 
    out <- rbind(outgamma, outalpha)
    
    out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
    out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)
    
    return(out)
  }
  
  
}