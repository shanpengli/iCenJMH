##' @title Print evaluation metrics for the fitted joint model
##' @name summary
##' @aliases summary.DynPredAcciCenJMMLSM
##' @param object object of class 'DynPredAcciCenJMMLSM'.
##' @param digits number of digits of decimal to be printed. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
##' @export
##' 

summary.DynPredAcciCenJMMLSM <- function (object, digits = 4, ...) {
  if (!inherits(object, "DynPredAcciCenJMMLSM"))
    stop("Use only with 'DynPredAcciCenJMMLSM' xs.\n") 
  
  if (is.null(object$metric.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    ### MAPE
    if(length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
      
      sum <- as.data.frame(matrix(0, nrow = length(object$horizon.time), ncol = 2))
      sum[, 1] <- object$horizon.time
      for (i in 1:length(object$horizon.time)) {
        for (j in 1:object$n.cv) {
          sum[i, 2] <- sum[i, 2] + mean(abs(object$metric.cv[[j]][[1]][[i]][, 1] - 
                                              object$metric.cv[[j]][[1]][[i]][, 2])) 
        }
      }
      sum[, -1] <- sum[, -1]/object$n.cv
      sum[, -1] <- round(sum[, -1], digits)
      colnames(sum) <- c("Horizon Time", "MAPE")
      ExpectedMAPE <- sum
      cat("\nMean of absolute error across quantiles of predicted risk scores at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      print(ExpectedMAPE)
    } else {
      cat("\nThe cross validation fails on MAPE. Please try using a different seed number.\n")
    }
    
    ### Brier Score
    if (length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$metric.cv[[j]][[2]]
      }
      sum <- sum/object$n.cv
      
      Brier <- round(sum[, 1], digits)
      ExpectedBrier <- data.frame(object$horizon.time, Brier)
      colnames(ExpectedBrier) <- c("Horizon Time", "Brier")
      cat("\nExpected Brier score at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      print(ExpectedBrier)
    } else {
      cat("\nThe cross validation fails on Brier score. Please try using a different seed number.\n")
    }
    
    ### AUC and Cindex
    if (length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$metric.cv[[j]][[3]]
      }
      sum <- sum/object$n.cv
      
      AUC <- round(sum[, 1], digits)
      ExpectedAUC <- data.frame(object$horizon.time, AUC)
      colnames(ExpectedAUC) <- c("Horizon Time", "AUC")
      cat("\nExpected AUC at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      print(ExpectedAUC)
    } else {
      stop("The cross validation fails on AUC. Please try using a different seed number.")
    }
    
    if (length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$metric.cv[[j]][[4]]
      }
      sum <- sum/object$n.cv
      
      Cindex <- round(sum[, 1], digits)
      ExpectedCindex <- data.frame(object$horizon.time, Cindex)
      colnames(ExpectedCindex) <- c("Horizon Time", "Cindex")
      cat("\nExpected Cindex at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      print(ExpectedCindex)
    } else {
      stop("The cross validation fails on Cindex. Please try using a different seed number.")
    }
    invisible(list(ExpectedMAPE = ExpectedMAPE,
                ExpectedBrier = ExpectedBrier,
                ExpectedAUC = ExpectedAUC,
                ExpectedCindex = ExpectedCindex))
  }
}