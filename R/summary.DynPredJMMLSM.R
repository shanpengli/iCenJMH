##' @title Print DynPredJMMLSM
##' @name summary.DynPredJMMLSM
##' @aliases summary.DynPredJMMLSM
##' @param object object of class 'DynPredJMMLSM'.
##' @param digits number of digits of decimal to be printed. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
##' @export
##' 

summary.DynPredJMMLSM <- function (object, digits = 4, ...) {
  if (!inherits(object, "DynPredJMMLSM"))
    stop("Use only with 'DynPredJMMLSM' xs.\n") 
  
  if (is.null(object$metric.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    metric <- object$metric
    
    if (metric == "Brier Score") {
      
      if (length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
        sum <- 0
        for (j in 1:object$n.cv) {
          sum <- sum + object$metric.cv[[j]]
        }
        sum <- sum/object$n.cv
        
        Brier <- round(sum[, 1], digits)
        ExpectedBrier <- data.frame(object$horizon.time, Brier)
        colnames(ExpectedBrier) <- c("Horizon Time", "Brier")
        
        cat("\nExpected Brier score at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
        return(ExpectedBrier)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    } else if (metric == "MAPE") {
      
      if(length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
        
        sum <- as.data.frame(matrix(0, nrow = length(object$horizon.time), ncol = 2))
        sum[, 1] <- object$horizon.time
        colnames(sum) <- c("Horizon Time", "SurvProb")
        for (i in 1:length(object$horizon.time)) {
          for (j in 1:object$n.cv) {
            sum[i, 2] <- sum[i, 2] + mean(abs(object$metric.cv[[j]]$AllSurv[[i]][, 1] - 
                                                object$metric.cv[[j]]$AllSurv[[i]][, 2])) 
          }
        }
        sum[, -1] <- sum[, -1]/object$n.cv
        
        sum[, -1] <- round(sum[, -1], digits)
        cat("\nMean of absolute error across quantiles of predicted risk scores at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
        return(sum)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    } else {
      
      if (length(object$metric.cv) == object$n.cv && sum(mapply(is.null, object$metric.cv)) == 0) {
        sum <- 0
        for (j in 1:object$n.cv) {
          sum <- sum + object$metric.cv[[j]]
        }
        sum <- sum/object$n.cv
        
        AUC <- round(sum[, 1], digits)
        ExpectedAUC <- data.frame(object$horizon.time, AUC)
        
        if (metric == "AUC") {
          colnames(ExpectedAUC) <- c("Horizon Time", "AUC")
          cat("\nExpected AUC at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
        } else {
          colnames(ExpectedAUC) <- c("Horizon Time", "Cindex")
          cat("\nExpected Cindex at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
        }
        return(ExpectedAUC)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    }
    
    
    
  }
  
  
  
  
}