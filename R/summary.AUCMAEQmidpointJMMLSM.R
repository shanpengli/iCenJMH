##' @export
##' 

summary.AUCMAEQmidpointJMMLSM <- function (object, digits = 3, ...) {
  if (!inherits(object, "AUCMAEQmidpointJMMLSM"))
    stop("Use only with 'AUCMAEQmidpointJMMLSM' xs.\n") 
  
  if (is.null(object$MAEQ.cv) && is.null(object$AUC.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else if (!is.null(object$MAEQ.cv) && length(object$MAEQ.cv) == object$n.cv) {
    if(sum(mapply(is.null, object$MAEQ.cv)) == 0) {
      
      sum <- as.data.frame(matrix(0, nrow = length(object$horizon.time), ncol = 2))
      sum[, 1] <- object$horizon.time
      colnames(sum) <- c("Horizon Time", "SurvProb")
      for (i in 1:length(object$horizon.time)) {
        for (j in 1:object$n.cv) {
          sum[i, 2] <- sum[i, 2] + sum(abs(object$MAEQ.cv[[j]]$AllSurv[[i]][, 1] - 
                                             object$MAEQ.cv[[j]]$AllSurv[[i]][, 2])) 
        }
      }
      sum[, -1] <- sum[, -1]/object$n.cv
      
      sum[, -1] <- round(sum[, -1], digits)
      cat("\nSum of absolute error across quantiles of predicted risk scores at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      return(sum)
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
  } else {
    if(length(object$AUC.cv) == object$n.cv && sum(mapply(is.null, object$AUC.cv)) == 0) {
      
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$AUC.cv[[j]]
      }
      sum <- sum/object$n.cv
      
      AUC <- sum[, 1]
      ExpectedAUC <- data.frame(object$horizon.time, AUC)
      colnames(ExpectedAUC) <- c("Horizon Time", "AUC")
      
      cat("\nExpected AUC at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      return(ExpectedAUC)
      
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
  }
  
}