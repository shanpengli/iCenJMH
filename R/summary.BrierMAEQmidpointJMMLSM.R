##' @export
##' 

summary.BrierMAEQmidpointJMMLSM <- function (object, digits = 3, ...) {
  if (!inherits(object, "BrierMAEQmidpointJMMLSM"))
    stop("Use only with 'BrierMAEQmidpointJMMLSM' xs.\n") 
  
  if (is.null(object$MAEQ.cv) && is.null(object$Brier.cv)) {
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
    if(length(object$Brier.cv) == object$n.cv && sum(mapply(is.null, object$Brier.cv)) == 0) {
      
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$Brier.cv[[j]]
      }
      sum <- sum/object$n.cv
      
      Brier <- sum[, 1]
      ExpectedBrier <- data.frame(object$horizon.time, Brier)
      colnames(ExpectedBrier) <- c("Horizon Time", "Brier")
      
      cat("\nExpected Brier score at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      return(ExpectedBrier)
      
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
  }
  
}