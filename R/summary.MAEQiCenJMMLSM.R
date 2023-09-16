##' @title Print MAEQiCenJMMLSM
##' @name summary.MAEQiCenJMMLSM
##' @aliases summary.MAEQiCenJMMLSM
##' @param object object of class 'MAEQiCenJMMLSM'.
##' @param digits number of decimal points to be rounded.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
##' @export
##' 

summary.MAEQiCenJMMLSM <- function (object, digits = 3, ...) {
  if (!inherits(object, "MAEQiCenJMMLSM"))
    stop("Use only with 'MAEQiCenJMMLSM' xs.\n") 
  
  if (is.null(object$MAEQ.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    if(length(object$MAEQ.cv) == object$n.cv && sum(mapply(is.null, object$MAEQ.cv)) == 0) {
      
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
  }
  
}