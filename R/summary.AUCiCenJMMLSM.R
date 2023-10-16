##' @title Print AUCiCenJMMLSM
##' @name summary.AUCiCenJMMLSM
##' @aliases summary.AUCiCenJMMLSM
##' @param object object of class 'AUCiCenJMMLSM'.
##' @param digits number of digits of decimal to be printed. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}, \link{survfitJMMLSM}}
##' @export
##' 

summary.AUCiCenJMMLSM <- function (object, digits = 4, ...) {
  if (!inherits(object, "AUCiCenJMMLSM"))
    stop("Use only with 'AUCiCenJMMLSM' xs.\n") 
  
  if (is.null(object$AUC.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    if (length(object$AUC.cv) == object$n.cv && sum(mapply(is.null, object$AUC.cv)) == 0) {
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