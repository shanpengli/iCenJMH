##' @title Print BrieriCenJMMLSM
##' @name summary.BrieriCenJMMLSM
##' @aliases summary.BrieriCenJMMLSM
##' @param object object of class 'BrieriCenJMMLSM'.
##' @param digits number of digits of decimal to be printed. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}, \link{survfitJMMLSM}}
##' @export
##' 

summary.BrieriCenJMMLSM <- function (object, digits = 4, ...) {
  if (!inherits(object, "BrieriCenJMMLSM"))
    stop("Use only with 'BrieriCenJMMLSM' xs.\n") 
  
  if (is.null(object$Brier.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    if (length(object$Brier.cv) == object$n.cv && sum(mapply(is.null, object$Brier.cv)) == 0) {
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