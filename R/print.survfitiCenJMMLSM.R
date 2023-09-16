##' @title Print survfitiCenJMMLSM
##' @name print.survfitiCenJMMLSM
##' @aliases print.survfitiCenJMMLSM
##' @param x x of class 'survfitiCenJMMLSM'.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
##' @export
##' 
print.survfitiCenJMMLSM <- function (x, ...) {
  if (!inherits(x, "survfitiCenJMMLSM"))
    stop("Use only with 'survfitiCenJMMLSM' xs.\n")
  
  f <- function (d, t) {
    a <- matrix(1, nrow = 1, ncol = 2)
    a[1, 1] <- t 
    a <- as.data.frame(a)
    colnames(a) <- colnames(d)
    d <- rbind(a, d)
    d
  }
  
  cat("\nPrediction of Conditional Probabilities of Event\nbased on the standard Guass-Hermite quadrature rule with", x$quadpoint,
      "quadrature points\n")
  
  print(mapply(f, x$Pred, x$Last.time[, 2], SIMPLIFY = FALSE))
  
  invisible(x)
  
}