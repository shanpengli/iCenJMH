##' @title A metric of prediction accuracy of joint model by comparing the predicted risk
##' with the empirical risks stratified on different predicted risk group.
##' @name MAEQiCenJMMLSM
##' @aliases MAEQiCenJMMLSM
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'iCenJMMLSM'.
##' @param landmark.time a numeric value of time for which dynamic prediction starts..
##' @param horizon.time a numeric vector of future times for which predicted probabilities are to be computed.
##' @param obs.time a character string of specifying a longitudinal time variable.
##' @param quadpoint the number of standard Gauss-Hermite quadrature points.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param quantile.width a numeric value of width of quantile to be specified. Default is 0.25.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
##' @export
##' 

MAEQiCenJMMLSM <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL, 
                           obs.time = NULL, quadpoint = NULL, maxiter = 1000, n.cv = 3, 
                           quantile.width = 0.25, ...) {
  
  if (!inherits(object, "iCenJMMLSM"))
    stop("Use only with 'iCenJMMLSM' xs.\n")
  if (is.null(landmark.time)) 
    stop("Please specify the landmark.time for dynamic prediction.")   
  if (!is.vector(horizon.time)) 
    stop("horizon.time must be vector typed.")
  if (is.null(quadpoint)) {
    quadpoint <- object$quadpoint
  }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from Ydata.")
  } else {
    if (!obs.time %in% colnames(object$Ydata)) {
      stop(paste0(obs.time, " is not found in Ydata."))
    }
  }
  groups <- 1/quantile.width
  if (floor(groups) != groups)
    stop("The reciprocal of quantile.width must be an integer.")
  set.seed(seed)
  Tdata <- object$Tdata
  Ydata <- object$Ydata
  long.formula <- object$long.formula
  surv.formula <- object$surv.formula
  surv.var <- all.vars(surv.formula)
  variance.formula <- object$variance.formula
  random <- all.vars(object$random) 
  ID <- random[length(random)]
  iCen.info <- object$iCen.info
  Sdata <- unique(iCen.info$iCen.data[, c(ID, iCen.info$iCen.tL, iCen.info$iCen.tR)])
  int.time.Var <- object$int.time.Var
  
  folds <- caret::groupKFold(c(1:nrow(Tdata)), k = n.cv)
  MAEQ.cv <- list()
  for (t in 1:n.cv) {
    
    train.Tdata <- Tdata[folds[[t]], ]
    train.Ydata <- Ydata[Ydata[, ID] %in% train.Tdata[, ID], ]
    train.Sdata <- Sdata[Sdata[, ID] %in% train.Tdata[, ID], ]
    
    train.iCen.info <- GetSupport(iCen.data = train.Sdata, iCen.tL = iCen.info$iCen.tL, iCen.tR = iCen.info$iCen.tR,
                                  ID = ID, S = iCen.info$S, weight = iCen.info$weight, weight.ID = iCen.info$weight.ID)
    
    fit <- try(iCenJMMLSM(Tdata = train.Tdata, Ydata = train.Ydata, 
                          long.formula = long.formula,
                          surv.formula = surv.formula,
                          variance.formula = variance.formula, 
                          timeVar = obs.time,
                          quadpoint = quadpoint, random = object$random,
                          iCen.info = train.iCen.info,
                          maxiter = maxiter,
                          epsilon = object$epsilon,
                          print.para = FALSE,
                          initial.para = FALSE,
                          pStol = object$pStol,
                          c = object$c,
                          hazard.kernel = object$hazard.kernel,
                          int.time.Var = int.time.Var), silent = TRUE)
    
    writeLines(paste0("The ", t, " th training is done!"))
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      MAEQ.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      MAEQ.cv[[t]] <- NULL
    } else {

      val.Tdata <- Tdata[-folds[[t]], ]
      val.Ydata <- Ydata[Ydata[, ID] %in% val.Tdata[, ID], ]
      val.Sdata <- Sdata[Sdata[, ID] %in% val.Tdata[, ID], ]
      ## Only consider the subjects who receive the first positive diagnosis of the initial event at the landmark time
      val.Sdata <- val.Sdata[val.Sdata[, iCen.info$iCen.tR] <= landmark.time, ]
      val.Tdata <- val.Tdata[val.Tdata[, ID] %in% val.Sdata[, ID], ]
      val.Ydata <- val.Ydata[val.Ydata[, ID] %in% val.Sdata[, ID], ]

      val.Tdata <- val.Tdata[val.Tdata[, surv.var[1]] > landmark.time, ]
      val.Ydata <- val.Ydata[val.Ydata[, ID] %in% val.Tdata[, ID], ]
      val.Ydata <- val.Ydata[val.Ydata[, obs.time] <= landmark.time, ]
      NewyID <- unique(val.Ydata[, ID])
      val.Tdata <- val.Tdata[val.Tdata[, ID] %in% NewyID, ]
      val.Sdata <- val.Sdata[val.Sdata[, ID] %in% NewyID, ]
      
      
      survfit <- try(survfitiCenJMMLSM(fit, seed = seed, Ynewdata = val.Ydata, 
                                       Tnewdata = val.Tdata, 
                                       iCennewdata = val.Sdata,
                                       u = horizon.time, 
                                       Last.time = rep(landmark.time, nrow(val.Tdata)), obs.time = obs.time, 
                                       pStol = 1e-2), silent = TRUE)

      if ('try-error' %in% class(survfit)) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        MAEQ.cv[[t]] <- NULL
      } else {

        AllSurv <- list()

        for (j in 1:length(horizon.time)) {
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.Tdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.Tdata[, ID]
          ## extract estimated survival prob
          for (k in 1:nrow(Surv)) {
            Surv[k, 2] <- survfit$Pred[[k]][j, 2]
          }
          Surv <- Surv[!is.nan(Surv$Surv), ]
          ## group subjects based on survival prob
          quant <- quantile(Surv$Surv, probs = seq(0, 1, by = quantile.width))
          EmpiricalSurv <- rep(NA, groups)
          PredictedSurv <- rep(NA, groups)
          for (i in 1:groups) {
            subquant <- Surv[Surv$Surv > quant[i] &
                               Surv$Surv <= quant[i+1], c(1, 2)]
            quantsubdata <- val.Tdata[val.Tdata[, ID] %in% subquant$ID, surv.var]
            colnames(quantsubdata) <- c("time", "status")
            fitKM <- survival::survfit(survival::Surv(time, status) ~ 1, data = quantsubdata)
            fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
            if ('try-error' %in% class(fitKM.horizon)) {
              EmpiricalSurv[i] <- summary(fitKM, times = max(quantsubdata$time))$surv
            } else {
              EmpiricalSurv[i] <- summary(fitKM, times = horizon.time[j])$surv
            }
            PredictedSurv[i] <-mean(subquant$Surv)
          }
          AllSurv[[j]] <- data.frame(EmpiricalSurv, PredictedSurv)
        }
        names(AllSurv) <- horizon.time
        result <- list(AllSurv = AllSurv)

        MAEQ.cv[[t]] <- result
        writeLines(paste0("The ", t, " th validation is done!"))
      }
    }
  }
  result <- list(MAEQ.cv = MAEQ.cv, n.cv = n.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, quadpoint = quadpoint, 
                 seed = seed)
  class(result) <- "MAEQiCenJMMLSM"
  return(result)
}
