##' @title Dynamic prediction accuracy assessment for joint models
##' @name DynPredAcciCenJMMLSM
##' @aliases DynPredAcciCenJMMLSM
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'iCenJMMLSM'.
##' @param landmark.time a numeric value of time for which dynamic prediction starts..
##' @param horizon.time a numeric vector of future times for which predicted probabilities are to be computed.
##' @param obs.time a character string of specifying a longitudinal time variable.
##' @param quadpoint the number of standard Gauss-Hermite quadrature points.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param metric a string to indicate which metric is used. 
##' The available options for assessing the prediction accuracy are \code{Brier Score}, \code{AUC}, \code{Cindex}, \code{MAPE}.
##' @param quantile.width a numeric value of width of quantile to be specified if \code{metric} is specified as \code{MAPE}. Default is 0.25.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}, \link{survfitiCenJMMLSM}}
##' @export
##' 

DynPredAcciCenJMMLSM <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL, 
                            obs.time = NULL, quadpoint = NULL, maxiter = 1000, n.cv = 3, 
                            metric = c("Brier Score", "AUC", "Cindex", "MAPE"), quantile.width = NULL, ...) {
  
  
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
  if (metric == "MAPE") {
    if (is.null(quantile.width)) {
      quantile.width <- 0.25
    }
    groups <- 1/quantile.width
    if (floor(groups) != groups)
      stop("The reciprocal of quantile.width must be an integer.")
  }
  
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
  metric.cv <- list()
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
                          method = object$method, 
                          random = object$random,
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
      metric.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      metric.cv[[t]] <- NULL
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
                                       quadpoint = fit$quadpoint,
                                       Last.time = rep(landmark.time, nrow(val.Tdata)), obs.time = obs.time, 
                                       pStol = 1e-2), silent = TRUE)
      
      if ('try-error' %in% class(survfit)) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        metric.cv[[t]] <- NULL
      } else {
        
        if (metric == "Brier Score") {
          
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.Tdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.Tdata[, ID]
          
          ## fit a Kalplan-Meier estimator
          New.surv.formula.out <- paste0("survival::Surv(", surv.var[1], ",", 
                                         surv.var[2], "==0)")
          New.surv.formula <- as.formula(paste(New.surv.formula.out, 1, sep = "~"))
          fitKM <- survival::survfit(New.surv.formula, data = val.Tdata)
          
          Gs <- summary(fitKM, times = landmark.time)$surv
          mean.Brier <- matrix(NA, nrow = length(horizon.time), ncol = 1)
          for (j in 1:length(horizon.time)) {
            fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
            if ('try-error' %in% class(fitKM.horizon)) {
              mean.Brier[j, 1] <- NA
            } else {
              Gu <- fitKM.horizon$surv
              ## true counting process
              N1 <- vector()
              Gt <- vector()
              W.IPCW <- vector()
              for (i in 1:nrow(Surv)) {
                if (val.Tdata[i, surv.var[1]] <= horizon.time[j] && val.Tdata[i, surv.var[2]] == 1) {
                  N1[i] <- 1
                  Gt[i] <- summary(fitKM, times = val.Tdata[i, surv.var[1]])$surv
                } else {
                  N1[i] <- 0
                  Gt[i] <- NA
                }
                
                if (val.Tdata[i, surv.var[1]] > horizon.time[j]) {
                  W.IPCW[i] <- 1/(Gu/Gs)
                } else if (val.Tdata[i, surv.var[1]] <= horizon.time[j] && val.Tdata[i, surv.var[2]] == 1) {
                  W.IPCW[i] <- 1/(Gt[i]/Gs)
                } else {
                  W.IPCW[i] <- NA
                }
              }
              ## extract estimated Survival probability
              for (k in 1:nrow(Surv)) {
                Surv[k, 2] <- survfit$Pred[[k]][j, 2]
              }
              
              RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
              colnames(RAWData.Brier)[1:2] <- c("ID", "Surv")
              RAWData.Brier$Brier <- RAWData.Brier$W.IPCW*
                abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^2
              mean.Brier[j, 1] <- sum(RAWData.Brier$Brier, na.rm = TRUE)/nrow(RAWData.Brier)
            }
            
            
          }
          metric.cv[[t]] <- mean.Brier
          
        } else if (metric == "MAPE") {
          
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
          
          metric.cv[[t]] <- result
          
        } else if (metric %in% c("AUC", "Cindex")) {
          
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.Tdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.Tdata[, ID]
          Surv$time <- val.Tdata[, surv.var[1]]
          Surv$status <- val.Tdata[, surv.var[2]]
          mean.AUC <- matrix(NA, nrow = length(horizon.time), ncol = 1)
          ## extract estimated Survival probability
          for (j in 1:length(horizon.time)) {
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit$Pred[[k]][j, 2]
            }
            
            if (metric == "AUC") {
              ROC <- timeROC::timeROC(T = Surv$time, delta = Surv$status,
                                      weighting = "marginal",
                                      marker = -Surv$Surv, cause = 1,
                                      times = horizon.time[j])
              
              mean.AUC[j, 1] <- ROC$AUC[2]
            } else {
              
              mean.AUC[j, 1] <- CindexCR(Surv$time, Surv$status, -Surv$Surv, 1)
              
            }
            
            
          }
          
          metric.cv[[t]] <- mean.AUC
          
        } else {
          stop("Please choose one of the following metrics: Brier Score, AUC, Cindex, MAPE.")
        }

      }
    }
    writeLines(paste0("The ", t, " th validation is done!"))
  }
  result <- list(metric.cv = metric.cv, n.cv = n.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, quadpoint = quadpoint, 
                 seed = seed, metric = metric)
  class(result) <- "DynPredAcciCenJMMLSM"
  return(result)
}