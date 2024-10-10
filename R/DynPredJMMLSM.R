##' @export
##' 

DynPredJMMLSM <- function(seed = 100, 
                               Ydata = NULL,
                               Tdata = NULL,
                               Sdata = NULL,
                               long.formula = NULL, 
                               surv.formula = NULL,
                               variance.formula = NULL,
                               random = NULL,
                               iCen.tL = NULL, iCen.tR = NULL,
                               landmark.time = NULL, horizon.time = NULL, 
                               obs.time = NULL, 
                               int.time.Var = NULL,
                               datatype = c("midpoint", "rightpoint", "uniform"),
                               quadpoint = NULL, maxiter = 1000, n.cv = 3, 
                               quantile.width = 0.25, 
                               metric = c("MAPE", "Brier Score", "AUC", "Cindex"), ...) {
  
  if (is.null(landmark.time)) 
    stop("Please specify the landmark.time for dynamic prediction.")   
  if (!is.vector(horizon.time)) 
    stop("horizon.time must be vector typed.")
  if (!is.null(int.time.Var)) {
    if (!int.time.Var %in% colnames(Ydata)) {
      stop(paste(int.time.Var, "does not exist in Ydata."))
    }
  }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from Ydata.")
  } else {
    if (!obs.time %in% colnames(Ydata)) {
      stop(paste0(obs.time, " is not found in Ydata."))
    }
  }
  groups <- 1/quantile.width
  if (floor(groups) != groups)
    stop("The reciprocal of quantile.width must be an integer.")
  set.seed(seed)

  long.var <- all.vars(long.formula)
  variance.var <- all.vars(variance.formula)
  surv.var <- all.vars(surv.formula)
  random.var <- all.vars(random) 
  ID <- random.var[length(random.var)]
  
  folds <- caret::groupKFold(c(1:nrow(Tdata)), k = n.cv)
  metric.cv <- list()
  data <- Sdata[, c(ID, iCen.tL, iCen.tR)]
  if (datatype == "midpoint") {
    data$Stime <- (data[, iCen.tL] + data[, iCen.tR])/2
  } else if (datatype == "uniform") {
    f <- Vectorize(function(iCen.tL, iCen.tR) runif(1, min = iCen.tL, max = iCen.tR), 
                   vectorize.args = c("iCen.tL", "iCen.tR"))
    
    Stime <- with(data, f(iCen.tL, iCen.tR))
    data <- data.frame(data, Stime)
    colnames(data) <- c(ID, iCen.tL, iCen.tR, Stime)
  } else if (datatype == "rightpoint") {
    data$Stime <- data[, iCen.tR]
  } else {
    stop("Please choose the right datatype.")
  }
  for (t in 1:n.cv) {
    
    train.Tdata <- Tdata[folds[[t]], ]
    train.Ydata <- Ydata[Ydata[, ID] %in% train.Tdata[, ID], ]
    train.Sdata <- Sdata[Sdata[, ID] %in% train.Tdata[, ID], ]
    
    train.Ydata <- dplyr::left_join(train.Ydata, data, by = ID)
    train.Ydata$time <- train.Ydata[, obs.time] - train.Ydata$Stime
    train.Tdata <- dplyr::left_join(train.Tdata, data, by = ID)
    train.Tdata$survtime <- train.Tdata[, surv.var[1]] - train.Tdata$Stime
    
    if (!is.null(int.time.Var)) {
      interaction <- paste("time", int.time.Var, sep = ":")
      long.fixed <- paste(c("Stime", "time", long.var[2:length(long.var)], interaction), collapse = "+")
    } else {
      long.fixed <- paste(c("Stime", "time", long.var[2:length(long.var)]), collapse = "+")
    }
    
    new.long.formula <- as.formula(paste(long.var[1], long.fixed, sep = "~"))
    
    if (!is.null(int.time.Var)) {
      interaction <- paste("time", int.time.Var, sep = ":")
      long.fixed <- paste(c("Stime", "time", variance.var, interaction), collapse = "+")
    } else {
      long.fixed <- paste(c("Stime", "time", variance.var), collapse = "+")
    }
    
    new.var.formula <- as.formula(paste0("~", long.fixed))
    
    surv.fixed <- paste(c("Stime", surv.var[3:length(surv.var)]), collapse = "+")
    surv.out <- paste0("survival::Surv(", "survtime", ", ", surv.var[2], ")")
    new.surv.formula <- as.formula(paste(surv.out, surv.fixed, sep = "~"))  
    
    train.Ydata <- train.Ydata[, c(ID, long.var[1], "Stime", "time", long.var[2:length(long.var)])]
    train.Tdata <- train.Tdata[, c(ID, surv.var[1:2], "Stime", surv.var[3:length(surv.var)])]
    colnames(train.Tdata)[2] <- "survtime"
    
    a <- proc.time()
    fit <- try(JMH::JMMLSM(cdata = train.Tdata, ydata = train.Ydata,
               long.formula = new.long.formula,
               variance.formula = new.var.formula,
               surv.formula = new.surv.formula,
               quadpoint = quadpoint,
               maxiter = maxiter,
               opt = "optim",
               random = random), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      
      metric.cv[[t]] <- NULL

    } else if (fit$iter == maxiter) {
      metric.cv[[t]] <- NULL

    } else {
      writeLines(paste0("The ", t, " th training is done!"))
      val.Tdata <- Tdata[-folds[[t]], ]
      val.Ydata <- Ydata[Ydata[, ID] %in% val.Tdata[, ID], ]
      val.Sdata <- Sdata[Sdata[, ID] %in% val.Tdata[, ID], ]
      ## Only consider the subjects who receive the first positive diagnosis of the initial event at the landmark time
      val.Sdata <- val.Sdata[val.Sdata[, iCen.tR] <= landmark.time, ]
      val.Tdata <- val.Tdata[val.Tdata[, ID] %in% val.Sdata[, ID], ]
      val.Ydata <- val.Ydata[val.Ydata[, ID] %in% val.Sdata[, ID], ]

      val.Tdata <- val.Tdata[val.Tdata[, surv.var[1]] > landmark.time, ]
      val.Ydata <- val.Ydata[val.Ydata[, ID] %in% val.Tdata[, ID], ]
      val.Ydata <- val.Ydata[val.Ydata[, obs.time] <= landmark.time, ]
      NewyID <- unique(val.Ydata[, ID])
      val.Tdata <- val.Tdata[val.Tdata[, ID] %in% NewyID, ]
      val.Sdata <- val.Sdata[val.Sdata[, ID] %in% NewyID, ]
      
      val.Ydata <- dplyr::left_join(val.Ydata, data, by = ID)
      val.Ydata$time <- val.Ydata[, obs.time] - val.Ydata$Stime
      val.Ydata <- val.Ydata[, c(ID, long.var[1], "Stime", "time", long.var[2:length(long.var)])]
      val.Tdata <- dplyr::left_join(val.Tdata, data, by = ID)
      val.Tdata <- val.Tdata[, c(ID, surv.var[1:2], "Stime", surv.var[3:length(surv.var)])]
      colnames(val.Tdata)[2] <- "survtime"
      
      survfit <- list()
      errormess <- vector()
      for (i in 1:length(NewyID)) {
        ynewdata <- val.Ydata[val.Ydata[, ID] == NewyID[i], ]
        cnewdata <- val.Tdata[val.Tdata[, ID] == NewyID[i], ]
        u <- horizon.time - cnewdata$Stime
        Last.time <- landmark.time - cnewdata$Stime
        survfit[[i]] <- try(JMH::survfitJMMLSM(fit, seed = seed, ynewdata = ynewdata, 
                                               cnewdata = cnewdata, u = u, 
                                               Last.time = Last.time, obs.time = "time", method = "GH"), silent = TRUE)
        
        if ('try-error' %in% class(survfit[[i]])) {
          errormess[i] <- 0
        } else {
            next
        }
      }
      
      if (prod(errormess) == 0) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        metric.cv[[t]] <- NULL
      } else {
        
        if (metric == "MAPE") {
          AllSurv <- list()
          for (j in 1:length(horizon.time)) {
            Surv <- as.data.frame(matrix(0, nrow = nrow(val.Tdata), ncol = 2))
            colnames(Surv) <- c("ID", "Surv")
            Surv$ID <- val.Tdata[, ID]
            ## extract estimated survival prob
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit[[k]]$Pred[[1]][j, 2]
            }
            Surv <- Surv[!is.nan(Surv$Surv), ]
            ## group subjects based on survival prob
            quant <- quantile(Surv$Surv, probs = seq(0, 1, by = quantile.width))
            EmpiricalSurv <- rep(NA, groups)
            PredictedSurv <- rep(NA, groups)
            for (i in 1:groups) {
              subquant <- Surv[Surv$Surv > quant[i] &
                                 Surv$Surv <= quant[i+1], c(1, 2)]
              quantsubdata <- val.Tdata[val.Tdata[, ID] %in% subquant$ID, 2:3]
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
        } else if (metric == "Brier Score") {
          
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.Tdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.Tdata[, ID]
          
          ## fit a Kalplan-Meier estimator
          New.surv.formula.out <- paste0("survival::Surv(survtime,", 
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
                if (val.Tdata[i, 2] <= horizon.time[j] && val.Tdata[i, 3] == 1) {
                  N1[i] <- 1
                  Gt[i] <- summary(fitKM, times = val.Tdata[i, 2])$surv
                } else {
                  N1[i] <- 0
                  Gt[i] <- NA
                }
                
                if (val.Tdata[i, 2] > horizon.time[j]) {
                  W.IPCW[i] <- 1/(Gu/Gs)
                } else if (val.Tdata[i, 2] <= horizon.time[j] && val.Tdata[i, 3] == 1) {
                  W.IPCW[i] <- 1/(Gt[i]/Gs)
                } else {
                  W.IPCW[i] <- NA
                }
              }
              ## extract estimated Survival probability
              for (k in 1:nrow(Surv)) {
                Surv[k, 2] <- survfit[[k]]$Pred[[1]][j, 2]
              }
              
              RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
              colnames(RAWData.Brier)[1:2] <- c("ID", "Surv")
              RAWData.Brier$Brier <- RAWData.Brier$W.IPCW*
                abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^2
              mean.Brier[j, 1] <- sum(RAWData.Brier$Brier, na.rm = TRUE)/nrow(RAWData.Brier)
            }
            
            
          }
          metric.cv[[t]] <- mean.Brier
          
        } else if (metric %in% c("AUC", "Cindex")) {
          
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.Tdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.Tdata[, ID]
          Surv$time <- val.Tdata[, "survtime"]
          Surv$status <- val.Tdata[, surv.var[2]]
          mean.AUC <- matrix(NA, nrow = length(horizon.time), ncol = 1)
          ## extract estimated Survival probability
          for (j in 1:length(horizon.time)) {
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit[[k]]$Pred[[1]][j, 2]
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
          stop("Please specify one of the following metrics for prediction assessment: MAPE, Brier Score, AUC.")
        }


        writeLines(paste0("The ", t, " th validation is done!"))
      }
    }
  }
  
  result <- list(metric.cv = metric.cv, n.cv = n.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, quadpoint = quadpoint, 
                 seed = seed, metric = metric)
  class(result) <- "DynPredJMMLSM"
  return(result)
}
