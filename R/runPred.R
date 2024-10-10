##' @export
##' 

runPred <- function(i = 1, seed = 100, metric, n, exact.obs = 0.7,
                    landmark.time, horizon.time, datatype = NULL, model = c("both", "iCenJMH", "JMH")) {
  
  library(iCenJMH)
  library(JMH)
  
  Pred <- Pred.JMMLSM <- NULL
  
  data <- SimJMdata(seed = seed, n = n,
                    nc = 100,
                    covbw = matrix(c(1, 0.1, 0.1, 0.5), nrow = 2, ncol = 2),
                    lambda = 0.3, lambdaC = 0.05,
                    Cmin = 0.1,
                    Cmax = 1,
                    gamma = c(-0.05, 0.2, -0.1),
                    alpha = c(0.5, -0.5),
                    beta = c(5, 2, -3, -3, 2),
                    tau = c(1, 0.05, 0.2, 0.1, -0.2),
                    increment = 0.3,
                    digits = 1, exact.obs = exact.obs)
  
  Ydata <- data$Ydata
  Sdata <- data$Sdata
  Tdata <- data$Tdata
  iCen.info <- GetSupport(iCen.data = Sdata, iCen.tL = "Li", iCen.tR = "Ri",
                          ID = "ID", S = "Stime", weight = "Stime.w", weight.ID = "w.ID")
  
  if (model == "both") {
    
    iCenfit <- iCenJMMLSM(Ydata = Ydata, Tdata = Tdata,
                          long.formula = Yij ~ X11 + X12,
                          surv.formula = Surv(Ei, status) ~ X21 + X22,
                          variance.formula = ~ X11 + X12,
                          random = ~ 1|ID,
                          timeVar = "Oij",
                          iCen.info = iCen.info,
                          maxiter = 1000, epsilon = 1e-04, print.para = FALSE,
                          initial.para = FALSE,
                          pStol = 1e-06, c = 0.95, quadpoint = 8,
                          hazard.kernel = "Epanechnikov", method = "adaptive")
    
    Pred <- DynPredAcciCenJMMLSM(seed = seed + i, object = iCenfit, landmark.time = landmark.time, 
                                 horizon.time = horizon.time,
                                 obs.time = "Oij", metric = metric, n.cv = 4)
    
    #sum <- summary(Pred)
    
    Pred.JMMLSM <- DynPredJMMLSM(seed = seed + i, 
                                 Ydata = Ydata,
                                 Tdata = Tdata,
                                 Sdata = Sdata,
                                 long.formula = Yij ~ X11 + X12, 
                                 surv.formula = Surv(Ei, status) ~ X21 + X22,
                                 variance.formula = ~ X11 + X12,
                                 random = ~ 1|ID,
                                 iCen.tL = "Li", iCen.tR = "Ri",
                                 landmark.time = landmark.time, 
                                 horizon.time = horizon.time, 
                                 obs.time = "Oij", 
                                 datatype = datatype,
                                 quadpoint = 8, maxiter = 1000, n.cv = 4, 
                                 quantile.width = 0.25,
                                 metric = metric)
    
    #sum.JMMLSM <- summary(Pred.JMMLSM)
    
  } else if (model == "iCenJMH") {
    
    iCenfit <- iCenJMMLSM(Ydata = Ydata, Tdata = Tdata,
                          long.formula = Yij ~ X11 + X12,
                          surv.formula = Surv(Ei, status) ~ X21 + X22,
                          variance.formula = ~ X11 + X12,
                          random = ~ 1|ID,
                          timeVar = "Oij",
                          iCen.info = iCen.info,
                          maxiter = 1000, epsilon = 1e-04, print.para = FALSE,
                          initial.para = FALSE,
                          pStol = 1e-06, c = 0.95, quadpoint = 8,
                          hazard.kernel = "Epanechnikov", method = "adaptive")
    
    Pred <- DynPredAcciCenJMMLSM(seed = seed + i, object = iCenfit, landmark.time = landmark.time, 
                                 horizon.time = horizon.time,
                                 obs.time = "Oij", metric = metric, n.cv = 4)
    
    #sum <- summary(Pred)
    
  } else {
    
    Pred.JMMLSM <- DynPredJMMLSM(seed = seed + i, 
                                 Ydata = Ydata,
                                 Tdata = Tdata,
                                 Sdata = Sdata,
                                 long.formula = Yij ~ X11 + X12, 
                                 surv.formula = Surv(Ei, status) ~ X21 + X22,
                                 variance.formula = ~ X11 + X12,
                                 random = ~ 1|ID,
                                 iCen.tL = "Li", iCen.tR = "Ri",
                                 landmark.time = landmark.time, 
                                 horizon.time = horizon.time, 
                                 obs.time = "Oij", 
                                 datatype = datatype,
                                 quadpoint = 8, maxiter = 1000, n.cv = 4, 
                                 quantile.width = 0.25,
                                 metric = metric)
    
    #sum.JMMLSM <- summary(Pred.JMMLSM)
    
  }
  

  return(list(Pred = Pred, Pred.JMMLSM = Pred.JMMLSM))
}