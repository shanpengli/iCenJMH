##' @export
##' 

bootsDynPrediCenJMMLSM <- function(i = 100, seed = 100, n = 2000, landmark.time = 4,
                                horizon.time = 5, increment = 0.3, exact.obs = 0,
                                n.cv = 4, 
                                datatype = c("midpoint", "rightpoint", "uniform"),
                                metric = c("Brier Score", "AUC", "MAPE")) {
  
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
                    increment = increment,
                    digits = 1, exact.obs = exact.obs)
  
  Ydata <- data$Ydata
  Sdata <- data$Sdata
  Tdata <- data$Tdata
  iCen.info <- GetSupport(iCen.data = Sdata, iCen.tL = "Li", iCen.tR = "Ri",
                          ID = "ID", S = "Stime", weight = "Stime.w", weight.ID = "w.ID")

  iCenfit <- iCenJMMLSM(Ydata = Ydata, Tdata = Tdata,
                        long.formula = Yij ~ X11 + X12,
                        surv.formula = Surv(Ei, status) ~ X21 + X22,
                        variance.formula = ~ X11 + X12,
                        random = ~ 1|ID,
                        timeVar = "Oij",
                        iCen.info = iCen.info,
                        maxiter = 1000, epsilon = 1e-04,
                        quadpoint = 8, print.para = FALSE,
                        initial.para = TRUE,
                        pStol = 1e-06, c = 0.95,
                        method = "adaptive",
                        hazard.kernel = "Epanechnikov")
  
  MAEQ <- DynPredAcciCenJMMLSM(seed = seed + i, iCenfit, landmark.time = landmark.time, 
                         horizon.time = horizon.time, obs.time = "Oij", metric = metric,
                         n.cv = n.cv)
  iCenRES <- summary(MAEQ, digits = 4)
  
  JMMLSM.fit <- DynPredJMMLSM(seed = seed, 
                                Ydata = Ydata,
                                Tdata = Tdata,
                                Sdata = Sdata,
                                long.formula = Yij ~ X11 + X12, 
                                surv.formula = Surv(Ei, status) ~ X21 + X22,
                                variance.formula = ~ X11 + X12,
                                random = ~ 1|ID,
                                iCen.tL = "Li", iCen.tR = "Ri",
                                landmark.time = landmark.time, horizon.time = horizon.time, 
                                obs.time = "Oij", 
                                datatype = datatype,
                                quadpoint = 8, maxiter = 1000, n.cv = n.cv, 
                                quantile.width = 0.25,
                                metric = metric)
  
  JMMLSMRES <- summary(JMMLSM.fit, digits = 4) 
  
  result <- list(iCenRES = iCenRES, JMMLSMRES = JMMLSMRES)
  
}