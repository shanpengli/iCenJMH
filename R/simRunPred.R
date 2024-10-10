##' @export
##' 

simRunPred <- function(sim = 10, n = 3000, seed = 100,
                       exact.obs = 0.7,
                       landmark.time = 6, horizon.time = c(7:12), 
                       datatype = NULL, model = c("both", "iCenJMH", "JMH"), mc.cores = 10) {
  
  if (is.null(datatype) & model != "iCenJMH") {
    stop("A datatype must be spercified.")
  }
  
  boots.BS <- boots.AUC <- boots.MAPE <- boots.Cindex <- list()
  for (i in 1:sim) {
    
    boots.BS[[i]] <- mclapply(1:sim, runPred,
                              seed = seed + i, metric = "Brier Score", 
                              n = n, exact.obs = exact.obs,
                              landmark.time = landmark.time, horizon.time = horizon.time, 
                              datatype = datatype, model = model, mc.cores = mc.cores)
    
    boots.AUC[[i]] <- mclapply(1:sim, runPred,
                               seed = seed + i, metric = "AUC", 
                               n = n, exact.obs = exact.obs,
                               landmark.time = landmark.time, horizon.time = horizon.time, 
                               datatype = datatype, model = model, mc.cores = mc.cores)
    
    boots.MAPE[[i]] <- mclapply(1:sim, runPred,
                                seed = seed + i, metric = "MAPE", 
                                n = n, exact.obs = exact.obs,
                                landmark.time = landmark.time, horizon.time = horizon.time, 
                                datatype = datatype, model = model, mc.cores = mc.cores)
    
    boots.Cindex[[i]] <- mclapply(1:sim, runPred,
                                  seed = seed + i, metric = "Cindex", 
                                  n = n, exact.obs = exact.obs,
                                  landmark.time = landmark.time, horizon.time = horizon.time, 
                                  datatype = datatype, model = model, mc.cores = mc.cores)
    
  }
  
  return(list(boots.BS = boots.BS, boots.AUC = boots.AUC,
              boots.MAPE = boots.MAPE, boots.Cindex = boots.Cindex,
              datatype = datatype, model = model))
  
}