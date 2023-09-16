getdummy <- function(long.formula, surv.formula, variance.formula, random, Ydata, Tdata) {
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  variance <- all.vars(variance.formula)
  cnames <- colnames(Tdata)
  ynames <- colnames(Ydata)
  ID <- all.vars(random)[length(all.vars(random))]
  
  random.var <- all.vars(random)
  m <- model.frame(long.formula, Ydata)
  Ydata2 <- model.matrix(long.formula, m)
  ynames <- colnames(Ydata2)
  Ydata2 <- data.frame(Ydata[, random.var[length(random.var)]], m[[1]], Ydata2)
  colnames(Ydata2) <- c(random.var[length(random.var)], names(m)[1], ynames)
  m <- model.frame(variance.formula, Ydata)
  Ydata3 <- model.matrix(variance.formula, m)
  ynames <- colnames(Ydata3)
  Ydata3 <- data.frame(Ydata[, random.var[length(random.var)]], Ydata3)
  colnames(Ydata3) <- c(random.var[length(random.var)], ynames)
  
  surv.formula <- surv.formula[-2]
  m <- model.frame(surv.formula, Tdata)
  Tdata2 <- model.matrix(surv.formula, m)
  cnames <- colnames(Tdata2)
  surv.var <- survival[c(1:2)]
  Tdata <- Tdata[, c(ID, surv.var)]
  Tdata2 <- data.frame(Tdata, Tdata2[, -1])
  colnames(Tdata2) <- c(ID, surv.var, cnames[-1])
  
  result <- list(Ydata2, Ydata3, Tdata2)
  names(result) <- c("Ydata.mean", "Ydata.variance", "Tdata")
  return(result)
  
}