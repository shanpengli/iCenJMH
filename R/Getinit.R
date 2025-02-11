Getinit <- function(Tdata = Tdata, Ydata = Ydata, long.formula = long.formula,
                    surv.formula = surv.formula, variance.formula = variance.formula,
                    model = model, RE = RE, random = random, timeVar = timeVar,
                    int.time.Var = int.time.Var, initial.para = initial.para, iCen.info = iCen.info) {
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  variance <- all.vars(variance.formula)
  cnames <- colnames(Tdata)
  ynames <- colnames(Ydata)
  ID <- iCen.info$ID
  yID <- unique(Ydata[, ID])
  cID <- Tdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in Ydata doesn't match with Tdata.")
  }
  
  ni <- as.data.frame(table(Ydata[, ID]))
  
  Impute.Sdata <- iCen.info$iCen.midpoint.data[, c(ID, iCen.info$S)]
  YdataS <- dplyr::left_join(Ydata, Impute.Sdata, by = ID)
  TdataS <- dplyr::left_join(Tdata, Impute.Sdata, by = ID)
    
  ## obtain initial guess for the parameters in the longitudinal sub-model
  YdataS$Ytime.aft.S <- unlist(YdataS[, timeVar] - YdataS[, iCen.info$S])
  longVar <- long[-1]
  if (!is.null(int.time.Var)) {
    YdataS$int.time.Var <- YdataS$Ytime.aft.S*YdataS[, int.time.Var]
    colnames(YdataS)[ncol(YdataS)] <- paste("Ytime.aft.S", int.time.Var, sep = ":")
    long.init.formula <- as.formula(paste(long.formula[2], 
                                          paste(c(iCen.info$S, "Ytime.aft.S", longVar,
                                                  colnames(YdataS)[ncol(YdataS)]), 
                                                collapse = "+"), sep = "~"))
  } else {
    long.init.formula <- as.formula(paste(long.formula[2], 
                                          paste(c(iCen.info$S, "Ytime.aft.S", longVar),
                                                collapse = "+"), sep = "~"))
  }

  longfit <- nlme::lme(fixed = long.init.formula, random = random, data = YdataS, method = "REML", 
                        control = nlme::lmeControl(opt = "optim"))
    
  ## obtain initial guess for the parameters in the survival sub-model
  TdataS <- as.data.frame(TdataS)
  TdataS$T.aft.S <- TdataS[, survival[1]] - TdataS[, iCen.info$S]
  survfmla.fixed <- c(iCen.info$S, survival[-c(1:2)])
  survfmla.fixed <- paste(survfmla.fixed, collapse = "+")
  survfmla.out1 <- paste0("survival::Surv(", "T.aft.S", ", ", survival[2], ")")
  survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
  fitSURV1 <- survival::coxph(formula = survfmla, data = TdataS, x = TRUE)
    
  if (!initial.para) {
    gamma <- as.vector(fitSURV1$coefficients)
    
    beta <- longfit$coefficients$fixed
    D <- as.matrix(nlme::getVarCov(longfit))
    resid <- resid(longfit)
    logResidsquare <- as.vector(log(resid^2))
    
    if (!is.null(int.time.Var)) {
      var.init.formula <- c(iCen.info$S, "Ytime.aft.S", variance, colnames(YdataS)[ncol(YdataS)])
    } else {
      var.init.formula <- c(iCen.info$S, "Ytime.aft.S", variance)
    }

    W <- as.matrix(data.frame(1, YdataS[, var.init.formula]))
    
    Tau <- OLS(W, logResidsquare)
    tau <- Tau$betahat
    
    if (model == "intercept") {
      alpha <- rep(0, 2)
      Sig <- matrix(0, nrow = 2, ncol = 2)
      Sig[1, 1] <- D
      Sig[2, 2] <- 1
      Sig[1, 2] <- 0
      Sig[2, 1] <- Sig[1, 2]
    } else {
      alpha = rep(0, 3)
      Sig <- matrix(0, nrow = 3, ncol = 3)
      Sig[1:2, 1:2] <- D
      Sig[3, 3] <- 1
      Sig[1:2, 3] <- 0
      Sig[3, 1:2] <- Sig[1:2, 3] 
    }
  } else {
    beta = c(5, 2, -3, -3, 2)
    tau = c(1, 0.05, 0.2, 0.1, -0.2)
    gamma = c(-0.05, 0.2, -0.1)
    alpha = c(0.5, -0.5)
    Sig = matrix(c(1, 0.1, 0.1, 0.5), nrow = 2, ncol = 2)
  }

  fixed.para <- list(beta = beta, tau = tau, gamma = gamma, alpha = alpha, Sig = Sig)
  
  ## obtain covariates of Tdata
  iCen.data <- iCen.info$iCen.data
  Tdata <- dplyr::left_join(Tdata, iCen.data[, c(ID, iCen.info$S, iCen.info$weight, 
                                                 iCen.info$weight.ID)], by = ID)
  Tdata$T.aft.S <- Tdata[, survival[1]] - Tdata[, iCen.info$S]
  Tdata <- Tdata[, c(ID, "T.aft.S", survival[2], iCen.info$S, survival[3:length(survival)], 
                     iCen.info$weight, iCen.info$weight.ID)]
  Tdata <- Tdata[order(-Tdata$T.aft.S), ]
  survtime <- as.vector(Tdata[, "T.aft.S"])
  status <- as.vector(Tdata[, survival[2]])
  X2 <- as.matrix(Tdata[, c(4:(3+length(gamma)))])
  sortTIDwID <- Tdata[, c(1, ncol(Tdata))]
  
  ## obtain covariates of Ydata
  Ydata <- dplyr::left_join(Ydata, iCen.data[, c(ID, iCen.info$S, iCen.info$weight, 
                                                 iCen.info$weight.ID)], by = ID)
  Ydata <- Ydata[order(Ydata[, ID], Ydata[, iCen.info$weight.ID]), ]
  Ydata$Ytime.aft.S <- Ydata[, timeVar] - Ydata[, iCen.info$S]
  if (!is.null(int.time.Var)) {
    Ydata$int.time.Var <- Ydata$Ytime.aft.S*Ydata[, int.time.Var]
    colnames(Ydata)[ncol(Ydata)] <- paste("Ytime.aft.S", int.time.Var, sep = ":")
    X1 <- as.matrix(cbind(1, Ydata[, c(iCen.info$S, "Ytime.aft.S", longVar, colnames(Ydata)[ncol(Ydata)])]))
    W <- as.matrix(cbind(1, Ydata[, c(iCen.info$S, "Ytime.aft.S", variance, colnames(Ydata)[ncol(Ydata)])]))
    WSfmla <- c(iCen.info$S, "Ytime.aft.S", variance, colnames(Ydata)[ncol(Ydata)]) 
  } else {
    X1 <- as.matrix(cbind(1, Ydata[, c(iCen.info$S, "Ytime.aft.S", longVar)]))
    W <- as.matrix(cbind(1, Ydata[, c(iCen.info$S, "Ytime.aft.S", variance)]))
    WSfmla <- c(iCen.info$S, "Ytime.aft.S", variance) 
  }

  WSfmla <- paste(WSfmla, collapse = "+")
  variance.formula <- as.formula(paste("log(sigma^2)", WSfmla, sep = "~"))
  
  random.var <- all.vars(random)
  ##random effect covariates
  if (model == "interslope") {
    if (prod(RE %in% ynames) == 0) {
      Fakename <- which(RE %in% ynames == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found in the longitudinal dataset.\n"))
    } else if (timeVar == random.var[1]) {
      Z <- Ydata[, "Ytime.aft.S"]
      Z <- cbind(1, Z)
      Z <- as.matrix(Z)
    } else if (prod(RE %in% long) == 0) {
      Fakename <- which(RE %in% long == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found in the long.formula argument. 
                  Please include this variable in the random argument.\n"))
    } else {
      Z <- Ydata[, RE]
      Z <- cbind(1, Z)
      Z <- as.matrix(Z)
    }
  } else if (model == "intercept") {
    if (!is.null(RE)) {
      stop("You are fitting a mixed effects model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
    }
    Z <- rep(1, nrow(X1))
    Z <- as.data.frame(Z)
    Z <- as.matrix(Z)
    
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }
  
  Y <- Ydata[, long[1]]
  
  YIDwID <- Ydata[, c(ID, iCen.info$weight.ID)] 
  
  
  ## obtain initial guess of baseline hazards
  subTdata <- Tdata[, c(ID, survival[2], iCen.info$weight, "T.aft.S")]
  H0 <- getBH(subTdata)
  
  ## obtain initial guess of hazard of initial event
  iCen.observed <- data.frame(unique(iCen.info$iCen.data[, iCen.info$ID]), iCen.info$iCen.observed)
  colnames(iCen.observed) <- c(iCen.info$ID, "iCen.observed")
  
  subiCendata <- iCen.data[, c(ID, iCen.info$S, iCen.info$weight)]
  subiCendata <- dplyr::left_join(subiCendata, iCen.observed, by = iCen.info$ID)
  
  subiCendata <- subiCendata[order(-subiCendata[, iCen.info$S]), ]
  #subiCendata <- subiCendata[subiCendata$iCen.observed, ]
  
  subiCendata <- as.matrix(subiCendata)
  phi <- GetrisksetS(subiCendata)
  phi <- as.data.frame(phi)
  colnames(phi) <- c(iCen.info$S, iCen.info$weight, "phi.Su")

  return(list(fixed.para = fixed.para, H0 = H0, Tdata = Tdata, Ydata = Ydata, phi = phi,
              X1 = X1, Z = Z, W = W, Y = Y, X2 = X2, sortTIDwID = sortTIDwID, YIDwID = YIDwID,
              ni = ni, survtime = survtime, status = status,
              long.formula = long.init.formula, surv.formula = survfmla, 
              variance.formula = variance.formula))
  
}