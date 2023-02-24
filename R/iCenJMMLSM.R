##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes
##' @param Ydata a longitudinal data frame in long format.
##' @param Tdata a survival data frame.
##' Each subject has one data entry.
##' @param long.formula a formula object with the response variable and time-independent fixed effects covariates
##' to be included in the longitudinal sub-model.
##' @param surv.formula a formula object with the survival time, event indicator, and the time-independent fixed effects covariates
##' to be included in the survival sub-model.
##' @param variance.formula an one-sided formula object with the fixed effects covariates to model the variance of longituidnal sub-model.
##' @param random a one-sided formula object describing the random effects part of the longitudinal sub-model.
##' For example, fitting a random intercept model takes the form ~ 1|ID.
##' Alternatively. Fitting a random intercept and slope model takes the form ~ x1 + ... + xn|ID.
##' @param timeVar a string of an interval-censored time covariate.
##' @param maxiter the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.
##' @param epsilon Tolerance parameter for parametric components. Default is 0.0001.
##' @param epsilonH0 Tolerance parameter for non-parametric Breslow estimator for the event of interest. Default is 0.0001.
##' @param pStol Tolerance parameter for the posterior probability of an initial event time Si. 
##' If the posterior probability is less than a specified value, then override the probability and move forward to 
##' the next probability in the E-step. Default is 1e-6. 
##' @param quadpoint the number of Gauss-Hermite quadrature points. Default is 20.
##' @param print.para Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.
##' @param initial.para Input initial estimate of parameters. Default is FALSE.
##' @export
##'

iCenJMMLSM <- function(Ydata = NULL, Tdata = NULL,
                       long.formula = NULL,
                       surv.formula = NULL,
                       variance.formula = NULL,
                       random = NULL,
                       timeVar = NULL,
                       iCen.info = NULL,
                       maxiter = 1000, 
                       epsilon = 1e-04, epsilonH0 = 1e-04,
                       pStol = 1e-6,
                       quadpoint = 10, print.para = FALSE,
                       initial.para = TRUE) {
  
  if (!inherits(long.formula, "formula") || length(long.formula) != 3) {
    stop("\nLinear mixed effects model must be a formula of the form \"resp ~ pred\".\n")
  }
  
  if (!inherits(surv.formula, "formula") || length(surv.formula) != 3) {
    stop("\nCox proportional hazards model must be a formula of the form \"Surv(.,.) ~ pred\".\n")
  }
  
  if (!inherits(iCen.info, "iCen.info.iCenJMMLSM")) {
    stop("Use only with 'iCen.info.iCenJMMLSM' objects.\n")
  }
  
  cnames <- colnames(Tdata)
  ynames <- colnames(Ydata)
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  variance <- all.vars(variance.formula)
  random.form <- all.vars(random)
  ID <- random.form[length(random.form)]
  if (length(random.form) == 1) {
    p1a <- 1
    RE <- NULL
    model <- "intercept"
  } else {
    RE <- random.form[-length(random.form)]
    model <- "interslope"
    p1a <- length(RE)
    if (p1a > 2)
      stop("The total dimension of random effects of the mean trajectory cannot exceed 2.")
  }
  
  ##variable check
  if (prod(long %in% ynames) == 0) {
    Fakename <- which(long %in% ynames == FALSE)
    stop(paste0("The variable ", long[Fakename], " not found in Ydata.\n"))
  }
  if (prod(survival %in% cnames) == 0) {
    Fakename <- which(survival %in% cnames == FALSE)
    stop(paste0("The variable ", survival[Fakename], " not found in Tdata.\n"))
  }
  if (!(ID %in% ynames)) {
    stop(paste0("ID column ", ID, " not found in Ydata.\n"))
  }
  if (!(ID %in% cnames)) {
    stop(paste0("ID column ", ID, " not found in Tdata.\n"))
  }
  if (!(timeVar %in% ynames)) {
    stop(paste0("The an interval-censored variable ", timeVar, " is not found in Ydata.\n"))
  }
  
  if (is.null(RE) & model == "interslope") {
    stop("Random effects covariates must be specified.\n")
  }
  if (prod(variance %in% ynames) == 0) {
    Fakename <- which(variance %in% ynames == FALSE)
    stop(paste0("The WS variables ", long[Fakename], " not found.\n"))
  }

  getinit <- Getinit(Tdata = Tdata, Ydata = Ydata, long.formula = long.formula,
                     surv.formula = surv.formula, variance.formula = variance.formula,
                     model = model, RE = RE, random = random, timeVar = timeVar,
                     initial.para = initial.para, iCen.info = iCen.info)
  
  
  ## initialize parameters
  beta <- getinit$fixed.para$beta
  tau <- getinit$fixed.para$tau
  gamma <- getinit$fixed.para$gamma
  alpha <- getinit$fixed.para$alpha
  Sig <- getinit$fixed.para$Sig
  p1a <- ncol(Sig) - 1
  if (p1a == 2) Sigb <- Sig[1:2, 1:2]
  if (p1a == 1) Sigb <- as.matrix(Sig[1, 1])
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sigb = Sigb)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  ## Extract data
  
  iter=0
  subiCendata <- unique(getinit$Ydata[, c(iCen.info$ID, iCen.info$S, iCen.info$weight)])
  phi <- getinit$phi
  phiname <- colnames(phi)
  
  Y <- getinit$Y
  X1 <- getinit$X1
  Z <- getinit$Z
  W <- getinit$W
  
  survtime <- getinit$survtime
  status <- getinit$status
  X2 <- getinit$X2
  H0 <- getinit$H0
  
  TID <- getinit$sortTIDwID
  YID <- getinit$YIDwID
  
  YS <- iCen.info$iCen.data[, c(iCen.info$ID, iCen.info$weight.ID, iCen.info$S)]
  S <- iCen.info$S
  iCen.ID <- iCen.info$ID
  iCen.wID <- iCen.info$weight.ID
  
  ni <- getinit$ni
  nt <- as.data.frame(table(YS[, iCen.info$ID]))

  repeat
  {
    iter <- iter + 1
    prebeta <- beta
    pretau <- tau
    pregamma <- gamma
    prealpha <- alpha
    preSig <- Sig
    preH0 <- H0
    prephi <- phi
    
    if (print.para) {
      writeLines("iter is:")
      print(iter)
      writeLines("beta is:")
      print(beta)
      writeLines("tau is:")
      print(tau)
      writeLines("gamma is:")
      print(gamma)
      writeLines("alpha is:")
      print(alpha)
      writeLines("Sig is:")
      print(Sig)
      }

    GetEfun <- GetE(beta, tau, gamma, alpha, H0, Sig, phi, Z, X1, W, Y,
                    X2, survtime, status, TID, YID, ni, nt, YS, xsmatrix, wsmatrix,
                    S, iCen.ID, iCen.wID, pStol)

    GetMpara <- GetM(GetEfun, beta, tau, gamma, alpha, Sig, Z, X1, W, Y, X2, 
                     survtime, status, TID, YID, ni, nt, YS, subiCendata, phiname, pStol)
    
      beta <- GetMpara$beta
      tau <- GetMpara$tau
      Sig <- GetMpara$Sig
      gamma <- GetMpara$gamma
      alpha <- GetMpara$alpha
      H0 <- GetMpara$H0
      phi <- GetMpara$phi
      
        if((Diff(beta, prebeta, tau, pretau, gamma, pregamma, alpha, prealpha, 
                 Sig, preSig, H0, preH0, phi, prephi, epsilon, epsilonH0) != 1)
            || (iter == maxiter) || (!is.list(GetEfun$AllFUN)) || (!is.list(GetMpara))) {
          break
        }
  }
  
  if (Diff(beta, prebeta, tau, pretau, gamma, pregamma, alpha, prealpha, 
           Sig, preSig, H0, preH0, phi, prephi, epsilon, epsilonH0) == 2) {
    writeLines("program stops because of numerical issues")
    convergence = 0
    beta <- NULL
    tau <- NULL
    gamma <- NULL
    alpha <- NULL
    H0 <- NULL
    Sig <- NULL
    phi <- NULL
    iter <- NULL
    result <- list(beta = beta, tau = tau, gamma = gamma, alpha = alpha, 
                   H0 = H0, Sig = Sig, phi = phi, iter = iter, convergence)
  } else if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta = beta, tau = tau, gamma = gamma, alpha = alpha, 
                     H0 = H0, Sig = Sig, phi = phi, iter = iter, convergence)
  } else if (!is.list(GetEfun$AllFUN)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      tau <- NULL
      gamma <- NULL
      alpha <- NULL
      H0 <- NULL
      Sig <- NULL
      phi <- NULL
      iter <- NULL
      result <- list(beta = beta, tau = tau, gamma = gamma, alpha = alpha, 
                     H0 = H0, Sig = Sig, phi = phi, iter = iter, convergence)
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      tau <- NULL
      gamma <- NULL
      alpha <- NULL
      H0 <- NULL
      Sig <- NULL
      phi <- NULL
      iter <- NULL
      result <- list(beta = beta, tau = tau, gamma = gamma, alpha = alpha, 
                     H0 = H0, Sig = Sig, phi = phi, iter = iter, convergence)
      return(result)
    } else {

        convergence = 1
        
        GetEfun <- GetE(beta, tau, gamma, alpha, H0, Sig, phi, Z, X1, W, Y,
                        X2, survtime, status, TID, YID, ni, nt, YS, xsmatrix, wsmatrix,
                        S, iCen.ID, iCen.wID, pStol)
        
        GetfunE <- GetEfunSE(GetEfun, Z, TID, YID, ni, nt, YS, subiCendata)
        
        Psl <- GetfunE$Psl
        PslT <- GetfunE$PslT
        FUNENW <- GetfunE$FUNENW
        FUNEBNW <- GetfunE$FUNEBNW
        FUNEBSNW <- GetfunE$FUNEBSNW
        FUNE <- GetfunE$FUNE
        FUNBW <- GetfunE$FUNBW
        FUNBWE <- GetfunE$FUNBWE
        FUNBWSE <- GetfunE$FUNBWSE
        FUNBWS <- GetfunE$FUNBWS
        
        ni <- as.vector(as.numeric(ni[, 2]))
        nt <- as.vector(as.numeric(nt[, 2]))
        
        GradY <- getGradY(beta, tau, Sig, Z, X1, W, Y, ni, nt, Psl, FUNENW, FUNEBNW, 
                          FUNEBSNW, FUNBWS, pStol)
        
        GradT <- getGradT(gamma, alpha, H0, X2, survtime, status, ni, nt, PslT,
                          FUNE, FUNBW, FUNBWE)
        
        GradT <- cbind(TID, GradT)
        GradY <- as.data.frame(cbind(unique(YID[, 1]), GradY))
        colnames(GradY)[1] <- iCen.info$ID
        
        GradT <- aggregate(GradT[, -c(1:2)], by=list(GradT[, 1]), FUN = sum)
        colnames(GradT)[1] <- iCen.info$ID
        
        GradS <- GetGradS(Psl, phi, iCen.info, iCen.data, nt, pStol)
        colnames(GradS)[1] <- iCen.info$ID
        GradS <- as.data.frame(GradS)
        zeroindex <- which(abs(colSums(GradS[, -1])) <= pStol) + 1
        GradS <- GradS[, -zeroindex]
        
        Grad <- dplyr::left_join(GradY, GradT, by = iCen.info$ID)
        Grad <- dplyr::left_join(Grad, GradS, by = iCen.info$ID)
        Grad <- Grad[, -1]
        Grad <- as.matrix(Grad)
        Cov <- GetCov(Grad) 
        
        nbeta <- length(beta)
        ntau <- length(tau)
        nsig <- ncol(Sig)
        ngamma <- length(gamma)
        nalpha <- length(alpha)
        
        SE <- GetSE(nbeta, ntau, nsig, ngamma, nalpha, Cov)
        
        sebeta <- SE$sebeta
        setau <- SE$setau
        seSig <- SE$seSig
        segamma <- SE$segamma
        sealpha <- SE$sealpha

        long <- all.vars(long.formula)
        survival <- all.vars(surv.formula)

        names(beta) <- c("intercept", long[-1])
        variance.form <- all.vars(variance.formula)
        names(tau) <- paste0(c("intercept", variance.form), "_var")

        names(gamma) <- paste0(survival[-(1:2)], "_1")

        PropComp <- as.data.frame(table(Tdata[, survival[2]]))

        LongOut <- long[1]
        LongX <- paste0(long[-1], collapse = "+")
        FunCall_long <- as.formula(paste(LongOut, LongX, sep = "~"))

        LongVarOut <- "log(sigma^2)"
        LongVarX <- paste0(variance.form, collapse = "+")
        FunCall_longVar <- as.formula(paste(LongVarOut, LongVarX, sep = "~"))

        SurvOut <- paste0("Surv(", survival[1], ",", survival[2], ")")
        SurvX <- paste0(survival[-(1:2)], collapse = "+")
        FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))

        ## return the joint modelling result
        mycall <- match.call()
        
        result <- list(beta = beta, tau = tau, gamma = gamma, alpha = alpha,
                       H0 = H0, Sig = Sig, vcov = Cov, phi = phi, sebeta = sebeta, setau = setau,
                       seSig = seSig, segamma = segamma, sealpha = sealpha, 
                       iter = iter, convergence = convergence,
                       quadpoint = quadpoint, Ydata, Tdata, PropComp, FunCall_long, 
                       FunCall_longVar, FunCall_survival, random, mycall, iCen.info = iCen.info)

        class(result) <- "iCenJMMLSM"

        return(result)
      }

}
