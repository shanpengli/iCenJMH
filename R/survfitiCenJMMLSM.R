##' @title Prediction in Joint Models
##' @name survfitiCenJMMLSM
##' @aliases survfitiCenJMMLSM
##' @description This function computes the conditional probability of 
##' surviving later times than the last observed time for which a longitudinal 
##' measurement was available.
##' @param object an object inheriting from class \code{iCenJMMLSM}.
##' @param seed a random seed number to proceed non-parametric bootstrap. Default is 100.
##' @param Ynewdata a data frame that contains the longitudinal and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param Tnewdata a data frame that contains the survival and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param iCennewdata a data frame that contains the interval-censored initial event information for the subjects.
##' @param u a numeric vector of times for which prediction survival probabilities are to be computed.
##' @param Last.time a numeric vector or character string. This specifies the known time at which each of 
##' the subjects in Tnewdata was known to be alive. If NULL, then this is automatically taken as the 
##' survival time of each subject. If a numeric vector, then it is assumed to be greater than or equals to the 
##' last available longitudinal time point for each subject. If a character string, then it should be 
##' a variable in Tnewdata.
##' @param obs.time a character string of specifying a longitudinal time variable in Ynewdata.
##' @param quadpoint number of quadrature points used for estimating conditional probabilities 
##' when \code{method = "GH"}. Default is NULL. If \code{method = "GH"}, then 15 is used.
##' @param pStol Tolerance parameter for the posterior probability of an initial event time Si.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{iCenJMMLSM}}
##' @export
##' 
survfitiCenJMMLSM <- function(object, seed = 100, Ynewdata = NULL, Tnewdata = NULL, iCennewdata = NULL,
                          u = NULL, Last.time = NULL, obs.time = NULL, 
                          quadpoint = NULL, pStol = 1e-2, ...) {
  if (!inherits(object, "iCenJMMLSM"))
    stop("Use only with 'iCenJMMLSM' objects.\n")
  if (is.null(Ynewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(Tnewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(u)) 
    stop("Please specify the future time for dynamic prediction.")
  if (!is.vector(u)) 
    stop("u must be vector typed.")
  if (is.null(quadpoint)) {
    quadpoint <- object$quadpoint
  }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from Ydatanew.")
  } else {
    if (!obs.time %in% colnames(Ynewdata)) {
      stop(paste0(obs.time, " is not found in Ynewdata."))
    }
  }
  
  bvar <- all.vars(object$random)
  ID <- bvar[length(bvar)]
  if (!(ID %in% colnames(Ynewdata)))
    stop(paste("The ID variable", ID, "is not found in Ynewdata."))
  if (!(ID %in% colnames(Tnewdata)))
    stop(paste("The ID variable", ID, "is not found in Tnewdata."))
  if (!(ID %in% colnames(iCennewdata)))
    stop(paste("The ID variable", ID, "is not found in Tnewdata."))
  
  Ynewdata <- Ynewdata[, colnames(object$Ydata)]
  Tnewdata <- Tnewdata[, colnames(object$Tdata)]
  
  Ydata2 <- rbind(object$Ydata, Ynewdata)
  Tdata2 <- rbind(object$Tdata, Tnewdata)
  
  getdum <- getdummy(long.formula = object$long.formula,
                     surv.formula = object$surv.formula,
                     variance.formula = object$variance.formula,
                     random = object$random, Ydata = Ydata2, Tdata = Tdata2)
  

  Ydata.mean <- getdum$Ydata.mean
  Ydata.variance <- getdum$Ydata.variance
  Tdata2 <- getdum$Tdata
  
  Yvar <- colnames(Ydata.mean)[-1]
  Cvar <- colnames(Tdata2)[-1]

  ny <- nrow(Ynewdata)
  nc <- nrow(Tnewdata)
  Ny <- nrow(Ydata2)
  Nc <- nrow(Tdata2)
  
  Sig <- object$Sig
  p1a <- ncol(Sig) - 1
  if (p1a == 2) Sigb <- Sig[1:2, 1:2]
  if (p1a == 1) Sigb <- as.matrix(Sig[1, 1])
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sigb = Sigb)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  nsig <- p1a + 1
  
  Ynewdata.mean <- Ydata.mean[c((Ny-ny+1):Ny), ]
  Ynewdata.variance <- Ydata.variance[c((Ny-ny+1):Ny), ]
  Tnewdata <- Tdata2[c((Nc-nc+1):Nc), ]
  
  if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
  yID <- unique(Ynewdata.mean[, ID])
  N.ID <- length(yID)
  cID <- Tnewdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in Ydata doesn't match with Tnewdata.")
  }
  
  if (!is.null(Last.time)) {
    if (is.character(Last.time)) {
      if (Last.time %in% colnames(Tnewdata)) {
        Last.time <- Tnewdata[, Last.time]
      } else {
        stop(paste(Last.time, "is not found in Tnewdata."))
      }
    }
    if (is.numeric(Last.time) && (length(Last.time) != nrow(Tnewdata)))
      stop("The last.time vector does not match Tnewdata.")
  } else {
    Last.time <- Tnewdata[, Cvar[1]]
  }
  
  Pred <- list()
  
  Predraw <- matrix(0, nrow = nrow(Tnewdata), ncol = length(u))
  beta <- object$beta
  tau <- object$tau
  gamma <- object$gamma
  alpha <- object$alpha
  H0 <- object$H0
  ## obtain phi before calculating p_i(sl|data)
  phi <- object$phi
  ## clean H0 with small weights
  H0 <- H0[H0[, 2] > pStol, ]
  
  y.obs <- list()
  lengthu <- length(u)
  
  for (j in 1:N.ID) {
      subNDy.mean <- Ynewdata.mean[Ynewdata.mean[, ID] == yID[j], ]
      subNDy.variance <- Ynewdata.variance[Ynewdata.variance[, ID] == yID[j], ]
      subNDc <- Tnewdata[Tnewdata[, ID] == yID[j], ]
      subobs.time <- Ynewdata[Ynewdata.mean[, ID] == yID[j], obs.time]
      y.obs[[j]] <- data.frame(subobs.time, subNDy.mean[, Yvar[1]])
      tL <- iCennewdata[iCennewdata[, ID] == yID[j], object$iCen.info$iCen.tL]
      tR <- iCennewdata[iCennewdata[, ID] == yID[j], object$iCen.info$iCen.tR]
      piSl <- phi[phi[, object$iCen.info$S] <= tR & phi[, object$iCen.info$S] > tL, ]
      if (nrow(piSl) == 0) {
        Pred[[j]] <- data.frame(u, Predraw[j, ])
        colnames(Pred[[j]]) <- c("times", "PredSurv")
        Pred[[j]][, 2] <- NaN
      } else {
        pSLR <- vector()
        phisu <- 0
        phisusum <- sum(piSl[, 3])
        for (sl in 1:nrow(piSl)) {
          pSLR[sl] <- exp(-phisu)*(1 - exp(-piSl[sl, 3]))/(1-exp(-phisusum))
          phisu <- phisu + piSl[sl, 3]
        }
        last.time.minus.sl <- Last.time[j] - piSl[, 1]
        CH0 <- vector()
        for (sl in 1:nrow(piSl)) CH0[sl] <- CH(H0, last.time.minus.sl[sl])
        CH0u <- matrix(NA, nrow = nrow(piSl), ncol = lengthu)
        for (jj in 1:lengthu) {
          u.minus.sl <- u[jj] - piSl[, 1]
          for (sl in 1:nrow(piSl)) CH0u[sl, jj] <- CH(H0, u.minus.sl[sl])
        }
        colnames(CH0u) <- u
        rownames(CH0u) <- piSl[, 1]
        
        Y <- subNDy.mean[, Yvar[1]]
        X <- subNDy.mean[, Yvar[2:length(Yvar)]]
        if (!is.null(object$int.time.Var)) {
          int.index.X <- which(colnames(X) %in% object$int.time.Var) 
        } else {
          int.index.X <- 0
        }
        X <- as.matrix(X)
        W <- subNDy.variance[, -1]
        if (!is.null(object$int.time.Var)) {
          int.index.W <- which(colnames(W) %in% object$int.time.Var) 
        } else {
          int.index.W <- 0
        }
        W <- as.matrix(W)
        if (nsig == 2) {
          Z <- matrix(1, ncol = 1, nrow = length(Y))
        } else {
          Z <- data.frame(1, subNDy.mean[, bvar1])
          Z <- as.matrix(Z)
        }
        X2 <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])
        
        ## add space for interval-censored covariates
        if (is.null(object$int.time.Var)) {
          if (nrow(X) == 1) {
            tX <- matrix(0, nrow = 1, ncol = 1+2+ncol(X)-1)
            tX[1, 1] <- 1
            tX[1, (1+2+1):(1+2+ncol(X)-1)] <- X[1, 2:ncol(X)]
            X <- tX
          } else {
            X <- cbind(X[, 1], 0, 0, X[, 2:ncol(X)])
          }
        } else {
          if (nrow(X) == 1) {
            tX <- matrix(0, nrow = 1, ncol = 1+3+ncol(X)-1)
            tX[1, 1] <- 1
            tX[1, (1+2+1):(1+2+ncol(X)-1)] <- X[1, 2:ncol(X)]
            X <- tX
          } else {
            X <- cbind(X[, 1], 0, 0, X[, 2:ncol(X)], 0)
          }
          int.index.X <- 3 + int.index.X - 1 
        }
        
        
        if (is.null(object$int.time.Var)) {
          if (nrow(W) == 1) {
            tW <- matrix(0, nrow = 1, ncol = 1+2+ncol(W)-1)
            tW[1, 1] <- 1
            tW[1, (1+2+1):(1+2+ncol(W)-1)] <- W[1, 2:ncol(W)]
            W <- tW
          } else {
            W <- cbind(W[, 1], 0, 0, W[, 2:ncol(W)])
          }
        } else {
          if (nrow(W) == 1) {
            tW <- matrix(0, nrow = 1, ncol = 1+3+ncol(W)-1)
            tW[1, 1] <- 1
            tW[1, (1+2+1):(1+2+ncol(W)-1)] <- W[1, 2:ncol(W)]
            W <- tW
          } else {
            W <- cbind(W[, 1], 0, 0, W[, 2:ncol(W)], 0)
          }
          int.index.W <- 3 + int.index.W - 1
        }
        
        X2 <- cbind(0, X2)
        
        for (jj in 1:lengthu) {
          Predraw[j, jj] <- getES(beta, tau, gamma, alpha, Sig, Z, X, W, Y,
                                  as.vector(X2), subobs.time, xsmatrix, wsmatrix, pSLR, 
                                  piSl[, 1], CH0, CH0u[, jj], int.index.X, int.index.W)
        }
        
        Pred[[j]] <- data.frame(u, Predraw[j, ])
        colnames(Pred[[j]]) <- c("times", "PredSurv")
      }

  }
  
  names(y.obs) <- names(Pred) <- yID
  Last.time <- data.frame(cID, Last.time)
  colnames(Last.time)[1] <- ID
  sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, method = "GH", quadpoint = quadpoint)
  class(sum) <- "survfitiCenJMMLSM"
  sum
  
}