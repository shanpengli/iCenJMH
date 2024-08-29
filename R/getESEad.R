GetESEad <- function(beta, tau, gamma, alpha, H0, Sig, phi, Z, X1, W, Y,
                   X2, survtime, status, TID, YID, ni, nt, YS, xsmatrix, wsmatrix,
                   S, iCen.ID, iCen.wID, pStol) {
  
  ## obtain p(S | L, R)
  YS <- dplyr::left_join(YS, phi[, c(1, 3)], by = S)
  idsum <- aggregate(YS[, 4] ~ YS[, 1], YS, sum)
  idsum <- as.vector(idsum[, 2])
  YS <- as.matrix(YS)
  mdata <- table(YS[, 1])
  mdata <- as.data.frame(mdata)
  mdata <- as.vector(mdata[, 2])
  n <- length(mdata)
  mdataS <- rep(0, n) 
  mdataS[1] <- 1
  mdataCum <- cumsum(mdata)
  mdata2 <- mdata - 1
  mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
  
  pSLR <- updatePSLR(YS, mdata, mdataS, idsum)
  colnames(pSLR) <- c(colnames(YS)[1:3], "pSLR")
  
  ## obtain H0(T)
  CumuH0 <- cumsum(H0[, 3])
  n <- nrow(TID)
  CUH0 <- rep(0, n)
  HAZ0 <- rep(0, n)
  
  getHazard(CumuH0, survtime, status, H0, CUH0, HAZ0)
  
  H0T <- cbind(TID, CUH0, HAZ0)
  YSID <- as.data.frame(YS[, c(iCen.ID, iCen.wID)])
  H0Y <- dplyr::left_join(YSID, H0T, by = c(iCen.ID, iCen.wID))
  H0Y <- as.matrix(H0Y)
  ## reorder Tdata based on Ydata
  Re.Tdata <- as.data.frame(cbind(TID, survtime, status, X2))
  Re.Tdata <- dplyr::left_join(as.data.frame(YS[, c(iCen.ID, iCen.wID)]), Re.Tdata, 
                               by = c(iCen.ID, iCen.wID))
  
  survtime <- Re.Tdata[, 3]
  status <- Re.Tdata[, 4]
  X2 <- as.matrix(Re.Tdata[, 5:ncol(Re.Tdata)])
  
  ## calculate Psl
  Psl <- YSID
  Psl$psl <- NA
  Psl <- as.matrix(Psl)
  
  ni <- as.vector(as.numeric(ni[, 2]))
  nt <- as.data.frame(table(YS[, 1]))
  nt <- as.vector(as.numeric(nt[, 2]))
  
  ### calculate rescaled quadrature points for aGH
  nsig <- ncol(Z) + 1
  n <- length(idsum)
  posterior.mode <- matrix(0, nrow = sum(nt), ncol = nsig)
  posterior.var <- matrix(0, nrow = nsig*sum(nt), ncol = nsig)
  
  counti <- 1
  countt <- 1
  for (i in 1:n) {
    
    nti <- nt[i]
    nii <- ni[i]
    for (t in 1:nti) {
      
      subY <- Y[counti:(counti+nii-1)]
      subW <- as.matrix(W[counti:(counti+nii-1), ], ncol = ncol(W))
      subX1 <- as.matrix(X1[counti:(counti+nii-1), ], ncol = ncol(X1))
      subZ <- as.matrix(Z[counti:(counti+nii-1), ], ncol = ncol(Z))
      subX2 <- t(as.matrix(X2[countt, ]))
      CH0 <- H0Y[countt, 3]
      HAZ0 <- H0Y[countt, 4]
      substatus <- status[countt]
      
      data <- list(subY, subX1, subZ, subW, subX2, CH0, HAZ0, beta, tau, gamma, 
                   alpha, Sig, substatus)
      names(data) <- c("Y", "X", "Z", "W", "X2", "CH0", "HAZ0", "beta", "tau",
                       "gamma", "alpha", "Sig", "status")
      opt <- optim(rep(0, nsig), logLik.learn, data = data, method = "BFGS", hessian = TRUE)
      posterior.mode[countt, ] <- opt$par
      posterior.var[(nsig*(countt-1) + 1):(countt * nsig), 1:nsig] <- solve(opt$hessian)
      
      counti <- counti + nii
      countt <- countt + 1
      
      
    }
    
    
  }
  
  FUNENW <- rep(0, nrow(Psl))
  FUNEBNW <- matrix(0, nrow = nrow(Psl), ncol = ncol(Z))
  FUNEBSNW <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+1)*ncol(Z)/2)
  FUNE <- rep(0, nrow(Psl))
  FUNBW <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+1))
  FUNBWE <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+1))
  FUNBWSE <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+2)*(ncol(Z)+1)/2)
  FUNBWS <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+2)*(ncol(Z)+1)/2)
  
  AllFUN <- getECad(beta, tau, gamma, alpha, H0Y, Sig, X1, Z, W, Y, X2, survtime,
                    status, ni, nt, xsmatrix, wsmatrix, pSLR, posterior.mode, posterior.var,
                    Psl, FUNENW, FUNEBNW,
                    FUNEBSNW, FUNE, FUNBW, FUNBWE, FUNBWSE, FUNBWS, pStol)
  
  if (AllFUN == 0) {
    AllFUN <- list(FUNENW = FUNENW, FUNEBNW = FUNEBNW, FUNEBSNW = FUNEBSNW,
                   FUNE = FUNE, FUNBW = FUNBW, FUNBWE = FUNBWE, FUNBWSE = FUNBWSE,
                   FUNBWS = FUNBWS)
  }
  
  res <- list(Psl = Psl, AllFUN = AllFUN, H0T = H0T, H0Y = H0Y)
  
  return(res)
  
}