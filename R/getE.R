GetE <- function(beta, tau, gamma, alpha, H0, Sig, phi, Z, X1, W, Y,
                 X2, survtime, status, TID, YID, ni, nt, YS, xsmatrix, wsmatrix,
                 S, iCen.ID, iCen.wID, pStol, c, hazard.kernel) {

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
  CUH0crude <- rep(0, n)
  HAZ0crude <- rep(0, n)
  
  getHazard(CumuH0, survtime, status, H0, CUH0crude, HAZ0crude)
  kernelH0 <- H0
  ## kernel estimate of smoothing hazards
  ### bandwidth
  b = c*H0[nrow(H0),1]/(8*nrow(H0)^0.2)
  uncensoredsurvtime <- cbind(survtime, status)
  uncensoredsurvtime <- uncensoredsurvtime[uncensoredsurvtime[, 2] == 1, ]
  uncensoredsurvtime <- unique(uncensoredsurvtime[order(uncensoredsurvtime[, 1]), ])
  for (i in 1:nrow(uncensoredsurvtime)) {
    x <- as.numeric(abs(rep(uncensoredsurvtime[i, 1], nrow(H0)) - H0[, 1])/b <= 1)
    if (hazard.kernel == "Epanechnikov") {
      kx <- 0.75*(1-(abs(rep(uncensoredsurvtime[i, 1], nrow(H0)) - H0[, 1])/b)^2)
    } else if (hazard.kernel == "uniform") {
      kx <- 0.5
    } else if (hazard.kernel == "biweight") {
      kx <- 15/16*(1-(abs(rep(uncensoredsurvtime[i, 1], nrow(H0)) - H0[, 1])/b)^2)^2
    } else {
      stop("Please specify one of the following kernels: Epanechnikov, uniform, biweight.")
    }
    H0smooth <- cbind(H0, x, kx)
    newHazard <- H0smooth[, 4]*H0smooth[, 5]*H0smooth[, 3]
    H0smooth <- cbind(H0smooth, newHazard)
    kernelH0[i, 3] <- sum(H0smooth[, 6])/b
  }
  CumukernelH0 <- cumsum(kernelH0[, 3])
  CUH0 <- rep(0, n)
  HAZ0 <- rep(0, n)
  getHazard(CumukernelH0, survtime, status, kernelH0, CUH0, HAZ0)
  
  ## recover HAZ0crude for exact observed initial event time
  HAZ <- cbind(TID, HAZ0crude, HAZ0)
  HAZ <- dplyr::left_join(HAZ, nt, by = colnames(nt)[1])
  HAZ$newHAZ0 <- ifelse(HAZ$n == 1, HAZ$HAZ0crude, HAZ$HAZ0)
  
  H0T <- cbind(TID, CUH0crude, HAZ$newHAZ0)
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
  
  FUNENW <- rep(0, nrow(Psl))
  FUNEBNW <- matrix(0, nrow = nrow(Psl), ncol = ncol(Z))
  FUNEBSNW <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+1)*ncol(Z)/2)
  FUNE <- rep(0, nrow(Psl))
  FUNBW <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+1))
  FUNBWE <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+1))
  FUNBWSE <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+2)*(ncol(Z)+1)/2)
  FUNBWS <- matrix(0, nrow = nrow(Psl), ncol = (ncol(Z)+2)*(ncol(Z)+1)/2)
  
  AllFUN <- getEC(beta, tau, gamma, alpha, H0Y, Sig, X1, Z, W, Y, X2, survtime, 
                  status, ni, nt, xsmatrix, wsmatrix, pSLR, Psl, FUNENW, FUNEBNW,
                  FUNEBSNW, FUNE, FUNBW, FUNBWE, FUNBWSE, FUNBWS, pStol)
  
  if (AllFUN == 0) {
    AllFUN <- list(FUNENW = FUNENW, FUNEBNW = FUNEBNW, FUNEBSNW = FUNEBSNW,
                   FUNE = FUNE, FUNBW = FUNBW, FUNBWE = FUNBWE, FUNBWSE = FUNBWSE,
                   FUNBWS = FUNBWS)
  }
  
  res <- list(Psl = Psl, AllFUN = AllFUN, H0T = H0T, H0Y = H0Y)
  
  return(res)
  
}