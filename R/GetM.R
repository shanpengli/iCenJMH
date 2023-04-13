GetM <- function(GetEfun, beta, tau, gamma, alpha, Sig, Z, X1, W, Y, X2, 
                 survtime, status, TID, YID, ni, nt, YS, subiCendata, phiname, 
                 pStol, iCen.observed) {
  
  Psl <- GetEfun$Psl
  htheta <- GetEfun$AllFUN
  ni <- as.vector(as.numeric(ni[, 2]))
  nt <- as.data.frame(table(YS[, 1]))
  nt <- as.vector(as.numeric(nt[, 2]))
  
  if(is.list(htheta)) {
    
    ## update phi
    phi <- subiCendata
    phi <- dplyr::left_join(phi, iCen.observed, by = colnames(phi)[1])
    phi[, 3] <- Psl[, 3]
    phi <- phi[order(-phi[, 2]), ]
    phi <- phi[phi$iCen.observed, ]
    phi <- as.matrix(phi)
    phi <- GetrisksetS(phi)
    phi <- as.data.frame(phi)
    colnames(phi) <- phiname
    
    FUNENW <- htheta$FUNENW
    FUNEBNW <- htheta$FUNEBNW
    FUNEBSNW <- htheta$FUNEBSNW
    FUNE <- htheta$FUNE
    FUNBW <- htheta$FUNBW
    FUNBWE <- htheta$FUNBWE
    FUNBWSE <- htheta$FUNBWSE
    FUNBWS <- htheta$FUNBWS
    
    IDwID <- colnames(YS)[1:2]
    ## order some htheta based on survival time
    FUNE <- cbind(unique(YID), FUNE)
    FUNBW <- cbind(unique(YID), FUNBW)
    FUNBWE <- cbind(unique(YID), FUNBWE)
    FUNBWSE <- cbind(unique(YID), FUNBWSE)
    colnames(FUNE)[1:2] <- IDwID
    colnames(FUNBW) <- c(IDwID, 1:(ncol(Z)+1))
    colnames(FUNBWE) <- c(IDwID, 1:(ncol(Z)+1))
    colnames(FUNBWSE) <- c(IDwID, 1:((ncol(Z)+2)*(ncol(Z)+1)/2))
    FUNE <- as.data.frame(FUNE)
    FUNBW <- as.data.frame(FUNBW)
    FUNBWE <- as.data.frame(FUNBWE)
    FUNBWSE <- as.data.frame(FUNBWSE)
    FUNE <- dplyr::left_join(TID, FUNE, by = IDwID)
    FUNBW <- dplyr::left_join(TID, FUNBW, by = IDwID)
    FUNBWE <- dplyr::left_join(TID, FUNBWE, by = IDwID)
    FUNBWSE <- dplyr::left_join(TID, FUNBWSE, by = IDwID)
    
    FUNE <- as.vector(FUNE[, -c(1:2)])
    FUNBW <- as.matrix(FUNBW[, -c(1:2)])
    FUNBWE <- as.matrix(FUNBWE[, -c(1:2)])
    FUNBWSE <- as.matrix(FUNBWSE[, -c(1:2)])
    
    ## update H0 jump sizes
    ID <- colnames(YS)[1]
    w.ID <- colnames(YS)[2]
    Psl <- as.data.frame(Psl)
    subTdata <- dplyr::left_join(TID, Psl, by = IDwID)
    subTdata <- cbind(subTdata, survtime, status)
    subTdata <- subTdata[, c(ID, "status", "psl", "survtime")]
    colnames(subTdata)[4] <- "T.aft.S"
    H0 <- getBH(subTdata)
    Psl <- as.matrix(Psl)
    PslT <- as.vector(subTdata$psl)
    
    getMpara <- getMC(beta, tau, gamma, alpha, H0, Sig, Z, X1, W, Y, X2, 
                      survtime, status, ni, nt, Psl, PslT, FUNENW, FUNEBNW, FUNEBSNW,
                      FUNE, FUNBW, FUNBWE, FUNBWSE, FUNBWS, pStol)
    
    getMpara$phi <- phi
    
    return(getMpara)
    
  } else {
    return(0);
  }
}