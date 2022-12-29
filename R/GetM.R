GetM <- function(GetEfun, beta, tau, gamma, alpha, Sig, Z, X1, W, Y, X2, 
                 survtime, status, TID, YID, ni, nt, YS, subiCendata, phiname) {
  
  Psl <- GetEfun$Psl
  htheta <- GetEfun$AllFUN
  ni <- as.vector(as.numeric(ni[, 2]))
  nt <- as.vector(as.numeric(nt[, 2]))
  
  if(is.list(htheta)) {
    
    ## update phi
    phi <- subiCendata
    phi[, 3] <- Psl[, 3]
    phi <- phi[order(-phi[, 2]), ]
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
    
    ID <- colnames(YS)[1]
    ## order some htheta based on survival time
    FUNE <- cbind(unique(YID[, 1]), FUNE)
    FUNBW <- cbind(unique(YID[, 1]), FUNBW)
    FUNBWE <- cbind(unique(YID[, 1]), FUNBWE)
    FUNBWSE <- cbind(unique(YID[, 1]), FUNBWSE)
    colnames(FUNE)[1] <- ID
    colnames(FUNBW) <- c(ID, 1:(ncol(Z)+1))
    colnames(FUNBWE) <- c(ID, 1:(ncol(Z)+1))
    colnames(FUNBWSE) <- c(ID, 1:((ncol(Z)+2)*(ncol(Z)+1)/2))
    FUNE <- as.data.frame(FUNE)
    FUNBW <- as.data.frame(FUNBW)
    FUNBWE <- as.data.frame(FUNBWE)
    FUNBWSE <- as.data.frame(FUNBWSE)
    FUNE <- dplyr::left_join(TID, FUNE, by = ID)
    FUNBW <- dplyr::left_join(TID, FUNBW, by = ID)
    FUNBWE <- dplyr::left_join(TID, FUNBWE, by = ID)
    FUNBWSE <- dplyr::left_join(TID, FUNBWSE, by = ID)
    
    FUNE <- as.vector(FUNE[, -c(1:2)])
    FUNBW <- as.matrix(FUNBW[, -c(1:2)])
    FUNBWE <- as.matrix(FUNBWE[, -c(1:2)])
    FUNBWSE <- as.matrix(FUNBWSE[, -c(1:2)])
    
    ## update H0 jump sizes
    w.ID <- colnames(YS)[2]
    Psl <- as.data.frame(Psl)
    subTdata <- dplyr::left_join(TID, Psl, by = c(ID, w.ID))
    subTdata <- cbind(subTdata, survtime, status)
    subTdata <- subTdata[, c(ID, "status", "psl", "survtime")]
    colnames(subTdata)[4] <- "T.aft.S"
    H0 <- getBH(subTdata)
    Psl <- as.matrix(Psl)
    PslT <- as.vector(subTdata$psl)
    
    getMpara <- getMC(beta, tau, gamma, alpha, H0, Sig, Z, X1, W, Y, X2, 
                      survtime, status, ni, nt, Psl, PslT, FUNENW, FUNEBNW, FUNEBSNW,
                      FUNE, FUNBW, FUNBWE, FUNBWSE, FUNBWS)
    
    getMpara$phi <- phi
    
    return(getMpara)
    
  } else {
    return(0);
  }
}