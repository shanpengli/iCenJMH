GetEfunSE <- function(GetEfun, Z, TID, YID, ni, nt, YS, subiCendata) {
  
  Psl <- GetEfun$Psl
  htheta <- GetEfun$AllFUN
  ni <- as.vector(as.numeric(ni[, 2]))
  nt <- as.vector(as.numeric(nt[, 2]))
  
  if(is.list(htheta)) {
    
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
    Psl <- as.matrix(Psl)
    PslT <- as.vector(subTdata$psl)
    
    
    return(list(Psl = Psl, PslT = PslT, FUNENW = FUNENW, FUNEBNW = FUNEBNW, 
                FUNEBSNW = FUNEBSNW, FUNE = FUNE, FUNBW = FUNBW,
                FUNBWE = FUNBWE, FUNBWSE = FUNBWSE, FUNBWS = FUNBWS))
    
  } else {
    return(0);
  }
}