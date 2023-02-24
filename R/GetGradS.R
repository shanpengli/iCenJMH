GetGradS <- function(Psl, phi, iCen.info, nt, pStol) {
  
  Psl2 <- cbind(Psl, iCen.info$iCen.data[, iCen.info$S])
  colnames(Psl2)[4] <- colnames(phi)[1]
  Psl2 <- as.data.frame(Psl2)
  Psl2 <- dplyr::left_join(Psl2, phi[, c(1, 3)], by = colnames(phi)[1])
  GradS <- matrix(0, nrow = length(nt), ncol = nrow(phi))
  colnames(GradS) <- paste(iCen.info$S, phi[, 1], sep = "_")
  Psl2 <- Psl2[order(Psl2[, iCen.info$ID], -Psl2[, iCen.info$S]), ]
  Psl2$Grad <- 0
  Sindex <- as.data.frame(cbind(phi[, 1], 1:nrow(phi)))
  colnames(Sindex) <- c(colnames(phi)[1], "index")
  Psl2 <- dplyr::left_join(Psl2, Sindex, by = colnames(phi)[1])
  Psl2 <- as.matrix(Psl2)
  status <- getGradS(Psl2, nt, pStol, GradS)
  Psl2 <- as.data.frame(Psl2)
  ID <- unique(Psl2[, iCen.info$ID])
  GradS <- cbind(ID, GradS)
  return(GradS)
  
}

