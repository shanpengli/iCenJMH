##' @export
##' 

GetSupport <- function(iCen.data = NULL, iCen.tL = NULL, iCen.tR = NULL,
                       ID = NULL, S = "S", weight = "weight", weight.ID = "weight.ID") {
  
  ## midpoint regression
  midpoint <- (iCen.data[, iCen.tL] + iCen.data[, iCen.tR])/2
  iCen.midpoint.data <- data.frame(iCen.data[, ID], iCen.data[, iCen.tL],
                                   iCen.data[, iCen.tR], midpoint)
  colnames(iCen.midpoint.data) <- c(ID, iCen.tL, iCen.tR, S)
  
  ## right-point regression
  rightpoint <- iCen.data[, iCen.tR]
  iCen.rightpoint.data <- data.frame(iCen.data[, ID], iCen.data[, iCen.tL],
                                   iCen.data[, iCen.tR], rightpoint)
  colnames(iCen.rightpoint.data) <- c(ID, iCen.tL, iCen.tR, S)
  
  ## uniform regression
  iCen.uniform.data <- iCen.data
  iCen.uniform.data$iCen.tL <- iCen.data[, iCen.tL]
  iCen.uniform.data$iCen.tR <- iCen.data[, iCen.tR]
  f <- Vectorize(function(iCen.tL, iCen.tR) runif(1, min = iCen.tL, max = iCen.tR), 
                 vectorize.args = c("iCen.tL", "iCen.tR"))
  
  uniform <- with(iCen.uniform.data, f(iCen.tL, iCen.tR))
  iCen.uniform.data <- data.frame(iCen.data[, ID], iCen.data[, iCen.tL],
                                  iCen.data[, iCen.tR], uniform)
  colnames(iCen.uniform.data) <- c(ID, iCen.tL, iCen.tR, S)
  
  iCen.observed <- ifelse(iCen.data[, iCen.tL] < iCen.data[, iCen.tR], TRUE, FALSE)
  
  # UniqueSi <- sort(unique(c(iCen.data[iCen.observed, iCen.tL], 
  #                           iCen.data[iCen.observed, iCen.tR])))
  
  UniqueSi <- sort(unique(c(iCen.data[, iCen.tL], 
                            iCen.data[, iCen.tR])))
  
  n <- nrow(iCen.data)
  out_list <- vector("list", n)
  for (i in seq_len(n)) {
    left  <- iCen.data[[iCen.tL]][i]
    right <- iCen.data[[iCen.tR]][i]
    idval <- iCen.data[[ID]][i]
    
    if (left < right) {
      SL <- left < UniqueSi
      SU <- right >= UniqueSi
      Si <- UniqueSi[SL & SU]
      
      subdata <- data.frame(
        id_tmp = rep(idval, length(Si)),
        s_tmp = Si,
        w_tmp = rep(1 / length(Si), length(Si)),
        wid_tmp = seq_along(Si)
      )
    } else {
      subdata <- data.frame(
        id_tmp = idval,
        s_tmp = right,
        w_tmp = 1,
        wid_tmp = 1
      )
    }
    
    names(subdata) <- c(ID, S, weight, weight.ID)
    out_list[[i]] <- subdata
  }
  SupportData <- do.call(rbind, out_list)
  rownames(SupportData) <- NULL
  
  iCen.data <- dplyr::left_join(iCen.data, SupportData, by = ID)
  
  result <- list(iCen.data = iCen.data, S = S, weight = weight, ID = ID, weight.ID = weight.ID,
                 iCen.tL = iCen.tL, iCen.tR = iCen.tR,
                 iCen.midpoint.data = iCen.midpoint.data,
                 iCen.uniform.data = iCen.uniform.data,
                 iCen.rightpoint.data = iCen.rightpoint.data,
                 iCen.observed = iCen.observed)
  class(result) <- "iCen.info.iCenJMMLSM"
  return(result)
}


