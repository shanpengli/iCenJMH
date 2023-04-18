##' @export
##' 

GetSupport <- function(iCen.data = NULL, iCen.tL = NULL, iCen.tR = NULL,
                       ID = NULL, S = "S", weight = "weight", weight.ID = "weight.ID") {
  
  ## midpoint regression
  midpoint <- (iCen.data[, iCen.tL] + iCen.data[, iCen.tR])/2
  iCen.midpoint.data <- data.frame(iCen.data[, ID], iCen.data[, iCen.tL],
                                   iCen.data[, iCen.tR], midpoint)
  colnames(iCen.midpoint.data) <- c(ID, iCen.tL, iCen.tR, S)
  
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
  SupportData <- NULL
  for (i in 1:nrow(iCen.data)) {
    if (iCen.data[i, iCen.tL] < iCen.data[i, iCen.tR]) {
      SL <- iCen.data[i, iCen.tL] < UniqueSi
      SU <- iCen.data[i, iCen.tR] >= UniqueSi
      Si <- as.logical(SL*SU)
      Si <- UniqueSi[Si]
      subdata <- data.frame(iCen.data[i, ID], Si, 1/length(Si))
      subdata$wID <- c(1:nrow(subdata))
      colnames(subdata) <- c(ID, S, weight, weight.ID)
    } else {
      subdata <- data.frame(iCen.data[i, ID], iCen.data[i, iCen.tR], 1)
      subdata$wID <- c(1:nrow(subdata))
      colnames(subdata) <- c(ID, S, weight, weight.ID)
    }
    SupportData <- rbind(SupportData, subdata)
  }
  iCen.data <- dplyr::left_join(iCen.data, SupportData, by = ID)
  
  result <- list(iCen.data = iCen.data, S = S, weight = weight, ID = ID, weight.ID = weight.ID,
                 iCen.tL = iCen.tL, iCen.tR = iCen.tR,
                 iCen.midpoint.data = iCen.midpoint.data,
                 iCen.uniform.data = iCen.uniform.data,
                 iCen.observed = iCen.observed)
  class(result) <- "iCen.info.iCenJMMLSM"
  return(result)
}


