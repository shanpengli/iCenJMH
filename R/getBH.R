getBH <- function(subTdata) {
  
  subTdata <- subTdata[order(-subTdata$T.aft.S), ]
  subTdata <- as.matrix(subTdata)
  riskset <- GetrisksetC(subTdata)
  H01 <- riskset$H01
  colnames(H01) <- c("T.aft.S", "n.event", "Hazard")
  return(H01)
  
}