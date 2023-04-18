Getriskset <- function(cdata, surv.formula) {
  
  survname <- all.vars(surv.formula)
  status <- as.vector(cdata[, survname[2]])
  
  subcdata <- cdata[, survname[1:2]]
  subcdata <- as.matrix(subcdata)
  
  riskset <- GetrisksetCSF(subcdata)
  tablerisk1 <- riskset$H01
  
  colnames(tablerisk1) <- c("survtime", "d", "hazard")
  a <- list(tablerisk1)
  names(a) <- c("tablerisk1")
  return(a)
  
}