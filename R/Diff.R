Diff <- function(beta, prebeta, tau, pretau, gamma, pregamma, alpha, prealpha, 
                 Sig, preSig, H0, preH0, phi, prephi, 
                 epsilon) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  tauAbsdiff <- max(abs(tau - pretau))
  gammaAbsdiff <- max(abs(gamma - pregamma))
  alphaAbsdiff <- max(abs(alpha - prealpha))
  SigAbsdiff <- max(abs(Sig - preSig))

  if (sum(is.na(betaAbsdiff)) || sum(is.na(tauAbsdiff)) || sum(is.na(gammaAbsdiff)) ||
      sum(is.na(alphaAbsdiff)) || sum(is.na(SigAbsdiff))) {
    return(2)
  } else {
    if ((betaAbsdiff > epsilon) || (tauAbsdiff > epsilon) || (gammaAbsdiff > epsilon)
        || (alphaAbsdiff > epsilon) || (SigAbsdiff  > epsilon)) {
      return(1)
    } else {
      return(0)
    }
  }
  
  
}