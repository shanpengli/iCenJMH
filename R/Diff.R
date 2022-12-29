Diff <- function(beta, prebeta, tau, pretau, gamma, pregamma, alpha, prealpha, 
                 Sig, preSig, H0, preH0, phi, prephi, 
                 epsilon, epsilonH0) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  tauAbsdiff <- max(abs(tau - pretau))
  gammaAbsdiff <- max(abs(gamma - pregamma))
  alphaAbsdiff <- max(abs(alpha - prealpha))
  SigAbsdiff <- max(abs(Sig - preSig))
  H0Absdiff <- max(abs(H0[, 3] - preH0[, 3]))

  if (sum(is.na(betaAbsdiff)) || sum(is.na(tauAbsdiff)) || sum(is.na(gammaAbsdiff)) ||
      sum(is.na(alphaAbsdiff)) || sum(is.na(SigAbsdiff)) || sum(is.na(H0Absdiff))) {
    return(2)
  } else {
    if ((betaAbsdiff > epsilon) || (tauAbsdiff > epsilon) || (gammaAbsdiff > epsilon)
        || (alphaAbsdiff > epsilon) || (SigAbsdiff  > epsilon) || (H0Absdiff > epsilonH0)) {
      return(1)
    } else {
      return(0)
    }
  }
  
  
}