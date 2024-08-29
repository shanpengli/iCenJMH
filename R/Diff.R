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

Diffrelative <- function(beta, prebeta, tau, pretau, gamma, pregamma, alpha, prealpha, 
                         Sig, preSig, H0, preH0, phi, prephi, epsilon) {
  
  betarediff <- max(abs(beta - prebeta)/(abs(prebeta) + 10*epsilon))
  taurediff <- max(abs(tau - pretau)/(abs(pretau) + 10*epsilon))
  gammarediff <- max(abs(gamma - pregamma)/(abs(pregamma) + 10*epsilon))
  alpharediff <- max(abs(alpha - prealpha)/(abs(prealpha) + 10*epsilon))
  #H0rediff <- max(abs(H0[, 3] - preH0[, 3])/(abs(preH0[, 3]) + 10*epsilon))
  Sigrediff <- max(abs(Sig - preSig))
  
  if (sum(is.na(betarediff)) || sum(is.na(taurediff)) || sum(is.na(gammarediff)) ||
      sum(is.na(alpharediff)) || sum(is.na(Sigrediff))) {
    return(2)
  } else {
    if ((betarediff > epsilon) || (taurediff > epsilon) || (gammarediff > epsilon)
        || (alpharediff > epsilon) || (Sigrediff  > epsilon)) {
      return(1)
    } else {
      return(0)
    }
  }
  
  
}