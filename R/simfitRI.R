##' @export
##'

simfitRI <- function(sim = 100, seed = 99, n = 100, phi = 0.04,
                     nc = 100,
                     covbw = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2),
                     lambda = 0.3, lambdaC = 0.05,
                     Cmin = 1,
                     Cmax = 5,
                     gamma = c(-0.05, 0.2, -0.1),
                     alpha = c(0.5, -0.5),
                     beta = c(5, 1, 2, -3, 3),
                     tau = c(-0.5, -0.1, -0.2, 0.8, 0.4),
                     increment = 2, maxiter = 1000,
                     quadpoint = 15,
                     ncores = 10) {
  
  ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfitRI,
                                      seed = seed, n = n, phi = phi,
                                      nc = nc,
                                      covbw = covbw,
                                      lambda = lambda, lambdaC = lambdaC,
                                      Cmin = Cmin,
                                      Cmax = Cmax,
                                      gamma = gamma,
                                      alpha = alpha,
                                      beta = beta,
                                      tau = tau,
                                      increment = increment, maxiter = maxiter,
                                      quadpoint = quadpoint,
                                      mc.cores = ncores)
  
  paramatrix <- as.data.frame(matrix(0, nrow = sim, ncol = 20))
  paramatrixSE <- as.data.frame(matrix(0, nrow = sim, ncol = 18))
  
  for (i in 1:sim) {
    paramatrix[i, ] <- ParaMatrixRaw[[i]]$coef
    paramatrixSE[i, ] <- ParaMatrixRaw[[i]]$coefSE
  }
  
  count <- 1
  for (i in 1:5) {
    colnames(paramatrix)[count] <- paste0("beta_", i-1)
    count <- count + 1
  }
  
  for (i in 1:5) {
    colnames(paramatrix)[count] <- paste0("tau_", i-1)
    count <- count + 1
  }
  
  for (i in 1:3) {
    colnames(paramatrix)[count] <- paste0("gamma_", i)
    count <- count + 1
  }
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("alpha_", 1)
    count <- count + 1
  }
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("Sig_", i, i)
    count <- count + 1
  }
  
  colnames(paramatrix)[count] <- paste0("Sig_12")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Time")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Iter")
  
  name <- colnames(paramatrix)[-(19:20)]
  colnames(paramatrixSE) <- paste0("se", name)
  
  result <- list(paramatrix, paramatrixSE)
  names(result) <- c("paramatrix", "paramatrixSE")
  
  return(result)
  
}