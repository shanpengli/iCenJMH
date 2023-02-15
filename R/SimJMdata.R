##' @export
##' 

SimJMdata <- function(seed = 99, n = 100, phi = 0.04,
                      nc = 100,
                      covbw = matrix(c(1, 0.1, 0.1, 0.5), nrow = 2, ncol = 2),
                      lambda = 0.3, lambdaC = 0.05,
                      Cmin = 1,
                      Cmax = 5,
                      gamma = c(-0.05, 0.2, -0.1),
                      alpha = c(0.5, -0.5),
                      beta = c(5, 1, 2, -3, 3),
                      tau = c(-0.5, -0.1, -0.2, 0.8, 0.4),
                      increment = 2) {
  
  set.seed(seed = seed)
  ### generate interval censored event time for S_i
  ## 1. generate actual event time S_i
  # Si <- rep(0, n)
  # for (i in 1:n) {
  #   Si[i] <- round(rexp(1, rate = phi), 1) + 0.1
  #   while (Si[i] > 20) {
  #     Si[i] <- round(rexp(1, rate = phi), 1) + 0.1
  #   }
  # }
  Si <- round(runif(n, 0.1, 5), 1)
  ## 2. generate multiple 20 inspection times C_ij
  Ci <- matrix(0, nrow = n, ncol = nc)
  Ci[, 1] <- 0
  for (i in 1:nrow(Ci)) {
    for (j in 2:ncol(Ci)) {
      Ci[i, j] <- round(Ci[i, j-1] + runif(1, min = Cmin, max = Cmax), 1)
    }
  }
  
  ## 3. determine which c_ij contain S_i
  Li <- vector()
  Ri <- vector()
  for (i in 1:nrow(Ci)) {
    for (j in 1:(ncol(Ci)-1)) {
      if (Ci[i, j] < Si[i] && Si[i] <= Ci[i, j+1]) {
        Li[i] <- round(Ci[i, j], 1)
        Ri[i] <- round(Ci[i, j+1], 1)
      } else next
    }
  }
  EventTime <- data.frame(Li, Ri, Si)
  
  if (sum(complete.cases(EventTime)) != n) stop("Not enough inspection times!")
  ID <- c(1:n)
  Sdata <- data.frame(ID, EventTime)
  ### generate random effects theta_i
  bwi <- MASS::mvrnorm(n = n, c(0, 0), covbw, tol = 1e-6, empirical = FALSE)
  bwi <- as.matrix(bwi)
  # 
  ### generate the target event time T_i
  # ## 1. generate the actual event time 
  X1 <- runif(n, min = -1, max = 1)
  X2 <- rnorm(n, mean = 1, sd = 2)
  Xsurv <- cbind(Sdata$Si, X1, X2)
  Ti <- vector()
  Ci <- vector()
  
  for (i in 1:n) {
    Ti[i] <- 0
    Ci[i] <- 0
    while (Ti[i] <= (Sdata$Ri[i] - Sdata$Si[i]) || Ti[i] >= 30) {
      Ti[i] <- rexp(1, rate = lambda*exp(Xsurv[i, ]%*%gamma + bwi[i, ]%*%alpha))
    }
    while (Ci[i] <= (Sdata$Ri[i] - Sdata$Si[i]) || Ci[i] >= 30) {
      Ci[i] <- rexp(1, rate = lambdaC)
    }
  }
  
  survtimeraw <- cbind(Ti, Ci)
  status <- vector()
  survtime <- vector()
  for (i in 1:n) {
    if (min(survtimeraw[i, ]) == survtimeraw[i, 1]) {
      status[i] <- 1
      survtime[i] <- survtimeraw[i, 1]
    } else {
      status[i] <- 0
      survtime[i] <- survtimeraw[i, 2]
    }
  }
  survtimeraw <- data.frame(survtimeraw, survtime, status)
  Tdata <- data.frame(ID, survtimeraw$survtime, (survtimeraw$survtime + EventTime$Si),
                      survtimeraw$status, Xsurv)
  colnames(Tdata) <- c("ID", "Ti", "Ei", "status", "Si", "X21", "X22")
  
  table <- as.data.frame(table(Tdata$status)/n*100)
  writeLines(paste0("The censoring rate is: ", table[1, 2], "%"))
  writeLines(paste0("The target event rate is: ", table[2, 2], "%"))
  
  ### generate the longitudinal outcome y_i
  X11 <- Tdata$X21
  X12 <- Tdata$X22
  Xlong <- cbind(Tdata$Si, X11, X12)
  colnames(Xlong) <- c("Si", "X11", "X12")
  YdataRaw <- NULL
  for (i in 1:n) {
    ni <- floor((Tdata[i, 3] - Sdata[i, 3])/increment)
    suby <- matrix(0, nrow = ni+1, ncol = 4)
    suby[, 1] <- i
    ti1 <- Sdata[i, 3] - Sdata[i, 4]
    suby[1, 4] <- ti1
    sd <- sqrt(exp(tau[1] + tau[2]*Xlong[i, 1] + tau[3]*ti1 + 
                   tau[4]*Xlong[i, 2] + tau[5]*Xlong[i, 3] + bwi[i, 2]))
    suby[1, 2] <- beta[1] + beta[2]*Xlong[i, 1] + beta[3]*ti1 + 
      beta[4]*Xlong[i, 2] + beta[5]*Xlong[i, 3] + bwi[i, 1]
    epsilon <- rnorm(1, mean = 0, sd = sd)
    while (abs(epsilon) > 50) {
      epsilon <- rnorm(1, mean = 0, sd = sd)
    }
    suby[1, 2] <- suby[1, 2] + epsilon
      
    suby[1, 3] <- Sdata[i, 3]
    if (ni == 0) {
      colnames(suby) <- c("ID", "Yij", "Oij", "tij")
    } else {
      for (j in 1:ni) {
        tij <- ti1 + j*increment
        suby[j+1, 4] <- tij
        sd <- sqrt(exp(tau[1] + tau[2]*Xlong[i, 1] + tau[3]*tij + 
                         tau[4]*Xlong[i, 2] + tau[5]*Xlong[i, 3] + bwi[i, 2]))
        suby[j+1, 2] <- beta[1] + beta[2]*Xlong[i, 1] + beta[3]*tij + 
          beta[4]*Xlong[i, 2] + beta[5]*Xlong[i, 3] + bwi[i, 1]
        epsilon <- rnorm(1, mean = 0, sd = sd)
        while (abs(epsilon) > 50) {
          epsilon <- rnorm(1, mean = 0, sd = sd)
        }
        suby[j+1, 2] <- suby[j+1, 2] + epsilon
        suby[j+1, 3] <- Sdata[i, 3] + j*increment
      }
      colnames(suby) <- c("ID", "Yij", "Oij", "tij")
    }
    YdataRaw <- rbind(YdataRaw, suby)
  }
  ID <- c(1:n)
  Xlong <- cbind(ID, Xlong)
  Xlong <- as.data.frame(Xlong)
  YdataRaw <- as.data.frame(YdataRaw)
  Ydata <- dplyr::left_join(YdataRaw, Xlong, by = "ID")
  
  return(list(Ydata = Ydata, Tdata = Tdata, Sdata = Sdata))
}