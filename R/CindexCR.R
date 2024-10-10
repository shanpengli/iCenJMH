CindexCR <- function(time, status, predicted, Cause_int = 1) {
  if (any(is.na(time))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(status))) {
    stop("The input vector cannot have NA")
  }
  if (any(!(status %in% c(0, 1, 2)))) {
    stop("The status must be 0 or 1 or 2")
  }
  if (any(is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input of Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }
  
  Time_survival <- time
  Censoring <- ifelse(status == 0, 0, 1)
  Cause <- ifelse(status == 2, 2, 1)
  Prediction <- predicted
  Time <- max(Time_survival) + 1
  
  n <- length(Prediction)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)
  Num <- 0
  Den <- 0
  for (i in 1:n) {
    A[i, which(Time_survival[i] < Time_survival)] <- 1
    B[i, intersect(intersect(which((
      Time_survival[i] >= Time_survival
    )), which(Cause != Cause_int)), which(Censoring == 1))] <- 1
    Q[i, which(Prediction[i] > Prediction)] <- 1
  }
  for (i in 1:n) {
    if (Time_survival[i] <= Time &&
        Cause[i] == Cause_int && Censoring[i] == 1) {
      N_t[i, ] <- 1
    }
  }
  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t
  Num <- sum(Num_mat)
  Den <- sum(Den_mat)
  return(Num / Den)
}