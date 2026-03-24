#' Check OpenMP status
#'
#' Returns the number of OpenMP threads available at compile time.
#' A value greater than 1 usually indicates OpenMP is enabled.
#' A value of 1 indicates serial execution.
#'
#' @return An integer.
#' @export
openmp_status <- function() {
  openmp_status_cpp()
}