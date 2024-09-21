
#' Example performance metrics
#'
#' Functions that evaluate return probabilities of outcomes from the simulations.
#'
#' - `PNI75` calculates the probability that PNI exceeds 0.75
#' - `WILD80` calculates the probability that at least 80 percent of natural spawners are wild
#'
#' @param SMSE SMSE object returned by [salmonMSE()]
#' @param Ref Threshold for the performance metric, used to calculate the probability that the metric exceeds this value
#' @param Yrs Numeric vector of length 2 to indicate the year range over which to summarize performance.
#' If NULL, the performance is summarized over all projection years.
#' @export
#' @return A vector of probabilites corresponding to population
PNI75 <- function(SMSE, Ref = 0.75, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@PNI[,,Yrs[1]:Yrs[2]] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI75
#' @export
WILD80 <- function(SMSE, Ref = 0.80, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@p_wild[,,Yrs[1]:Yrs[2]] >= Ref, 2, mean, na.rm = TRUE)
}
