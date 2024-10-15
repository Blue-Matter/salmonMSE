
#' Example performance metrics
#'
#' Functions that evaluate return probabilities of outcomes from the simulations.
#'
#' - `PNI75` calculates the probability that PNI exceeds 0.75
#' - `WILD50` calculates the probability that at least 50 percent of natural spawners are wild
#' - `SMSY85` calculates the probability that NOS/SMSY exceeds 0.85
#' - `Sgen100` calculates the probability that NOS/Sgen exceeds 1
#'
#' @param SMSE SMSE object returned by [salmonMSE()]
#' @param Ref Threshold for the performance metric, used to calculate the probability that the metric exceeds this value
#' @param Yrs Numeric vector of length 2 to indicate the year range over which to summarize performance.
#' If NULL, the performance is summarized over all projection years.
#' @export
#' @return A vector of probabilities corresponding to population
PNI75 <- function(SMSE, Ref = 0.75, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@PNI[, , Yrs[1]:Yrs[2], drop = FALSE] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI75
#' @export
WILD50 <- function(SMSE, Ref = 0.50, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@p_wild[, , Yrs[1]:Yrs[2], drop = FALSE] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI75
#' @export
SMSY85 <- function(SMSE, Ref = 0.85, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  NOS <- SMSE@NOS[, , Yrs[1]:Yrs[2], drop = FALSE]
  SMSY <- SMSE@Misc$Ref["Spawners_MSY", ] # Length = 1 because ns = 1
  apply(NOS/SMSY >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI75
#' @export
Sgen100 <- function(SMSE, Ref = 1, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  NOS <- SMSE@NOS[, , Yrs[1]:Yrs[2], drop = FALSE]
  .Sgen <- SMSE@Misc$Ref["Sgen", ] # Length = 1 because ns = 1
  apply(NOS/.Sgen >= Ref, 2, mean, na.rm = TRUE)
}

