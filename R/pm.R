
#' Example performance metrics
#'
#' Functions that evaluate return probabilities of outcomes from the simulations.
#'
#' - `PNI50` calculates the probability that PNI exceeds 0.50 (threshold for an integrated-transition population, Withler et al. 2018)
#' - `PNI80` calculates the probability that PNI exceeds 0.80 (threshold for an integrated-wild population, Withler et al. 2018)
#' - `WILD50` calculates the probability that at least 50 percent of natural spawners are wild
#' - `SMSY85` calculates the probability that NOS/SMSY exceeds 0.85
#' - `Sgen100` calculates the probability that NOS/Sgen exceeds 1
#'
#' @param SMSE SMSE object returned by [salmonMSE()]
#' @param Ref Threshold for the performance metric, used to calculate the probability that the metric exceeds this value
#' @param Yrs Numeric vector of length 2 to indicate the year range over which to summarize performance.
#' If NULL, the performance is summarized over all projection years.
#' @references
#' Withler et al. 2018. Genetically Based Targets for Enhanced Contributions to Canadian Pacific Chinook Salmon Populations.
#' DFO Can. Sci. Advis. Sec. Res. Doc. 2018/019. xii + 88 p.
#' @export
#' @return A vector of probabilities corresponding to population
PNI50 <- function(SMSE, Ref = 0.5, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@PNI[, , Yrs[1]:Yrs[2], drop = FALSE] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI50
#' @export
PNI80 <- PNI50
formals(PNI80)$Ref <- 0.8

#' @rdname PNI50
#' @export
WILD50 <- function(SMSE, Ref = 0.50, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@p_wild[, , Yrs[1]:Yrs[2], drop = FALSE] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI50
#' @export
SMSY85 <- function(SMSE, Ref = 0.85, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  y <- seq(Yrs[1], Yrs[2])
  NOS <- SMSE@NOS[, , y, drop = FALSE]
  SMSY <- sapply(1:SMSE@nstocks, function(s) SMSE@Misc$Ref[[s]]["Spawners_MSY", ]) %>%
    array(c(SMSE@nsim, SMSE@nstocks, length(y)))
  ratio <- NOS/SMSY
  apply(ratio >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname PNI50
#' @export
Sgen100 <- function(SMSE, Ref = 1, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  y <- seq(Yrs[1], Yrs[2])
  NOS <- SMSE@NOS[, , y, drop = FALSE]
  Sgen <- sapply(1:SMSE@nstocks, function(s) SMSE@Misc$Ref[[s]]["Sgen", ]) %>%
    array(c(SMSE@nsim, SMSE@nstocks, length(y)))
  ratio <- NOS/Sgen
  apply(ratio >= Ref, 2, mean, na.rm = TRUE)
}

