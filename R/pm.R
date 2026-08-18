
#' Example performance metrics
#'
#' Functions that evaluate return probabilities of outcomes from the simulations.
#'
#' - `P_PNI50` calculates the probability that PNI exceeds 0.50 (threshold for an integrated-transition population, Withler et al. 2018)
#' - `P_PNI80` calculates the probability that PNI exceeds 0.80 (threshold for an integrated-wild population, Withler et al. 2018)
#' - `P_WILD50` calculates the probability that at least 50 percent of natural spawners are wild
#' - `P_SMSY80` calculates the probability that NOS exceeds 0.80 SMSY
#' - `P_SMSY85` calculates the probability that NOS/SMSY exceeds 0.85 SMSY
#' - `P_Sgen100` calculates the probability that NOS exceeds Sgen
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
P_PNI50 <- function(SMSE, Ref = 0.5, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@PNI[, , Yrs[1]:Yrs[2], drop = FALSE] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname P_PNI50
#' @export
P_PNI80 <- P_PNI50
formals(P_PNI80)$Ref <- 0.8

#' @rdname P_PNI50
#' @export
P_WILD50 <- function(SMSE, Ref = 0.50, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  apply(SMSE@p_wild[, , Yrs[1]:Yrs[2], drop = FALSE] >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname P_PNI50
#' @export
P_SMSY85 <- function(SMSE, Ref = 0.85, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  y <- seq(Yrs[1], Yrs[2])
  NOS <- apply(SMSE@NOS[, , , y, drop = FALSE], c(1, 2, 4), sum)
  SMSY <- sapply(1:SMSE@nstocks, function(s) SMSE@Misc$Ref[[s]]["Spawners_MSY", ]) %>%
    array(c(SMSE@nsim, SMSE@nstocks, length(y)))
  ratio <- NOS/SMSY
  apply(ratio >= Ref, 2, mean, na.rm = TRUE)
}

#' @rdname P_PNI50
#' @export
P_SMSY80 <- P_SMSY85
formals(P_PNI80)$Ref <- 0.8

#' @rdname P_PNI50
#' @export
P_Sgen100 <- function(SMSE, Ref = 1, Yrs = NULL) {
  if (is.null(Yrs)) Yrs <- c(1, SMSE@proyears)
  y <- seq(Yrs[1], Yrs[2])
  NOS <- apply(SMSE@NOS[, , , y, drop = FALSE], c(1, 2, 4), sum)
  Sgen <- sapply(1:SMSE@nstocks, function(s) SMSE@Misc$Ref[[s]]["Sgen", ]) %>%
    array(c(SMSE@nsim, SMSE@nstocks, length(y)))
  ratio <- NOS/Sgen
  apply(ratio >= Ref, 2, mean, na.rm = TRUE)
}

#' @name Deprecated
#' @title Deprecated functions
#' @description Deprecated performance metric functions. These functions are now replaced with functions of same name
#' but prefixed with `P_`.
#' @param ... Same arguments as [P_PNI50()] and similar functions
#' @return A vector of probabilities
#' @seealso [P_PNI50()]
#' @keywords internal
#' @export
PNI50 <- function(...) {
  .Deprecated("P_PNI50")
  P_PNI50(...)
}

#' @name Deprecated
#' @export
PNI80 <- function(...) {
  .Deprecated("P_PNI80")
  P_PNI80(...)
}

#' @name Deprecated
#' @export
WILD50 <- function(...) {
  .Deprecated("P_WILD50")
  P_WILD50(...)
}

#' @name Deprecated
#' @export
SMSY85 <- function(...) {
  .Deprecated("P_SMSY85")
  P_SMSY85(...)
}

#' @name Deprecated
#' @export
Sgen100 <- function(...) {
  .Deprecated("P_Sgen100")
  P_Sgen100(...)
}
