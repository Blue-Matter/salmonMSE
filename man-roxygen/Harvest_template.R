
#' @slot u_preterminal Numeric. The exploitation rate of the immature stock in the pre-terminal fishery.
#' @slot u_terminal Numeric. The exploitation rate of the return in the terminal fishery.
#' @slot MSF Logical. Whether to implement mark-selective fishing, with no retention on unmarked fish.
#' @slot m Numeric. The mark rate of hatchery origin fish, which affects fishery selectivity.
#' @slot release_mort Vector length 2. The proportion of released fish that die after release, in the pre-terminal and terminal fishery.
#' Implemented to model mark-selective fishing. Not used if `MSF = FALSE`.
#' @slot vulPT Vector length `maxage`. Vulnerability schedule (between 0-1) in the preterminal fishery.
#' @slot vulT Vector length `maxage`. Vulnerability schedule (between 0-1) in the terminal fishery.
