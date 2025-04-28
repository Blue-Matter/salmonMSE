
#' @slot u_preterminal Numeric. The exploitation rate of the immature stock in the pre-terminal fishery.
#' @slot u_terminal Numeric. The exploitation rate of the return in the terminal fishery.
#' @slot MSF_PT Logical. Whether to implement mark-selective fishing in the preterminal fishery, with no retention on unmarked fish.
#' @slot MSF_T Logical. Whether to implement mark-selective fishing in the terminal fishery, with no retention on unmarked fish.
#' @slot release_mort Vector length 2. The proportion of released fish that die after release, in the pre-terminal and terminal fishery.
#' Implemented to model mark-selective fishing. Not used if either `MSF_PT` or `MSF_T` is ` FALSE`.
#' @slot vulPT Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability schedule (between 0-1) in the preterminal fishery.
#' @slot vulT Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability schedule (between 0-1) in the terminal fishery.
