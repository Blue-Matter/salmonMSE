
#' @slot type_PT Character. Whether to manage preterminal fishery catch from exploitation rate ("u") or catch target ("catch"). Default is "u",
#' @slot type_T Character. Whether to manage terminal fishery catch from exploitation rate ("u") or catch target ("catch"). Default is "u",
#' @slot u_preterminal Numeric. Can also be a matrix `[nsim, proyears]`. If `type_PT = "u"`, the harvest rate of the immature component of the population in the pre-terminal fishery.
#' The harvest rate is the ratio to kept AEQ catch to (kept AEQ catch + return), where AEQ are adult equivalents.
#' @slot u_terminal Numeric. Can also be a matrix `[nsim, proyears]`. If `type_T = "u"`, the harvest rate (ratio of kept catch to return) of the terminal marine fishery.
#' @slot K_PT Numeric. If `type_PT = "catch"`, the catch target of the immature component of the population in the pre-terminal fishery.
#' @slot K_T Numeric. If `type_T = "catch"`, the catch target of the return in the terminal fishery.
#' @slot MSF_PT Logical. Whether to implement mark-selective fishing in the preterminal fishery, with no retention on unmarked fish.
#' @slot MSF_T Logical. Whether to implement mark-selective fishing in the terminal fishery, with no retention on unmarked fish.
#' @slot release_mort Vector length 2. The proportion of released fish that die after release, in the pre-terminal and terminal fishery.
#' Implemented to model mark-selective fishing. Not used if either `MSF_PT` or `MSF_T` is ` FALSE`.
#' @slot vulPT Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability schedule (between 0-1) in the preterminal fishery. Values indicate
#' the proportion of fishing mortality experienced by each age class, where `F_preterminal = -log(1 - u_preterminal)`.
#' @slot vulT Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability schedule (between 0-1) in the terminal fishery. Values indicate
#' the proportion of fishing mortality experienced by each age class, where `F_terminal = -log(1 - u_terminal)`.
