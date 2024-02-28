
#' @slot Name Identifying name.
#' @slot nsim Integer. Number of simulations.
#' @slot maxage Integer. The maximum age of the population age structure.
#' @slot p_mature Vector of length `maxage`. The proportion mature by age.
#' @slot capacity_smolt Vector length `nsim`. The asymptote of the Beverton-Holt stock-recruit function for
#'  density-dependent natural smolt production from fry life stage.
#' @slot prod_smolt Vector length `nsim`. The productivity parameter of the Beverton-Holt stock-recruit function.
#'  Natural per-capita smolt production as the population approaches zero (density-independent component).
#' @slot SAR Vector length `nsim`. Survival to adult return (density-independent, time-invariant).
#' @slot fec Vector length `nsim`. Female fecundity of natural origin spawners.
#' @slot p_female Vector length `nsim`. The proportion of females in the spawning population.
