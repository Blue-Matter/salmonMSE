
#' @slot nsim Integer. Number of simulations.
#' @slot maxage Integer. The maximum age of the population age structure.
#' @slot p_mature Vector of length `maxage`. The proportion mature by age.
#' @slot SRrel Character, stock-recruit relationship for density-dependent smolt production. Either "BH" (Beverton-Holt) or "Ricker"
#' @slot capacity_smolt Vector length `nsim`. The asymptote of the Beverton-Holt stock-recruit function, or the Ricker maximum for
#'  density-dependent natural smolt production from fry life stage. Units of smolts.
#' @slot prod_smolt Vector length `nsim`. The productivity parameter of the Beverton-Holt stock-recruit function.
#'  Natural per-capita smolt production as the population approaches zero (density-independent component). Units of smolts per spawner,
#'  frequently equivalent to recruits per spawner divided by marine survival.
#' @slot SAR_NOS Vector length `nsim`. Survival to adult return of natural origin fish (density-independent, time-invariant).
#' @slot fec Numeric. Female fecundity of natural origin spawners.
#' @slot p_female Numeric. The proportion of females in the spawning population.
