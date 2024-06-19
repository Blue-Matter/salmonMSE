
#' @slot nsim Integer. Number of simulations.
#' @slot maxage Integer. The maximum age of the population age structure.
#' @slot p_mature Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears]`. The proportion mature by age.
#' @slot SRrel Character, stock-recruit relationship for density-dependent smolt production. Either "BH" (Beverton-Holt) or "Ricker"
#' @slot capacity_smolt Vector length `nsim`. The asymptote of the Beverton-Holt stock-recruit function, or the Ricker maximum for
#'  density-dependent natural smolt production from fry life stage. **Units of smolts.**
#' @slot prod_smolt Vector length `nsim`. The productivity parameter of the Beverton-Holt stock-recruit function.
#'  Natural per-capita smolt production as the population approaches zero (density-independent component). **Units of smolts per spawner**,
#'  frequently equivalent to recruits per spawner divided by marine survival.
#' @slot a Vector length `nsim`. The alpha parameter of the Ricker stock-recruit function. **Units of smolts per egg.**
#' @slot Smax Vector length `nsim`. The Smax parameter of the Ricker stock-recruit function. **Units of eggs.**
#' @slot Mocean_NOS Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears]`.
#'  Natural mortality of immature natural origin fish.
#'  To replicate the SAR parameter of a stage-specific model, set `Mocean_NOS[a] = -log(SAR)` where `a` is the age class prior to maturation (and zero for all other ages).
#' @slot fec Vector by age (length `maxage`). Female fecundity of natural origin spawners.
#' @slot p_female Numeric. The proportion of females in the spawning population.
