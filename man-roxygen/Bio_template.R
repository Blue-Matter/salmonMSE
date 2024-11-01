
#' @slot maxage Integer. The maximum age of the population age structure.
#' @slot p_mature Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears]`. The proportion mature by age.
#' @slot SRrel Character, stock-recruit relationship for density-dependent smolt production. Either "BH" (Beverton-Holt) or "Ricker"
#' @slot capacity_smolt Vector length `nsim`. The asymptote of the Beverton-Holt stock-recruit function, or the Ricker maximum for
#'  density-dependent natural smolt production from egg production. **Units of smolts.**
#' @slot kappa Vector length `nsim`. The compensation ratio for the stock-recruit function, a.k.a. adult productivity.
#'  Natural per-capita production of recruits as the population approaches zero (density-independent component). **Units of recruits per egg**.
#'  In stage-based models, equivalent to the product of smolt productivity (smolts per spawner) and marine survival.
#' @slot Smax Vector length `nsim`. The spawning output that maximizes smolt production in the Ricker stock-recruit function. **Units of egg production.**
#' @slot phi Optional parameter, vector length `nsim`. Unfished egg production rate per smolt.
#'  The `alpha` parameter of the stock-recruit function will be the ratio of `kappa` and `phi`. In stage-based models,
#'  the product of marine survival, fecundity, and proportion female. If not provided, `phi` will be calculated from `Mjuv_NOS` corresponding
#'  to the first year.
#' @slot Mjuv_NOS Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears]`.
#'  Natural mortality of immature natural origin fish.
#'  To replicate the SAR parameter of a stage-specific model, set `Mjuv_NOS[a] = -log(SAR)` where `a` is the age class prior to maturation (and zero for all other ages).
#' @slot fec Vector by age (length `maxage`). Female fecundity of natural origin spawners.
#' @slot p_female Numeric. The proportion of females in the spawning population. Default is 0.5.
#' @slot s_enroute Numeric. Survival of escapement to the spawning grounds (for spawning and for broodtake). Default is 1.
