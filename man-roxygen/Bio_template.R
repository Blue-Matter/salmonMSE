
#' @slot maxage Integer. The maximum age of the population age structure.
#' @slot n_g Integer. Number of life history groups within a cohort.
#' Life history groups (LHGs) are sub-units of a cohort that have different biological parameters, e.g., survival,
#' but the egg production and smolt production in the next generation is calculated from the sum across life history groups.
#' Default is 1.
#' @slot p_LHG The proportion of the total egg production assigned to each life history group within a cohort. For example,
#' if `Bio@n_g <- 2`, then `Bio@p_LHG <- c(0.9, 0.1)`, then 90 percent of the egg production in the first population is assigned to the first life history group and
#' ten percent to the second LHG. Default is `rep(1/Bio@n_g, Bio@n_g)`
#' @slot p_mature Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears]`. The proportion mature by age.
#' @slot SRrel Character, stock-recruit relationship for density-dependent smolt production from fry. Either "BH" (Beverton-Holt) or "Ricker".
#'  Not used if habitat component is used. See `Habitat` object.
#' @slot capacity Vector length `nsim`. The asymptote of the Beverton-Holt stock-recruit function, or the Ricker maximum for
#'  density-dependent natural smolt production from egg production. **Units of smolts.**
#'  Not used if habitat component is used.
#' @slot kappa Vector length `nsim`. The adult productivity ratio for the stock-recruit function. **Units of recruits per spawner.**
#'  Natural per-capita production of recruits as the population approaches zero (density-independent component).
#'  In stage-based models, equivalent to the product of smolt productivity (smolts per spawner) and marine survival.
#'  Not used if habitat component is used.
#' @slot Smax Vector length `nsim`. The egg production that maximizes smolt production in the Ricker stock-recruit function. **Units of eggs.**
#' Equivalent to units of spawners if `fec = 1` for all spawners. Not used if habitat component is used.
#' @slot phi Optional parameter, vector length `nsim`. Unfished egg production per smolt. **Units of egg per smolt**.
#'  The `alpha` parameter of the stock-recruit function will be the ratio of `kappa` and `phi`. In stage-based models,
#'  this is the product of marine survival, fecundity, and proportion female. If not provided, `phi` will be calculated from `Mjuv_NOS`, `p_mature`, `s_enroute`, `p_female`, `fec`, and `p_LHG` corresponding
#'  to the first year and weighted by life history groups. Not used if habitat component is used.
#' @slot Mjuv_NOS Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears, n_g]`.
#'  Natural mortality of immature natural origin fish.
#'  To replicate the SAR parameter of a stage-specific model, set `Mjuv_NOS[a] = -log(SAR)` where `a` is the age class prior to maturation (and zero for all other ages).
#' @slot fec Vector by age (length `maxage`). Female fecundity of natural origin spawners.
#' @slot p_female Numeric. The proportion of females in the spawning population. Default is 0.5.
#' @slot s_enroute Numeric. Survival of escapement to the spawning grounds (for spawning and for broodtake). Default is 1.
