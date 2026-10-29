
#' @slot maxage Integer. The maximum age of the population age structure. Default is 2.
#' @slot n_g Integer. Number of life history groups within a cohort.
#' Life history groups (LHGs) are sub-units of a cohort that have different marine survival.
#' Default is 1.
#' @slot p_LHG The proportion of the total egg production assigned to each life history group within a cohort. For example,
#' if `Bio@n_g <- 2`, then `Bio@p_LHG <- c(0.9, 0.1)`, then 90 percent of the egg production in the first population is assigned to the first life history group and
#' ten percent to the second LHG. Default is `rep(1/Bio@n_g, Bio@n_g)`
#' @slot p_mature Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, proyears]`. The proportion mature by age.
#' Default is 1 at the maximum age, 0 otherwise.
#' @slot SRrel Character, stock-recruit relationship for density-dependent smolt production from fry. Either "Ricker" (default) or "BH".
#'  Not used if habitat component is used. See `Habitat` object.
#' @slot capacity Vector length `nsim`. Only used if `SRrel = "BH"`. The asymptote, i.e., maximum juvenile production, of the Beverton-Holt stock-recruit function for density-dependent natural smolt production from egg production.
#'  **Units of smolts.** Default is 1000. Not used if habitat component is used.
#' @slot kappa Vector length `nsim`. The adult productivity ratio for the stock-recruit function. **Units of recruits per spawner.**
#'  Natural per-capita production of recruits as the population approaches zero (density-independent component). Default is 3.
#'  Not used if habitat component is used.
#' @slot Smax Vector length `nsim`. Only used if `SRrel = "Ricker"`. The spawner abundance that maximizes smolt production in the Ricker stock-recruit function. **Units of spawners.**
#'  Default is 1000. Not used if habitat component is used.
#' @slot phi Optional, vector length `nsim`. Egg production per smolt at unfished replacement. **Units of egg per smolt**.
#'  Converts productivity from units of recruits/spawner (`kappa`) to smolts/egg (`alpha`), where `alpha = kappa/phi`. In simple models,
#'  `phi` is the product of marine survival, fecundity, and proportion female. If not provided, `phi` will be calculated from `Mjuv_NOS`, `p_mature`, `s_enroute`, `p_female`, `fec`, and `p_LHG` corresponding
#'  to the first year and weighted by life history groups. Not used if habitat component is used.
#' @slot tau Optional, vector length `nsim`. Spawner per smolt at unfished replacement, only used if `SRrel = "Ricker"`. **Units of spawner per smolt**.
#'  Used to convert `Smax` to `Emax`, the corresponding egg production that maximizes smolt production, where `Emax = Smax * phi/tau` and the Ricker parameter `beta = 1/Emax`. In simple models,
#'  `tau` is the product of marine survival and proportion female. If not provided, `tau` will be calculated from `Mjuv_NOS`, `p_mature`, `s_enroute`, `p_female`, and `p_LHG` corresponding
#'  to the first year and weighted by life history groups. Not used if habitat component is used.
#' @slot SRdev Optional, matrix `[nsim, proyears]` of multiplicative deviations from the stock-recruit relationship (default is 1 for no deviations).
#' @slot Mjuv_NOS Either vector by age (length `maxage-1`) or an array with dimension `[nsim, maxage-1, proyears, n_g]`.
#'  Instantaneous natural mortality of immature natural-origin fish, the value for the first age represents natural mortality from age 1 to 2, second age is mortality from age 2 to 3, and so on.
#'  Age class survival is `exp(-Mjuv_NOS[a])`, e.g., 1 percent survival should be `Mjuv_NOS[1] <- -log(0.01)`.
#'  Default is 0 for all age classes to model spawner-recruit dynamics to bypass juvenile life stages.
#' @slot fec Vector by age (length `maxage`) or an array with dimension `[nsim, maxage, proyears]`. Female fecundity of natural spawners. Default is 1 at the maximum age, 0 otherwise (used to model spawner-recruit dynamics instead of egg production).
#' @slot p_female Numeric. The proportion of females in the spawning population. Can also be a vector `[maxage]` (for situations where older spawners are predominantly female). Default is 1 to model spawner-recruit dynamics instead of egg production.
#' @slot s_enroute Numeric. Survival of escapement to the spawning grounds (for spawning and for broodtake). Default is 1.
