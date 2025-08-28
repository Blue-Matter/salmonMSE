
#' @details
#' A description of the fitness parameters and hatchery dynamics is available in the
#' [online documentation](https://docs.salmonmse.com/articles/equations.html#fitness-effects-on-survival).
#' @slot n_r Integer. Number of release strategies, sub-groups of fish with different survival schedules. Default is one.
#' @slot n_yearling Vector length `n_r`. The number of yearlings released by release strategy. No hatchery is modeled if `sum(n_yearling, n_subyearling) = 0`.
#' are zero. Yearlings survival is density-independent after release. Default is zero.
#' @slot n_subyearling Vector length `n_r`. The number of subyearlings released. No hatchery is modeled if `sum(n_yearling, n_subyearling) = 0` are zero.
#' Subyearlings experience density-dependent survival in competition with natural origin young. Default is zero.
#' @slot s_prespawn Numeric. The survival of broodtake prior to egg production. `1 - s_prespawn` is the proportion of fish not used for hatchery purposes, e.g., mortality or other resesarch purposes.
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`. Default is 1.
#' @slot s_egg_smolt Numeric. The survival of eggs to the smolt life stage (for yearling release).
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`. Default is 1.
#' @slot s_egg_subyearling Numeric. The survival of eggs to subyearling life stage (for subyearling release).
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`. Default is 1.
#' @slot brood_import Matrix by age (length `maxage`) for the number of annual imported hatchery origin broodstock. Egg production is weighted by `fec_brood`.
#'  To meet hatchery production target releases, imported brood and local marked brood are used indiscriminately. Default imported brood is zero.
#' @slot Mjuv_HOS Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears, n_r]`.
#'  Natural mortality of immature hatchery origin fish.
#'  To replicate the SAR parameter of a stage-specific model, set `Mjuv_HOS[a] = -log(SAR)` for the age class prior to maturation (and zero for all other ages).
#' @slot p_mature_HOS Vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears, n_r]` for the maturity of hatchery spawners.
#' Default is set equal to `Bio@p_mature` for all release strategies.
#' @slot stray_external Matrix by age (length `maxage`) and release strategy `n_r` that denotes the annual number of hatchery origin strays from other populations/systems
#' not included in the operating model. Default is zero. External strays are added at the escapement life stage, assumed unmarked. For multi-population models with straying within the system, see also `SOM@stray` matrix.
#' @slot gamma Numeric. The relative reproductive success of hatchery origin spawners (relative to natural origin spawners). Default is 1.
#' @slot m Numeric. The mark rate of hatchery origin fish, which affects selective broodtake and fishery retention if mark-selective fishing is utilized.
#' Set m = 1 for AHA compatibility with `ptarget_NOB`. Default is zero.
#' @slot pmax_esc Numeric. The maximum proportion of total escapement (after en route mortality) that could be used as broodtake. Set to 1 for AHA compatibility. Default is 0.75.
#' @slot pmax_NOB Numeric. The maximum proportion of the natural origin escapement (after en route mortality and `pmax_esc`) to be used as broodtake. If broodstock
#'  is limited by `pmax_esc < 1`, then this parameter should be 1. Default is 1.
#' @slot ptarget_NOB Numeric. The target proportion of the natural origin broodtake relative to the overall broodtake, assuming the mark rate is 1 and natural origin fish
#' can be identified in the hatchery. The realized proportion may be lower if there are insufficient natural origin escapement. If the mark rate < 1, then
#' this target proportion identifies the proportion of unmarked fished in the broodtake. If mark rate = 0, then pNOB is equal to the proportion in the escapement. Default is 0.9.
#' @slot phatchery Numeric. Optional parameter (default is `NA`). If set to a numeric between 0-1, this value is the proportion of the hatchery origin escapement that return to the hatchery, for example, by removal from spawning grounds
#'  or swim-in facilities. These fish are available for broodtake. None of these fish will spawn in the natural environment.
#'  With the default option, `NA` allows all hatchery origin escapement that is not used brood to go to the spawning grounds.
#' @slot premove_HOS Numeric. The target proportion of the hatchery origin escapement to be removed from the spawning grounds (in order to
#'  ensure a high proportion of NOS). The effective removal is discounted by the mark rate, i.e., `premove_HOS * m`.
#'  The removed hatchery origin fish do not spawn and are not available for broodtake. A value less than one can represent
#'  imperfect implementation of weir removal. Default is zero.
#' @slot fec_brood Vector of length `maxage` or an array with dimension `[nsim, maxage, nyears+proyears]`. The fecundity schedule of broodtake to calculate the total egg production for the hatchery. If missing, uses `Bio@fec`.
#' @slot fitness_type Character vector length 2. The fitness function to apply in the natural and hatchery environment, respectively. For each, either "Ford" or "none".
#' @slot theta Vector length 2. The optimum phenotype value for the natural and hatchery environments.
#' @slot rel_loss Vector length 3. The loss in fitness apportioned among the egg, fry, and smolt life stages which reduces survival. Theoretically, the three values should sum to 1. Alternatively,
#' set to zero to set fitness loss to zero for that specific life stage (survival is one).
#' @slot zbar_start Vector length 2. The mean phenotype value in the natural and hatchery populations at the start of the projection. Alternatively,
#' an array by dimension `[nsim, maxage, 2]`, where the age slot corresponds to cohort.
#' @slot fitness_variance Numeric. The variance (omega-squared) of the fitness function. Assumed identical between the natural and hatchery environments. Default is 100.
#' @slot phenotype_variance Numeric. The variance (sigma-squared) of the phenotypic trait (theta). Assumed identical between the natural and hatchery environments. Default is 10.
#' @slot heritability Numeric or vector length `[nsim]`. The heritability (h-squared) of the phenotypic trait. Between 0-1. Default is 0.5
#' @slot fitness_floor Numeric. The minimum fitness value in the natural and hatchery environments, i.e., fitness cannot drop below this threshold. Default is 0.5.

