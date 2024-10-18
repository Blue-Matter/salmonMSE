
#' @details
#' A description of the fitness parameters and hatchery dynamics is available in the
#' [online documentation](https://docs.salmonmse.com/articles/equations.html#fitness-effects-on-survival).
#'
#' @slot n_yearling Numeric. The number of yearlings released. No hatchery is modeled if both `n_yearling` and
#' `n_subyearling` are zero.
#' @slot n_subyearling Numeric. The number of subyearlings released. No hatchery is modeled if both `n_yearling` and
#' `n_subyearling` are zero.
#' @slot s_prespawn Numeric. The survival of broodtake prior to egg production. `1 - s_prespawn` is the proportion of fish not used for hatchery purposes, e.g., mortality or other resesarch purposes.
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`.
#' @slot s_egg_smolt Numeric. The survival of eggs to the smolt life stage (for yearling release).
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`.
#' @slot s_egg_subyearling Numeric. The survival of eggs to subyearling life stage (for subyearling release).
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`.
#' @slot Mjuv_HOS Either vector by age (length `maxage`) or an array with dimension `[nsim, maxage, nyears+proyears]`.
#'  Natural mortality of immature hatchery origin fish.
#'  To replicate the SAR parameter of a stage-specific model, set `Mjuv_HOS[a] = -log(SAR)` for the age class prior to maturation (and zero for all other ages).
#' @slot gamma Numeric. The relative reproductive success of hatchery origin spawners (relative to natural origin spawners).
#' @slot m Numeric. The mark rate of hatchery origin fish, which affects selective broodtake and fishery retention if mark-selective fishing is utilized.
#' @slot pmax_esc Numeric. The maximum proportion of total escapement (after en route mortality) that could be used as broodtake. Set to 1 for AHA compatibility.
#' @slot pmax_NOB Numeric. The maximum proportion of the natural origin escapement (after en route mortality and `pmax_esc`) to be used as broodtake. If broodstock
#' is limited by `pmax_esc < 1`, then this parameter should be 1.
#' @slot ptarget_NOB Numeric. The target proportion of the natural origin broodtake relative to the overall broodtake.
#' @slot phatchery Numeric. The proportion of the hatchery origin escapement that return to the hatchery, for example, by removal from spawning grounds
#'  or swim-in facilities. These fish are available for broodtake.
#' @slot premove_HOS Numeric. The proportion of the hatchery origin escapement to be removed from the spawning grounds (in order to
#'  ensure a high proportion of NOS). These fish are not available for broodtake. For example, a value less than one can represent
#'  imperfect implementation of weir removal.
#' @slot fec_brood Vector of length `maxage`. The fecundity schedule of broodtake to calculate the total egg production for the hatchery.
#' @slot fitness_type Character vector length 2. The fitness function to apply in the natural and hatchery environment, respectively. For each, either "Ford" or "none".
#' @slot theta Vector length 2. The optimum phenotype value for the natural and hatchery environments.
#' @slot rel_loss Vector length 3. The loss in fitness apportioned among the egg, fry, and smolt life stages. The three values should sum to 1.
#' @slot zbar_start Vector length 2. The mean phenotype value in the natural and hatchery populations at the start of the projection.
#' @slot fitness_variance Numeric. The variance of the phenotypic trait. Assumed identical between the natural and hatchery environments.
#' @slot selection_strength Numeric. The ratio between the fitness standard deviation and the phenotype standard deviation.
#' @slot heritability Numeric. The heritability of the phenotypic trait. Between 0-1.
#' @slot fitness_floor Numeric. The minimum fitness value in the natural and hatchery environments.

