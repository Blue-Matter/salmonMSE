
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
#' @slot SAR_HOS Vector length `nsim`. Survival to adult return of hatchery origin fish (density-independent, time-invariant).
#' @slot gamma Numeric. The relative fecundity of hatchery origin spawners relative to natural origin spawners.
#' @slot pmax_NOB Numeric. The maximum proportion of the natural origin escapement to be used as broodtake.
#' @slot ptarget_NOB Numeric. The target proportion of the natural origin broodtake relative to the overall broodtake.
#' @slot premove_HOS Numeric. The proportion of the hatchery origin spawners to be removed from the spawning grounds (in order to
#'  ensure a high proportion of NOS). This removal is also used to calculate the available hatchery origin fish for broodtake.
#' @slot fec_brood Numeric. The fecundity of broodtake.
#' @slot fitness_type Character. Fitness function for hatchery origin spawners. Either "Ford" or "none".
#' @slot theta Vector length 2. The optimum phenotype value for the natural and hatchery environments.
#' @slot rel_loss Vector length 3. The loss in fitness apportioned among the egg, fry, and smolt life stages. The three values should sum to 1.
#' @slot pbar_start Vector length 2. The mean phenotype value in the the natural and hatchery population at the start of the projection.
#' @slot fitness_variance Numeric. The variance of the phenotype in the population. Assumed identical between the natural and hatchery environments.
#' @slot selection_strength Numeric. The ratio between the fitness standard deviation and the phenotype standard deviation.
#' @slot heritability Numeric. The heritability of the phenotypic trait. Between 0-1.
#' @slot fitness_floor Numeric. The minimum fitness value in the natural environment.

