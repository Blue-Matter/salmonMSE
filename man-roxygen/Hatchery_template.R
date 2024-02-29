
#' @slot n_yearling Numeric. The number of yearlings released. No hatchery is modeled if both `n_yearling` and
#' `n_subyearling` are zero.
#' @slot n_subyearling Numeric. The number of subyearlings released. No hatchery is modeled if both `n_yearling` and
#' `n_subyearling` are zero.
#' @slot s_prespawn Numeric. The survival of broodtake prior to egg production.
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`.
#' @slot s_egg_smolt Numeric. The survival of eggs to the smolt life stage (for yearling release).
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`.
#' @slot s_egg_subyearling Numeric. The survival of eggs to subyearling life stage (for subyearling release).
#' Used to back-calculate the broodtake from `n_yearling` and `n_subyearling`.
#' @slot gamma Numeric. The relative fecundity of hatchery origin spawners relative to natural origin spawners.
#' @slot pmax_NOB Numeric. The maximum ratio of natural origin broodtake to the natural origin escapement.
#' @slot ptarget_NOB Numeric. The target proportion of the natural origin broodtake relative to the overall broodtake.
#' @slot premove_HOS Numeric. The proportion of the hatchery origin spawners removed from the spawning grounds (in order to
#'  ensure a high proportion of NOS).
#' @slot fec_brood Numeric. The fecundity of broodtake.
#' @slot fitness_type Character. Fitness function for hatchery origin spawners. Either "Ford" or "none".
#' @slot theta Vector length 3.
#' @slot rel_loss Vector length 3.
#' @slot pbar_start Vector length 2.
#' @slot fitness_variance Numeric.
#' @slot selection_strength Numeric.
#' @slot heritability Numeric.
#' @slot fitness_floor Numeric.

