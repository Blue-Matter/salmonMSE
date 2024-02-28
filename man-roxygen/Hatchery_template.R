
#' @slot Name Identifying name.
#' @slot n_yearling Numeric. The number of yearlings released. No hatchery is modeled if `n_yearling` and
#' `n_subyearling` are zero.
#' @slot n_subyearling Numeric. The number of subyearlings released. No hatchery is modeled if `n_yearling` and
#' `n_subyearling` are zero.
#' @slot s_prespawn Vector length `nsim`. The survival of broodtake prior to egg production. Default is 1.
#' @slot s_egg_smolt Vector length `nsim`. The survival of eggs to the smolt life stage (for yearling release).
#' @slot s_egg_subyearling Vector length `nsim`. The survival of eggs to subyearling life stage (for subyearling release).
#' @slot gamma Vector length `nsim`. The relative fecundity of hatchery origin spawners relative to natural origin spawners.
#' @slot pmax_NOB Numeric. The maximum ratio of natural origin broodtake to the natural origin escapement.
#' @slot ptarget_NOB Numeric. The target proportion of the natural origin broodtake relative to the overall broodtake.
#' @slot premove_HOS Numeric. The proportion of the hatchery origin spawners removed from the spawning grounds (in order to
#'  ensure a high proportion of NOS).
#' @slot fec_brood Numeric. The fecundity of broodtake.
#' @slot fitness_type Character. Fitness function for hatchery origin spawners. Either "Ford", "Busack", or "none".
#' @slot theta Vector length 3.
#' @slot rel_loss Vector length 3.
#' @slot Zpop_start Vector length 3.
#' @slot fitness_variance Numeric.
#' @slot selection_strength Numeric.
#' @slot heritability Numeric.
#' @slot fitness_floor Numeric.

