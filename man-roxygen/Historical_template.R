
#' @slot InitNOS Natural origin spawners at age at the start of the projection. Either a numeric to specify the total natural origin spawners (in the oldest age class)
#'   otherwise, an array by `[nsim, maxage, n_g]`. Default is 1,000 spawners in the oldest age class.
#' @slot InitHOS Hatchery origin spawners at age at the start of the projection. Either a numeric to specify the total hatchery origin spawners (in the oldest age class)
#'   otherwise, an array by `[nsim, maxage, n_r]`.
#'   Default is 1,000 spawners in the oldest age class if there is hatchery production or zero otherwise.
#' @slot InitNjuv_NOS Array by `[nsim, maxage-1, n_g]`. The abundance of immature natural origin fish at the beginning of the projection (first age class in array represents age 2).
#' Default assumes zero, which creates a population with single brood year returns.
#' @slot InitNjuv_HOS Array by `[nsim, maxage-1, n_r]`. The abundance of immature hatchery origin fish at the beginning of the projection (first age class in array represents age 2.
#' Default assumes zero, which creates a population with single brood year returns.
