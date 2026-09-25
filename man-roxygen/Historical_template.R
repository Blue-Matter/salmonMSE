

#' @slot InitNjuv_NOS Numeric. The abundance of immature natural-origin fish at the beginning of the projection.
#' Default assumes 1000 in the oldest age class, which creates a population with single brood year returns.
#' Alternatively, provide a full array by `[nsim, maxage, n_g]` for multiple-brood year returns.
#' @slot InitNjuv_HOS Numeric. The abundance of immature hatchery origin fish at the beginning of the projection.
#' Default assumes 1000 in the oldest age class, which creates a population with single brood year returns.
#' Alternatively, provide a full array by `[nsim, maxage, n_r]` for multiple-brood year returns.
#' @slot InitEsc_NOS *Optional* Numeric length 1 or `nsim`.
#' The natural-origin escapement (oldest age class) in the first year of the projection,
#' to be used in lieu of `InitNjuv_NOS`. Currently only for single brood year returns.
#' @slot InitEsc_HOS *Optional* Numeric length 1 or `nsim`.
#' The hatchery-origin escapement (oldest age class) in the first year of the projection,
#' to be used in lieu of `InitNjuv_HOS`. Currently only for single brood year returns.
