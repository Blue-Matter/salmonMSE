
#' @slot HistSpawner_NOS Array by `[nsim, maxage, nyears, n_g]`. Natural origin spawners at age in the natural environment. *Optional*
#' @slot HistSpawner_HOS Array by `[nsim, maxage, nyears, n_r]`. Hatchery origin spawners at age in the natural environment. *Optional*
#' @slot HistNjuv_NOS Array by `[nsim, maxage, nyears+1, n_g]`. The abundance of immature natural origin fish at the beginning of the annual time step.
#' Default assumes 1000 smolts (age-1) fish annually.
#' @slot HistNjuv_HOS Array by `[nsim, maxage, nyears+1, n_r]`. The abundance of immature hatchery origin fish at the beginning of the annual time step.
#' Default assumes 1000 smolts (age-1) fish annually.
#' @slot HistFPT Vector by historical years (`nyears`) or an array by dimension `[nsim, nyears, 2]`. The instantaneous fishing mortality in the preterminal fishery.
#' The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish. Default is zero.
#' @slot HistFT Vector by historical years (`nyears`) or an array by dimension `[nsim, nyears, 2]`. The instantaneous fishing mortality in the terminal fishery.
#' The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish. Default is zero.
