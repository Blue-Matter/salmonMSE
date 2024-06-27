
#' @slot HistSpawner Array by `[nsim, maxage, nyears, 2]`. Spawners at age in the natural environment. For the fourth dimension, the first array corresponds
#' to natural origin spawners and the second array corresponds to hatchery origin spawners. *Optional*
#' @slot HistN Array by `[nsim, maxage, nyears, 2]`. The abundance of immature fish at the beginning of the annual time step. For the fourth dimension, the first array corresponds
#' to natural origin fish and the second array corresponds to hatchery origin fish. Default assumes 1000 smolts (age-1) fish annually.
#' @slot HistFPT Vector by historical years (`nyears`) or an array by dimension `[nsim, nyears, 2]`. The instantaneous fishing mortality in the preterminal fishery.
#' The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish. Default is zero.
#' @slot HistFT Vector by historical years (`nyears`) or an array by dimension `[nsim, nyears, 2]`. The instantaneous fishing mortality in the terminal fishery.
#' The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish. Default is zero.
