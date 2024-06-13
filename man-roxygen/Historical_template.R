
#' @slot HistSmolt Vector by historical years (`nyears`) or a matrix by `[nsim, nyears]`. The natural smolt production.
#' @slot HistSpawner Array by `[nsim, maxage, nyears]`. Total spawners (natural and hatchery origin combined) in the natural environment.
#' @slot HistN Array by `[nsim, maxage, nyears, 2]`. The abundance of immature fish at the beginning of the annual time step. The first array slice
#' corresponds to natural origin fish and the second array slice corresponds to hatchery origin fish.
#' @slot HistYearling Vector by historical years. Hatchery releases (number of smolts).
#' @slot HistFPT Vector by historical years (`nyears`) or an array by dimension `[nsim, nyears, 2]`. The instantaneous fishing mortality in the preterminal fishery.
#' The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish.
#' @slot HistFT Vector by historical years (`nyears`) or an array by dimension `[nsim, nyears, 2]`. The instantaneous fishing mortality in the terminal fishery.
#' The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish.
