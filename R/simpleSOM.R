
# ---- simpleSOM Class -----
#' Class \code{"simpleSOM"}
#'
#' A simple operating model for modeling recruit-spawners and a terminal marine fishery, no hatchery production.
#'
#' Various parameters can be stochastic (length `nsim`) or input as a single numeric
#' (value identical across all simulations).
#'
#' @name simpleSOM-class
#' @docType class
#' @slot Name Character. Identifying name
#' @slot maxage Integer. The maximum age of the population age structure.
#' @slot nsim Integer. Number of simulations.
#' @slot ngen Integer. The number of generations (life cycles) to run the projection.
#' @slot kappa Vector length `nsim`. The adult productivity ratio for the stock-recruit function. **Units of recruits per spawner.**
#'  Natural per-capita production of recruits as the population approaches zero (density-independent component).
#' @slot Smax Vector length `nsim`. The spawner abundance that maximizes smolt production in the Ricker stock-recruit function. **Units of spawners.**
#' @slot sigmaR Vector length `nsim`. The lognormal standard deviation in the Ricker stock-recruit relationship. Recruitment anomalies
#'  will be simulated with `[stats::rlnorm()]`.
#' @slot Recdev Optional matrix `[nsim, ngen-1]` of recruitment anomalies in the Ricker stock-recruit relationship. Overrides the `sigmaR` argument
#'  and provides more user flexibility, for example, simulating with auto-correlation.
#' @slot type_T Character. Whether to manage terminal fishery catch from exploitation rate ("u") or catch target ("catch").
#' Default is "u".
#' @slot u_terminal Numeric, matrix `[nsim, proyears]`, or function. If `type_T = "u"`, the harvest rate (ratio of kept catch to
#' of the terminal marine fishery. Function should be of the form `function(NO, HO, m) return(u)`.
#' @slot K_T Numeric or function. If `type_T = "catch"`, the catch target of the return in the terminal fishery.
#' Function should be of the form `function(NO, HO, m) return(K)`.
#' #' @slot InitReturn Vector `[nsim]`. The return at the beginning of the projection. Default assumes 1000.
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("simpleSOM")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("simpleSOM")
simpleSOM <- setClass(
  "simpleSOM",
  slots = c(
    Name = "character",
    nsim = "numeric",
    maxage = "numeric",
    ngen = "numeric",
    kappa = "numeric",
    Smax = "numeric",
    sigmaR = "numeric",
    Recdev = "matrix",
    type_T = "character",
    u_terminal = "num.matrix.function",
    K_T = "num.function",
    InitReturn = "num.array"
  )
)

check_simpleSOM <- function(simpleSOM) {

  simpleSOM <- check_numeric(simpleSOM, "nsim")
  simpleSOM <- check_numeric(simpleSOM, "maxage")
  simpleSOM <- check_numeric(simpleSOM, "ngen")

  nsim <- simpleSOM@nsim

  simpleSOM <- check_numeric2nsim(simpleSOM, "kappa", nsim)
  simpleSOM <- check_numeric2nsim(simpleSOM, "Smax", nsim)

  if (length(simpleSOM@Recdev)) {
    simpleSOM <- check_array(simpleSOM, "Recdev", c(nsim, simpleSOM@ngen-1))
  } else {
    simpleSOM <- check_numeric2nsim(simpleSOM, "sigmaR", nsim)
  }

  simpleSOM <- check_numeric(simpleSOM, "type_T", default = "u")
  simpleSOM@type_T <- match.arg(simpleSOM@type_T, choices = c("u", "catch"))

  if (simpleSOM@type_T == "u") {
    if (length(simpleSOM@u_terminal) == 1) {
      simpleSOM <- check_numeric(simpleSOM, "u_terminal", default = 0)
    } else if (is.matrix(simpleSOM@u_terminal)) {
      simpleSOM <- check_array(simpleSOM, "u_terminal", c(nsim, proyears))
    }
    simpleSOM <- check_numeric(simpleSOM, "K_T", default = NA_real_)
  } else {
    simpleSOM <- check_numeric(simpleSOM, "u_terminal", default = NA_real_)
    simpleSOM <- check_numeric(simpleSOM, "K_T", default = 0)
  }

  if (!length(simpleSOM@InitReturn)) simpleSOM@InitReturn <- 1000
  simpleSOM <- check_numeric2nsim(simpleSOM, "InitReturn", nsim)

  return(simpleSOM)
}

#' @name salmonMSE
#' @description `simple_salmonMSE()` is helper function that converts a simple operating model (modeling only recruit-spawners
#' with Ricker relationship, terminal marine fishery, without hatchery production) to a full operating model and runs the projection.
#' @param simpleSOM An object of class \linkS4class{simpleSOM}
#' @param ... Other arguments to pass to `salmonMSE()`
#' @export
simple_salmonMSE <- function(simpleSOM, ...) {

  simpleSOM <- check_simpleSOM(simpleSOM)

  nsim <- simpleSOM@nsim
  maxage <- simpleSOM@maxage
  ngen <- simpleSOM@ngen

  if (maxage < 2) stop("maxage should be 2 or greater")
  proyears <- maxage * (ngen - 1) + 1

  if (length(simpleSOM@Recdev)) {
    Recdev <- simpleSOM@Recdev
  } else {
    sigmaR <- simpleSOM@sigmaR
    Recdev <- sapply(1:nsim, function(x) {
      rlnorm(ngen-1, -0.5 * sigmaR[x]^2, sigmaR[x])
    }) |> t()
  }

  # Recruitment anomalies
  Mjuv <- array(0, c(nsim, maxage-1, proyears))
  Mjuv[, 1, seq(2, proyears, maxage)] <- -log(Recdev)

  p_mature <- c(rep(0, maxage - 1), 1)
  Bio <- new(
    "Bio",
    maxage = maxage,
    p_mature = p_mature,
    SRrel = "Ricker",
    kappa = simpleSOM@kappa,
    Smax = simpleSOM@Smax,
    phi = rep(1, nsim),
    tau = rep(1, nsim),
    Mjuv_NOS = Mjuv,
    fec = p_mature,
    p_female = 1,
    s_enroute = 1
  )

  Habitat <- new(
    "Habitat",
    use_habitat = FALSE
  )

  Hatchery <- new(
    "Hatchery",
    n_yearling = 0,
    n_subyearling = 0
  )

  Harvest <- new(
    "Harvest",
    type_T = simpleSOM@type_T,
    u_preterminal = 0,
    u_terminal = simpleSOM@u_terminal,
    K_PT = 0,
    K_T = simpleSOM@K_T,
    MSF_PT = FALSE,
    MSF_T = FALSE,
    vulPT = rep(0, maxage),
    vulT = rep(1, maxage)
  )

  InitNjuv_NOS <- array(0, c(nsim, maxage, 1))
  InitNjuv_NOS[, maxage, ] <- simpleSOM@InitReturn
  Historical <- new("Historical", InitNjuv_NOS = InitNjuv_NOS)

  SOM <- new(
    "SOM",
    nsim = nsim,
    proyears = proyears,
    Bio = Bio,
    Habitat = Habitat,
    Hatchery = Hatchery,
    Harvest = Harvest,
    Historical = Historical
  )

  SMSE <- salmonMSE(SOM, ...)
  return(SMSE)
}

