
#' Reference points
#'
#' @description
#' Calculate MSY and Sgen reference points for the operating model. Uses the biological parameters (maturity, natural mortality)
#' in the last year of the projection.
#' - `calc_MSY()` calculates the MSY reference points from a set of biological and fishery parameters
#' - `calc_Sgen()` calculates the Sgen, the spawner abundance that would reach the spawner abundance at MSY after one generation without fishing
#' - `calc_ref()` is a wrapper function that calculates MSY and Sgen for an operating model
#'
#' @param SOM An object of class [salmonMSE::SOM-class]
#' @param rel_F Numeric length 2, indicates the relative effort in the preterminal and terminal fisheries, with a maximum value of 1.
#' The default is `c(0, 1)` which indicates a yield calculation with only the terminal fishery.
#' @param check Logical, whether to check the SOM object using [check_SOM()]
#' @param maximize Character, whether the MSY calculation is the optimum that maximizes catch (`"MSY"`) or excess recruitment (`"MER"`). The two
#' methods should be equivalent when `rel_F = c(0, 1)`.
#' @returns
#' - `calc_MSY` returns a vector of various state variables (catch, exploitation rate, egg production, spawners) at MSY
#' - `calc_Sgen` returns a numeric
#' - `calc_ref` returns a list by stock, each containing a matrix of MSY state variables and Sgen by simulation
#' @seealso [calc_Smsy_Ricker()]
#' @importFrom stats optimize
#' @export
calc_ref <- function(SOM, rel_F, check = TRUE, maximize = c("MSY", "MER")) {
  maximize <- match.arg(maximize)

  if (check) SOM <- check_SOM(SOM)

  ns <- length(SOM@Bio)

  y <- SOM@nyears + SOM@proyears

  if (missing(rel_F)) rel_F <- c(0, 1)
  rel_F_s <- lapply(1:ns, function(...) rel_F)

  ref_s <- lapply(1:ns, function(s) {

    if (SOM@Habitat[[s]]@use_habitat) return(matrix(0, 0, 0))

    SRR <- make_SRR(SOM@Bio[[s]])

    Mjuv_NOS <- SOM@Bio[[s]]@Mjuv_NOS
    p_female <- SOM@Bio[[s]]@p_female
    vulPT <- SOM@Harvest[[s]]@vulPT
    vulT <- SOM@Harvest[[s]]@vulT
    maxage <- SOM@Bio[[s]]@maxage
    s_enroute <- SOM@Bio[[s]]@s_enroute
    n_g <- SOM@Bio[[s]]@n_g
    p_LHG <- SOM@Bio[[s]]@p_LHG

    fec <- matrix(SOM@Bio[[s]]@fec, maxage, n_g)

    val <- sapply(1:SOM@nsim, function(x) {

      Mjuv <- matrix(Mjuv_NOS[x, , y, ], maxage, n_g)
      p_mature <- matrix(SOM@Bio[[s]]@p_mature[x, , y], maxage, n_g)

      ref <- calc_MSY(Mjuv, fec, p_female, rel_F_s[[s]], vulPT[x, ], vulT[x, ], p_mature, s_enroute, n_g, p_LHG,
                      SRR$SRRpars[x, ], maximize)

      if (any(ref < 0, na.rm = TRUE)) {

        ref[!is.na(ref)] <- 0
        ref["Sgen"] <- 0

      } else {
        Sgen <- calc_Sgen(Mjuv, fec, p_female, rel_F_s[[s]], vulPT[x, ], vulT[x, ], p_mature, s_enroute, n_g, p_LHG,
                          SRR$SRRpars[x, ], SMSY = ref["Spawners_MSY"], maximize)
        ref["Sgen"] <- Sgen
      }
      return(ref)
    })

    return(val)
  })

  return(ref_s)
}

#' @rdname calc_ref
#' @param Mjuv Numeric `maxage` for juvenile natural mortality. Can be a matrix `[maxage, n_g]`.
#' @param fec Numeric `maxage` for fecundity. Can be a matrix `[maxage, n_g]`.
#' @param p_female Numeric for proportion female spawners
#' @param vulPT Numeric `maxage` for preterminal vulnerability at age
#' @param vulT Numeric `maxage` for terminal vulnerability at age
#' @param p_mature Numeric `maxage` for maturity proportions at age. Can be a matrix `[maxage, n_g]`.
#' @param s_enroute Numeric for en-route survival of escapement to spawning grounds
#' @param n_g Integer, number of life history groups within a cohort
#' @param p_LHG Numeric `n_g` for proportion of the total egg production assigned to each life history group within a cohort
#' @param SRRpars Data frame, one row, that contains the stock recruit parameters that predicts density-dependent survival at the egg-smolt life stage
#' @param F_search Numeric, length 2 for the range of F values to search for the instantaneous fishing mortality that produces MSY
#' @export
calc_MSY <- function(Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature, s_enroute, n_g = 1, p_LHG = 1,
                     SRRpars, maximize = c("MSY", "MER"), F_search = c(1e-8, 5)) {
  maximize <- match.arg(maximize)

  if (!is.matrix(Mjuv)) Mjuv <- matrix(Mjuv, length(Mjuv), n_g)
  if (!is.matrix(p_mature)) p_mature <- matrix(p_mature, length(p_mature), n_g)
  if (!is.matrix(fec)) fec <- matrix(fec, length(fec), n_g)

  if (missing(p_LHG)) p_LHG <- rep(1/n_g, n_g)

  opt <- optimize(
    .calc_eq, interval = F_search,
    Mjuv = Mjuv,
    fec = fec,
    p_female = p_female, rel_F = rel_F,
    vulPT = vulPT, vulT = vulT,
    p_mature = p_mature,
    s_enroute = s_enroute,
    n_g = n_g,
    p_LHG = p_LHG,
    SRRpars = SRRpars,
    maximize = maximize,
    maximum = TRUE
  )

  ref <- .calc_eq(
    .F = opt$maximum,
    Mjuv = Mjuv, fec = fec, p_female = p_female, rel_F = rel_F,
    vulPT = vulPT, vulT = vulT, p_mature = p_mature, s_enroute = s_enroute,
    n_g = n_g, p_LHG = p_LHG,
    SRRpars = SRRpars, maximize = maximize, opt = FALSE
  ) %>% unlist()

  names(ref) <- paste0(names(ref), "_MSY")
  if (ref["KPT_MSY"] == 0) ref["UPT_MSY"] <- 0
  if (ref["KT_MSY"] == 0) ref["UT_MSY"] <- 0

  return(ref)
}

#' @rdname calc_ref
#' @param SMSY Numeric, spawning abundance at MSY
#' @param nyears Integer, number of years to project the population with no fishing to reach `SMSY`. Default
#' is the minimum age of maturity.
#' @export
calc_Sgen <- function(Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature, s_enroute, n_g = 1, p_LHG = 1,
                      SRRpars, SMSY, F_search = c(1e-8, 100), nyears) {

  if (!is.matrix(Mjuv)) Mjuv <- matrix(Mjuv, length(Mjuv), n_g)
  if (!is.matrix(p_mature)) p_mature <- matrix(p_mature, length(p_mature), n_g)
  if (!is.matrix(fec)) fec <- matrix(fec, length(fec), n_g)

  if (missing(p_LHG)) p_LHG <- rep(1/n_g, n_g)
  if (missing(nyears)) nyears <- which(rowSums(p_mature) > 0)[1]

  opt_Sgen <- try(
    uniroot(
      .calc_Sgen, interval = F_search,
      Mjuv = Mjuv, fec = fec, p_female = p_female, rel_F = rel_F,
      vulPT = vulPT, vulT = vulT, p_mature = p_mature, s_enroute = s_enroute,
      n_g = n_g, p_LHG = p_LHG,
      SRRpars = SRRpars, SMSY = SMSY, nyears = nyears
    ),
    silent = TRUE
  )

  if (is.character(opt_Sgen)) {
    Sgen <- NA
  } else {
    Sgen_out <- .calc_Sgen(
      .F = opt_Sgen$root,
      Mjuv = Mjuv, fec = fec, p_female = p_female, rel_F = rel_F,
      vulPT = vulPT, vulT = vulT, p_mature = p_mature, s_enroute = s_enroute,
      n_g = n_g, p_LHG = p_LHG,
      SRRpars = SRRpars, SMSY = SMSY, nyears = nyears, opt = FALSE
    )
    Sgen <- sum(Sgen_out$Sgen$Spawners)
    attr(Sgen, "Sgen") <- Sgen_out$Sgen
    attr(Sgen, "proj") <- Sgen_out$proj
  }

  return(Sgen)
}

.calc_Sgen <- function(.F, Mjuv, fec, p_female, gamma = 1, rel_F, vulPT, vulT, p_mature, s_enroute = 1,
                       n_g, p_LHG, SRRpars, popt = 1, SMSY, nyears, opt = TRUE) {

  # Equilibrium quantities at .F
  Sgen_ref <- .calc_eq(
    .F, Mjuv, fec, p_female, gamma, rel_F, vulPT, vulT, p_mature, s_enroute, n_g, p_LHG,
    SRRpars, opt = FALSE, aggregate_age = FALSE
  )

  if (Sgen_ref$Egg <= 0) {
    Sgen_plusone <- 0
  } else {
    maxage <- length(fec)
    if (missing(nyears)) nyears <- maxage

    Njuv <- Return <- Spawner <- array(0, c(maxage, nyears, n_g))
    Egg <- rep(0, nyears)
    Njuv[, 1, ] <- Sgen_ref$Njuv

    for (y in 1:nyears) {
      Return[, y, ] <- Njuv[, y, ] * p_mature * s_enroute
      Spawner[, y, ] <- Return[, y, ] * s_enroute

      Egg[y] <- sum(Spawner[, y, ] * gamma * p_female * fec) # sum across LHG

      if (y < nyears) {
        Smolt <- calc_smolt(
          Egg[y],
          kappa = SRRpars["kappa"], capacity = SRRpars["capacity"], Smax = SRRpars["Smax"], phi = SRRpars["phi"],
          fitness_loss = 1, SRrel = as.character(SRRpars["SRrel"])
        ) %>%
          as.numeric()
        Njuv[1, y+1, ] <- p_LHG * Smolt
        Njuv[2:maxage, y+1, ] <- Njuv[2:maxage - 1, y, ] * exp(-Mjuv[2:maxage - 1, ]) * (1 - p_mature[2:maxage - 1, ])
      }
    }

    Sgen_plusone <- sum(Spawner[, nyears, ])
  }

  if (opt) {
    return(Sgen_plusone - popt * SMSY)
  } else {
    return(list(Sgen = Sgen_ref, proj = list(Njuv = Njuv, Return = Return, Spawner = Spawner, Egg = Egg)))
  }
}


.calc_eq <- function(.F, Mjuv, fec, p_female, gamma = 1, rel_F, vulPT, vulT, p_mature, s_enroute = 1, n_g, p_LHG,
                     SRRpars, maximize = c("MSY", "MER"), opt = TRUE, aggregate_age = TRUE) {

  maximize <- match.arg(maximize)

  if (maximize == "MER") rel_F <- c(0, 1)

  FPT <- vulPT * .F * rel_F[1]
  FT <- vulT * .F * rel_F[2]

  surv_juv <- sapply(1:n_g, function(g) p_LHG[g] * calc_survival(Mjuv[, g] + FPT, p_mature[, g])) # First semester due to exploitation
  surv_return <- surv_juv * p_mature
  surv_esc <- surv_return * exp(-FT)
  surv_spawn <- surv_esc * s_enroute

  EPR <- sum(surv_spawn * gamma * p_female * fec) # Egg per smolt

  Smolt <- calc_smolt(
    EPR,
    kappa = SRRpars["kappa"], capacity = SRRpars["capacity"], Smax = SRRpars["Smax"], phi = SRRpars["phi"],
    fitness_loss = 1, SRrel = as.character(SRRpars["SRrel"]), per_recruit = TRUE
  ) %>%
    as.numeric()

  Njuv <- Smolt * surv_juv
  Return <- Smolt * surv_return

  KPT <- Njuv * (1 - exp(-FPT))
  KT <- Return * (1 - exp(-FT))

  Spawners <- Smolt * surv_spawn

  if (opt) {

    if (maximize == "MSY") {
      return(sum(KPT, KT))
    } else {
      return(sum(Return) - sum(Spawners))
    }


  } else {

    output <- list(
      KPT = KPT,
      KT = KT,
      FPT = max(FPT),
      FT = max(FT),
      UPT = max(1 - exp(-FPT)),
      UT = max(1 - exp(-FT)),
      Smolt = Smolt,
      Njuv = Njuv,
      Return = Return,
      Escapement = Smolt * surv_esc,
      Spawners = Spawners,
      Egg = Smolt * EPR
    )

    if (aggregate_age) output <- lapply(output, sum)

    return(output)
  }
}


#' Ricker reference points
#'
#' Compute reference points (Umsy, Smsy, and Sgen) from Ricker stock-recruit function
#' based on Scheuerell (2016).
#'
#' @param loga Numeric, alpha parameter (returns per spawner) in the Ricker function: \eqn{R=S\exp(\log(a)-bS)}
#' where `S` is the number of spawners and `R` is the return
#' @param b Numeric, beta parameter
#'
#' @importFrom gsl lambert_W0
#'
#' @returns All three functions return a numeric
#' @references Scheuerell, M.D. 2016. An explicit solution for
#' calculating optimum spawning stock size from Ricker’s stock recruitment model. PeerJ 4:e1623. \doi{10.7717/peerj.1623}
#' @seealso [calc_ref()]
#' @export
calc_Smsy_Ricker <- function(loga,b) {
  #gives the min Ricker log-likelihood
  Smsy <- (1 - gsl::lambert_W0(exp(1 - loga))) /b

  return(Smsy)
}

#' @rdname calc_Smsy_Ricker
#' @export
calc_Umsy_Ricker <- function(loga) {
  #gives the min Ricker log-likelihood
  umsy <- (1 - gsl::lambert_W0(exp(1 - loga)))

  return(umsy)
}



#' @rdname calc_Smsy_Ricker
#' @export
calc_Sgen_Ricker <- function(loga, b){
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  a <- exp(loga)

  return(-1/b*gsl::lambert_W0(-b*sMSY/a))
}

