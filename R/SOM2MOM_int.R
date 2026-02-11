
#' @rdname salmonMSE-int
#' @param NOS Logical, whether the Stock or Fleet object corresponds to natural origin or hatchery origin fish
#' @param s Integer, the population integer for which to create the Stock or Fleet object
#' @param g Integer, the life history group for which to create the Stock object. Not relevant if `NOS = FALSE`
#' @param r Integer, the hatchery release group for which to create the Stock object. Not relevant if `NOS = TRUE`
#' @param stage Character indicating the corresponding salmon life stage of the Stock or Fleet object
#' @return
#' `make_Stock`: List containing a [MSEtool::Stock-class] object and accompanying custom parameters list
#' @export
make_Stock <- function(SOM, s = 1, g = 1, r = 1, NOS = TRUE, stage = c("immature", "return", "escapement")) {
  stage <- match.arg(stage)

  Bio <- SOM@Bio[[s]]
  Habitat <- SOM@Habitat[[s]]
  Harvest <- SOM@Harvest[[s]]
  Hatchery <- SOM@Hatchery[[s]]
  Historical <- SOM@Historical[[s]]

  #### Stock object ----
  Stock <- new("Stock")
  cpars_bio <- list()

  Stock@Name <- paste("Population", s, ifelse(NOS, "NOS:", "HOS:"), stage)

  # Fish spawn at the beginning of the terminal age + 1
  n_age <- 2 * Bio@maxage + 1
  Stock@maxage <- n_age - 1
  nyears_real <- 2
  nyears <- 2 * nyears_real
  proyears <- 2 * SOM@proyears

  # Time steps for first and second half of the year
  t1 <- seq(1, nyears, 2)
  t2 <- seq(2, nyears, 2)

  t1_proj <- nyears + seq(1, proyears, 2)
  t2_proj <- nyears + seq(2, proyears, 2)

  all_t1 <- seq(1, nyears + proyears, 2)
  all_t2 <- seq(2, nyears + proyears, 2)

  a2 <- seq(2, n_age, 2)

  # Arbitrary placeholder, survival to be modeled in cpars
  Stock@M <- Stock@Msd <- c(0, 0)

  # Stock-recruit relationship (custom option)
  Stock@SRrel <- 3
  Stock@h <- rep(0.99, 2) # Arbitrary placeholder, not used
  Stock@R0 <- 1 # Arbitrary placeholder, not used

  # Custom SRR
  if (NOS && !Habitat@use_habitat) {
    cpars_bio$SRR <- make_SRR(Bio)
  } else {
    SRRpars <- data.frame(x = 1:SOM@nsim)
    SRRfun <- function(SB, SRRpars) if(sum(SB, na.rm = TRUE)) rep(1, length(SB)) else 0
    relRfun <- function(SSBpR, SRRpars) return(1)
    SPRcrashfun <- function(SSBpR0, SRRpars) return(0)
    cpars_bio$SRR <- list(SRRfun = SRRfun, SRRpars = SRRpars, relRfun = relRfun, SPRcrashfun = SPRcrashfun)
  }

  # Ignore depletion setup (use cpars$qs = 1, and Eff)
  Stock@D <- c(0, 0)

  # Rec devs to be updated in MICE rel
  Stock@Perr <- Stock@AC <- c(0, 0)

  # Ignore size dynamics
  Stock@Linf <- c(1, 1)
  Stock@K <- c(0.1, 0.1)
  Stock@t0 <- c(0, 0)
  Stock@LenCV <- c(0.05, 0.05)
  Stock@Linfsd <- Stock@Ksd <-
    Stock@L50 <- Stock@L50_95 <- c(0, 0)
  Stock@a <- Stock@b <- 1

  # Ignore area dynamics
  Stock@Size_area_1 <- Stock@Frac_area_1 <- Stock@Prob_staying <- c(0.5, 0.5)

  # Proportion of releases that die (only if there is discarding from mark-selective fishing)
  if (stage == "immature") {
    Stock@Fdisc <- rep(Harvest@release_mort[1], 2)
  } else {
    Stock@Fdisc <- rep(Harvest@release_mort[2], 2)
  }

  # Custom pars
  # openMSE catch is in biomass so we need to set weight at age = 1 so that catch is in numbers, i.e., pieces
  cpars_bio$Wt_age <- array(1, c(SOM@nsim, n_age, nyears + proyears))
  cpars_bio$Wa <- cpars_bio$Wb <- 1 # Arbitrary placeholders, avoids lm fitting in MSEtool::SampleStockPars()

  # Spawning at the last age class, solely for necessary but irrelevant per recruit calculations
  # Proportion maturity that activates return phase is actually modeled by Herm feature
  # Unfished spawners per recruit must be identical between immature and mature NOS population
  cpars_bio$Mat_age <- array(0, c(SOM@nsim, n_age, nyears + proyears))
  if (stage == "escapement") { # All escapement are mature and will spawn in the beginning of the first half of the year
    cpars_bio$Mat_age[, , all_t1] <- 1
  } else {
    cpars_bio$Mat_age[, n_age, ] <- 1 # No spawning, n_age should be an empty age class
  }

  # Ocean survival occurs in the second time step of the year
  cpars_bio$M_ageArray <- array(0, c(SOM@nsim, n_age, nyears + proyears))
  if (stage == "immature") {
    if (NOS) {
      cpars_bio$M_ageArray[, a2[-length(a2)], t2_proj] <- Bio@Mjuv_NOS[, , , g]
    } else {
      cpars_bio$M_ageArray[, a2[-length(a2)], t2_proj] <- Hatchery@Mjuv_HOS[, , , r]
    }
  } else if (stage == "return") {
    cpars_bio$M_ageArray[, n_age, ] <- 0.1 # Return does not experience M, n_age should be an empty age class
  } else if (stage == "escapement") {
    cpars_bio$M_ageArray[, , nyears + seq(1, proyears)] <- 1e8  # Escapement either removed (broodtake, etc.) or spawns at beginning of year and then dies
  }
  cpars_bio$M_ageArray[cpars_bio$M_ageArray < .Machine$double.eps] <- .Machine$double.eps

  # Unfished spawners per recruit must be identical between immature and mature NOS populations (?)
  cpars_bio$Fec_age <- array(0, c(SOM@nsim, n_age, nyears + proyears))
  if (stage == "escapement") {
    age_fec <- seq(3, n_age, 2) # We have to shift the age up once integer age class due to delay in the
    p_female <- array(Bio@p_female, c(Bio@maxage, SOM@nsim, SOM@proyears)) %>% aperm(c(2, 1, 3))
    if (NOS) {
      cpars_bio$Fec_age[, age_fec, t1_proj] <- p_female * Bio@fec
    } else {
      cpars_bio$Fec_age[, age_fec, t1_proj] <- p_female * Hatchery@gamma * Bio@fec
    }
  } else {
    cpars_bio$Fec_age[, n_age, ] <- 1
  }

  #cpars_bio$spawn_time_frac <- rep(0, SOM@nsim) # Spawn at the beginning of the year, i.e., same time step when age-1 smolts appear

  # No plus group - salmon are terminal spawners
  cpars_bio$plusgroup <- 0L

  # Generate recruitment devs from Historical object. Specify future rec dev = 1, updated by hatchery Rel
  cpars_bio$Perr_y <- matrix(0, SOM@nsim, Stock@maxage + nyears + proyears)
  #if (stage == "immature" && NOS) cpars_bio$Perr_y[, Stock@maxage + all_t1] <- 1
  if (stage == "immature" && NOS) cpars_bio$Perr_y[, Stock@maxage + t1_proj] <- 1

  # Escapement needs an initial abundance = 1 in order to generate the smolt production in the first historical year
  #if (stage == "escapement" && NOS) cpars_bio$Perr_y[, Stock@maxage + 1] <- 1

  # Main historical ts, get deviations in fecundity (escapement) and Perr_y (immature)
  #if (stage %in% c("immature", "escapement")) {
  #  pind <- ifelse(NOS, 1, 2)
  #  if (NOS) {
  #    HistPars <- lapply(1:SOM@nsim, function(x) {
  #      calc_Perr_main(
  #        Njuv = Historical@HistNjuv_NOS[x, , , g],
  #        Fjuv = outer(Harvest@vulPT[x, ], Historical@HistFPT[x, , pind]),
  #        Fterm = outer(Harvest@vulT[x, ], Historical@HistFT[x, , pind]),
  #        p_mature = Bio@p_mature[x, , ],
  #        nyears = SOM@nyears,
  #        HistSpawner = if (is.array(Historical@HistSpawner_NOS) && length(Historical@HistSpawner_NOS)) Historical@HistSpawner_NOS[x, , , g] else NULL,
  #        fec = Bio@fec[x, , seq(1, SOM@nyears)],
  #        p_female = Bio@p_female,
  #        gamma = 1,
  #        NOS = TRUE,
  #        SRRfun = cpars_bio$SRR$SRRfun,
  #        SRRpars = cpars_bio$SRR$SRRpars[x, ]
  #      )
  #    })
  #  } else {
  #    HistPars <- lapply(1:SOM@nsim, function(x) {
  #      calc_Perr_main(
  #        Njuv = Historical@HistNjuv_HOS[x, , , r],
  #        Fjuv = outer(Harvest@vulPT[x, ], Historical@HistFPT[x, , pind]),
  #        Fterm = outer(Harvest@vulT[x, ], Historical@HistFT[x, , pind]),
  #        p_mature = Hatchery@p_mature_HOS[x, , , r],
  #        nyears = SOM@nyears,
  #        HistSpawner = if (is.array(Historical@HistSpawner_HOS) && length(Historical@HistSpawner_HOS)) Historical@HistSpawner_HOS[x, , , r] else NULL,
  #        fec = Bio@fec[x, , seq(1, SOM@nyears)],
  #        p_female = Bio@p_female,
  #        gamma = Hatchery@gamma,
  #        NOS = FALSE
  #      )
  #    })
  #  }
#
  #  # Specify rec devs
  #  if (stage == "immature") {
  #    # Initial abundance vector. With custom SRR, we set R0 = 1!
  #    survOM <- sapply(1:SOM@nsim, function(x) {
  #      calc_survival(cpars_bio$M_ageArray[x, , 1], p_mature = rep(0, Stock@maxage))
  #    }) # maxage x nsim
  #    NPR <- t(survOM)
  #    age_init <- seq(1, Stock@maxage, 2)[-1]
  #    if (NOS) {
  #      Ninit <- Historical@HistNjuv_NOS[, -1, 1, g]
  #    } else {
  #      Ninit <- Historical@HistNjuv_HOS[, -1, 1, r]
  #    }
  #    Perr_init <- Ninit/NPR[, rev(age_init)]
  #    cpars_bio$Perr_y[, rev(age_init)] <- Perr_init
#
  #    Perr_main <- sapply(HistPars, getElement, "dev")
  #    cpars_bio$Perr_y[, Stock@maxage + c(t1, max(t1) + 2)] <- t(Perr_main)
#
  #  } else if (stage == "escapement") {
  #    # Specify fecundity adjustment to match escapement to spawners
  #    fec_adjust <- replicate(length(age_fec), sapply(HistPars, getElement, "fec_ratio")) %>%
  #      aperm(c(2, 3, 1))
  #    cpars_bio$Fec_age[, age_fec, t1 + 2] <- cpars_bio$Fec_age[, age_fec, t1 + 2] * fec_adjust
  #  }
  #}

  cpars_bio$AC <- rep(0, SOM@nsim)

  return(list(Stock = Stock, cpars_bio = cpars_bio))
}

#' @rdname salmonMSE-int
#' @return
#' `make_Stock`: List containing a [MSEtool::Fleet-class] object and accompanying custom parameters list
#' @export
make_Fleet <- function(SOM, s, NOS = TRUE, stage = c("immature", "return", "escapement")) {
  stage <- match.arg(stage)

  Bio <- SOM@Bio[[s]]
  Harvest <- SOM@Harvest[[s]]
  Hatchery <- SOM@Hatchery[[s]]
  Historical <- SOM@Historical[[s]]

  Fleet <- new("Fleet")
  Fleet@Name <- switch(
    stage,
    "immature" = "Preterminal fishery",
    "return" = "Terminal fishery",
    "escapement" = "Escapement placeholder"
  )

  n_age <- 2 * Bio@maxage + 1
  maxage <- n_age - 1
  nyears_real <- 2
  nyears <- 2 * nyears_real
  proyears <- 2 * SOM@proyears

  Fleet@nyears <- nyears
  Fleet@CurrentYr <- nyears
  Fleet@EffYears <- 1:Fleet@nyears

  # Time step for the first and second half of the year
  t1 <- seq(1, nyears, 2)
  t2 <- seq(2, nyears, 2)

  a1 <- seq(1, 2 * Bio@maxage, 2)
  a2 <- seq(2, n_age, 2)

  Fleet@EffLower <- c(0, 0.1)
  Fleet@EffUpper <- c(0, 0.1)

  # No change in fishing efficiency, projected F proportional to effort
  Fleet@Esd <- Fleet@qinc <- Fleet@qcv <- c(0, 0)

  # Fishery selectivity - placeholder, will use V slot
  Fleet@L5 <- c(0.7, 0.7)
  Fleet@LFS <- c(0.71, 0.71)
  Fleet@Vmaxlen <- c(1, 1)

  Fleet@LR5 <- Fleet@LFR <- c(-100, -100)
  Fleet@Rmaxlen <- c(1, 1)
  Fleet@DR <- c(0, 0) # Discard rate
  Fleet@Spat_targ <- c(1, 1)
  Fleet@MPA <- FALSE

  Fleet@isRel <- FALSE

  cpars_fleet <- list()
  cpars_fleet$qs <- rep(1, SOM@nsim)

  cpars_fleet$Find <- matrix(0, SOM@nsim, Fleet@nyears) # Escapement experiences zero F
  if (stage %in% c("immature", "return")) { # FinF must be greater than zero
    cpars_fleet$Find[, nyears] <- 1e-4
  }

  cpars_fleet$V <- array(0, c(SOM@nsim, n_age, nyears + proyears))
  if (stage == "immature") {

    if (all(!Harvest@vulPT)) {
      # Define dummy selectivity to remove excessive messages about selectivity < 0.01.
      cpars_fleet$V[, n_age, ] <- 1
    } else {
      vulPT <- array(Harvest@vulPT, c(SOM@nsim, Bio@maxage, nyears + proyears))
      cpars_fleet$V[, a1, ] <- vulPT
    }

  } else if (stage == "return") {

    if (all(!Harvest@vulT)) {
      cpars_fleet$V[, n_age, ] <- 1
    } else {
      vulT <- array(Harvest@vulT, c(SOM@nsim, Bio@maxage, nyears + proyears))
      cpars_fleet$V[, a2, ] <- vulT
    }

  } else {
    cpars_fleet$V[, n_age, ] <- 1
  }

  MSF_PT <- stage == "immature" && Harvest@MSF_PT
  MSF_T <- stage == "return" && Harvest@MSF_T

  if (MSF_PT || MSF_T) {
    cpars_fleet$retA <- cpars_fleet$V
    retAproj <- cpars_fleet$retA[, , nyears + seq(1, proyears)]
    if (NOS) {
      retAproj[] <- .Machine$double.eps
    } else {
      retAproj[retAproj > 0] <- Hatchery@m
    }
    cpars_fleet$retA[, , nyears + seq(1, proyears)] <- retAproj
  }

  return(list(Fleet = Fleet, cpars_fleet = cpars_fleet))
}


calc_survival <- function(Mjuv, p_mature, maxage = length(p_mature)) {
  surv <- rep(NA_real_, maxage)
  surv[1] <- 1
  for (a in 2:maxage) {
    surv[a] <- surv[a-1] * exp(-Mjuv[a-1]) * (1 - p_mature[a-1])
  }
  return(surv)
}


#' Calculate equilibrium quantities with life history groups
#'
#' Calculate eggs/smolt or spawners/smolt based on life history parameters (survival, maturity, fecundity)
#'
#' @param Mjuv Matrix `[maxage, n_g]`, but can be a vector if `n_g = 1`. Juvenile natural mortality
#' @param p_mature Matrix `[maxage, n_g]`, but can be a vector if `n_g = 1`. Maturity at age
#' @param p_female Numeric. Proportion female
#' @param fec Matrix `[maxage, n_g]`, but can be a vector if `n_g = 1`. Fecundity at age. Only used if `output = "egg"`
#' @param s_enroute Numeric, en-route survival of escapement to spawning grounds
#' @param n_g Integer. Number of life history groups
#' @param p_LHG Vector length `n_g` of proportion of life history groups per recruit. Default is `rep(1/n_g, n_g)`
#' @param output Character to indicate the output units, e.g., "egg" returns eggs per smolt, and "spawner" returns spawners per smolt
#' @keywords internal
#' @return Numeric, units depend on `"output"` argument
calc_phi <- function(Mjuv, p_mature, p_female, fec, s_enroute = 1, n_g = 1, p_LHG, output = c("egg", "spawner")) {
  output <- match.arg(output)

  if (n_g == 1) {
    if (!is.matrix(Mjuv)) Mjuv <- matrix(Mjuv, ncol = 1)
    if (!is.matrix(p_mature)) p_mature <- matrix(p_mature, ncol = 1)
    if (!is.matrix(fec)) fec <- matrix(fec, ncol = 1)
  }
  if (missing(p_LHG)) p_LHG <- rep(1/n_g, n_g)

  surv_juv <- sapply(1:n_g, function(g) p_LHG[g] * calc_survival(Mjuv[, g], p_mature[, g])) # age x g
  Esc <- p_mature * surv_juv # escapement per smolt

  if (output == "egg") {
    x <- p_female * Esc * s_enroute * fec
  } else {
    x <- p_female * Esc * s_enroute
  }

  return(sum(x))
}

make_Stock_objects <- function(SOM, s = 1) {
  do_hatchery <- sum(SOM@Hatchery[[s]]@n_yearling, SOM@Hatchery[[s]]@n_subyearling) > 0
  has_strays <- any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)

  Stocks_g <- lapply(1:SOM@Bio[[s]]@n_g, function(g) {
    S <- list()
    S[[1]] <- make_Stock(SOM, s, g = g, NOS = TRUE, stage = "immature")   # NOS_juv
    S[[2]] <- make_Stock(SOM, s, g = g, NOS = TRUE, stage = "return")     # NOS_recruitment
    S[[3]] <- make_Stock(SOM, s, g = g, NOS = TRUE, stage = "escapement") # NOS_escapement
    return(S)
  })
  Stocks_g <- do.call(c, Stocks_g)

  if (do_hatchery || has_strays) {
    Stocks_r <- lapply(1:SOM@Hatchery[[s]]@n_r, function(r) {
      S <- list()
      S[[1]] <- make_Stock(SOM, s, r = r, NOS = FALSE, stage = "immature")   # HOS_juv
      S[[2]] <- make_Stock(SOM, s, r = r, NOS = FALSE, stage = "return")     # HOS_recruitment
      S[[3]] <- make_Stock(SOM, s, r = r, NOS = FALSE, stage = "escapement") # HOS_escapement
      return(S)
    })
    Stocks_r <- do.call(c, Stocks_r)
  } else {
    Stocks_r <- list()
  }

  c(Stocks_g, Stocks_r)
}

make_Fleet_objects <- function(SOM, s = 1) {
  do_hatchery <- sum(SOM@Hatchery[[s]]@n_yearling, SOM@Hatchery[[s]]@n_subyearling) > 0
  has_strays <- any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)

  Fleets_g <- lapply(1:SOM@Bio[[s]]@n_g, function(g) {
    FF <- list()
    FF[[1]] <- make_Fleet(SOM, s, NOS = TRUE, stage = "immature")   # NOS_juv
    FF[[2]] <- make_Fleet(SOM, s, NOS = TRUE, stage = "return")     # NOS_recruitment
    FF[[3]] <- make_Fleet(SOM, s, NOS = TRUE, stage = "escapement") # NOS_escapement
    return(FF)
  })
  Fleets_g <- do.call(c, Fleets_g)

  if (do_hatchery || has_strays) {
    Fleets_r <- lapply(1:SOM@Hatchery[[s]]@n_r, function(r) {
      FF <- list()
      FF[[1]] <- make_Fleet(SOM, s, NOS = FALSE, stage = "immature")   # HOS_juv
      FF[[2]] <- make_Fleet(SOM, s, NOS = FALSE, stage = "return")     # HOS_recruitment
      FF[[3]] <- make_Fleet(SOM, s, NOS = FALSE, stage = "escapement") # HOS_escapement
      return(FF)
    })
    Fleets_r <- do.call(c, Fleets_r)
  } else {
    Fleets_r <- list()
  }
  c(Fleets_g, Fleets_r)
}

make_stock_index <- function(SOM, check = FALSE) {
  if (check) SOM <- check_SOM(SOM)
  ns <- length(SOM@Bio)
  n_g <- sapply(1:ns, function(s) SOM@Bio[[s]]@n_g)
  n_r <- sapply(1:ns, function(s) SOM@Hatchery[[s]]@n_r)
  np_s <- sapply(1:ns, function(s) {
    do_hatchery <- sum(SOM@Hatchery[[s]]@n_yearling, SOM@Hatchery[[s]]@n_subyearling) > 0
    has_strays <- any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)
    3 * (n_g[s] + ifelse(do_hatchery || has_strays, n_r[s], 0))
  })
  .make_stock_index(ns, n_g, n_r, np_s)
}

.make_stock_index <- function(ns, n_g, n_r, np_s) {

  dat <- lapply(1:ns, function(s) {

    g <- rep(1:n_g[s], each = 3)
    n_natural <- n_g[s] * 3
    n_hatchery <- np_s[s] - n_natural

    df_natural <- data.frame(
      s = s,
      g = g,
      r = NA,
      origin = "natural",
      stage = rep(c("juvenile", "recruitment", "escapement"), times = n_g[s])
    )

    if (n_hatchery > 0) {
      r <- rep(1:n_r[s], each = 3)
      df_hatchery <- data.frame(
        s = s,
        g = NA,
        r = r,
        origin = "hatchery",
        stage = rep(c("juvenile", "recruitment", "escapement"), times = n_r[s])
      )
    } else {
      df_hatchery <- data.frame()
    }

    rbind(df_natural, df_hatchery)

  }) %>%
    bind_rows()
  dat$p <- seq_len(nrow(dat))
  return(dat)
}

#calc_Perr_main <- function(Njuv, Fjuv, Fterm, p_mature, nyears, HistSpawner, fec, p_female, gamma = 1, NOS = TRUE,
#                           SRRfun, SRRpars) {
#
#  Recruit <- Njuv[, seq(1, nyears)] * exp(-Fjuv) * p_mature[, seq(1, nyears)]
#  SpawnerOM <- Recruit * exp(-Fterm) # Escapement
#
#  if (!length(HistSpawner)) HistSpawner <- SpawnerOM
#
#  EggOM <- colSums(SpawnerOM * gamma * fec * p_female)
#  EggHist <- colSums(HistSpawner * gamma * fec * p_female)
#
#  fec_ratio <- EggHist/EggOM
#  fec_ratio[is.na(fec_ratio)] <- 1 # No spawners
#
#  if (NOS) {
#    Smolt_pred <- SRRfun(EggHist, SRRpars)
#    Smolt_actual <- Njuv[1, ]
#    dev <- Smolt_actual/c(1, Smolt_pred) # Length nyears + 1
#  } else {
#    dev <- Njuv[1, ]
#  }
#  dev[is.na(dev)] <- 0
#  list(dev = dev, fec_ratio = fec_ratio)
#}
