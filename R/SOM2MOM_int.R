
#' @rdname salmonMSE-int
#' @param NOS Logical, whether the Stock or Fleet object corresponds to natural origin or hatchery origin fish
#' @param s Integer, the population integer for which to create the Stock or Fleet object
#' @param g Integer, the life history group for which to create the Stock object. Not relevant if `NOS = FALSE`
#' @param stage Character indicating the corresponding salmon life stage of the Stock or Fleet object
#' @return
#' `make_Stock`: List containing a \linkS4class{Stock} object and accompanying custom parameters list
#' @export
make_Stock <- function(SOM, s = 1, g = 1, NOS = TRUE, stage = c("immature", "return", "escapement")) {
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
  nyears <- 2 * SOM@nyears
  proyears <- 2 * SOM@proyears

  # Time steps for first and second half of the year
  t1 <- seq(1, nyears, 2)
  t2 <- seq(2, nyears, 2)

  all_t1 <- seq(1, nyears + proyears, 2)
  all_t2 <- seq(2, nyears + proyears, 2)

  a2 <- seq(2, n_age, 2)

  # Survival to be modeled in cpars
  Stock@M <- Stock@Msd <- c(0, 0)

  # Stock-recruit relationship (custom option)
  Stock@SRrel <- 3
  Stock@h <- rep(0.99, 2)
  Stock@R0 <- 1
  if (NOS) {
    cpars_bio$SRR <- make_SRR(Bio, Habitat)
  } else {
    # Custom SRR
    Stock@SRrel <- 3

    Stock@h <- rep(0.99, 2)
    Stock@R0 <- 1

    SRRpars <- data.frame(x = 1:SOM@nsim)
    SRRfun <- function(SB, SRRpars) ifelse(sum(SB), 1, 0)
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
  # openMSE catch is in biomass so we need to set weight at age = 1 to operate as numbers
  cpars_bio$Wt_age <- array(1, c(SOM@nsim, n_age, nyears + proyears))
  cpars_bio$Wa <- cpars_bio$Wb <- 1 # Avoids lm fitting in MSEtool::SampleStockPars()

  # Spawning at the last age class, solely for necessary but irrelevant per recruit calculations
  # Maturity that activates return phase modeled by Herm feature
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
      cpars_bio$M_ageArray[, a2, all_t2] <- Bio@Mjuv_NOS[, , , g]
    } else {
      cpars_bio$M_ageArray[, a2, all_t2] <- Hatchery@Mjuv_HOS
    }
  } else if (stage == "return") {
    cpars_bio$M_ageArray[, n_age, ] <- 0.1 # Return does not experience M, n_age should be an empty age class
  } else if (stage == "escapement") {
    cpars_bio$M_ageArray[] <- 1e8  # Escapement either removed (broodtake, etc.) or spawns at beginning of year and then dies
  }
  cpars_bio$M_ageArray[cpars_bio$M_ageArray < .Machine$double.eps] <- .Machine$double.eps

  # Unfished spawners per recruit must be identical between immature and mature NOS populations (?)
  cpars_bio$Fec_age <- array(0, c(SOM@nsim, n_age, nyears + proyears))
  if (stage == "escapement") {
    fec_array <- array(Bio@fec, c(Bio@maxage, SOM@nsim, SOM@nyears + SOM@proyears)) %>%
      aperm(c(2, 1, 3))
    age_fec <- seq(3, n_age, 2)
    if (NOS) {
      cpars_bio$Fec_age[, age_fec, all_t1] <- Bio@p_female * fec_array
    } else {
      cpars_bio$Fec_age[, age_fec, all_t1] <- Bio@p_female * Hatchery@gamma * fec_array
    }
  } else {
    cpars_bio$Fec_age[, n_age, ] <- 1
  }

  #cpars_bio$spawn_time_frac <- rep(0, SOM@nsim) # Spawn at the beginning of the year, i.e., same time step when age-1 smolts appear

  # No plus group
  cpars_bio$plusgroup <- 0L

  # Generate recruitment devs from Historical object. Specify future rec dev = 1, updated by hatchery Rel
  cpars_bio$Perr_y <- matrix(0, SOM@nsim, Stock@maxage + nyears + proyears)
  if (stage == "immature" && NOS) cpars_bio$Perr_y[, Stock@maxage + all_t1] <- 1

  # Escapement needs an initial abundance = 1 in order to generate the smolt production in the first historical year
  if (stage == "escapement" && NOS) cpars_bio$Perr_y[, Stock@maxage + 1] <- 1

  # Main historical ts, get deviations in fecundity (escapement) and Perr_y (immature)
  if (stage %in% c("immature", "escapement")) {
    pind <- ifelse(NOS, 1, 2)
    if (NOS) {
      HistPars <- lapply(1:SOM@nsim, function(x) {
        calc_Perr_main(
          Njuv = Historical@HistNjuv_NOS[x, , , g],
          Fjuv = outer(Harvest@vulPT, Historical@HistFPT[x, , pind]),
          Fterm = outer(Harvest@vulT, Historical@HistFT[x, , pind]),
          p_mature = Bio@p_mature[x, , ],
          nyears = SOM@nyears,
          HistSpawner = if (length(Historical@HistSpawner_NOS)) Historical@HistSpawner_NOS[x, , , g] else NULL,
          fec = Bio@fec,
          p_female = Bio@p_female,
          gamma = 1,
          NOS = TRUE,
          SRRfun = cpars_bio$SRR$SRRfun,
          SRRpars = cpars_bio$SRR$SRRpars[x, ]
        )
      })
    } else {
      HistPars <- lapply(1:SOM@nsim, function(x) {
        calc_Perr_main(
          Njuv = Historical@HistNjuv_HOS[x, , ],
          Fjuv = outer(Harvest@vulPT, Historical@HistFPT[x, , pind]),
          Fterm = outer(Harvest@vulT, Historical@HistFT[x, , pind]),
          p_mature = Bio@p_mature[x, , ],
          nyears = SOM@nyears,
          HistSpawner = if (length(Historical@HistSpawner_HOS)) Historical@HistSpawner_HOS[x, , ] else NULL,
          fec = Bio@fec,
          p_female = Bio@p_female,
          gamma = Hatchery@gamma,
          NOS = FALSE
        )
      })
    }

    # Specify rec devs
    if (stage == "immature") {
      # Initial abundance vector. With custom SRR, we set R0 = 1!
      survOM <- sapply(1:SOM@nsim, function(x) {
        calc_survival(cpars_bio$M_ageArray[x, , 1], p_mature = rep(0, Stock@maxage))
      }) # maxage x nsim
      NPR <- t(survOM)
      age_init <- seq(1, Stock@maxage, 2)[-1]
      if (NOS) {
        Ninit <- Historical@HistNjuv_NOS[, -1, 1, g]
      } else {
        Ninit <- Historical@HistNjuv_HOS[, -1, 1]
      }
      Perr_init <- Ninit/NPR[, rev(age_init)]
      cpars_bio$Perr_y[, rev(age_init)] <- Perr_init

      Perr_main <- sapply(HistPars, getElement, "dev")
      cpars_bio$Perr_y[, Stock@maxage + c(t1, max(t1) + 2)] <- t(Perr_main)

    } else if (stage == "escapement") {
      # Specify fecundity adjustment to match escapement to spawners
      fec_adjust <- replicate(length(age_fec), sapply(HistPars, getElement, "fec_ratio")) %>%
        aperm(c(2, 3, 1))
      cpars_bio$Fec_age[, age_fec, t1 + 2] <- cpars_bio$Fec_age[, age_fec, t1 + 2] * fec_adjust
    }
  }

  cpars_bio$AC <- rep(0, SOM@nsim)

  return(list(Stock = Stock, cpars_bio = cpars_bio))
}

#' @rdname salmonMSE-int
#' @return
#' `make_Stock`: List containing a \linkS4class{Fleet} object and accompanying custom parameters list
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
  nyears <- 2 * SOM@nyears
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
  Find <- ifelse(NOS, 1, 2)
  if (stage == "immature") {
    if (length(Historical@HistFPT)) {
      F_hist <- Historical@HistFPT[, , Find]
      cpars_fleet$Find[, t1] <- F_hist
      cpars_fleet$Find[, t2] <- 0
    }
  } else if (stage == "return") {
    if (length(Historical@HistFT)) {
      F_hist <- Historical@HistFT[, , Find]
      cpars_fleet$Find[, t1] <- 0
      cpars_fleet$Find[, t2] <- F_hist
    }
  }
  if (stage %in% c("immature", "return")) { # FinF must be greater than zero
    FinF <- cpars_fleet$Find[, nyears] < 1e-4
    cpars_fleet$Find[FinF, nyears] <- 1e-4
  }

  cpars_fleet$V <- array(0, c(SOM@nsim, n_age, nyears + proyears))
  if (stage == "immature") {

    if (all(!Harvest@vulPT)) {
      # Define dummy selectivity to remove excessive messages about selectivity < 0.01.
      cpars_fleet$V[, n_age, ] <- 1
    } else {
      vulPT <- array(Harvest@vulPT, c(Bio@maxage, SOM@nsim, nyears + proyears)) %>%
        aperm(c(2, 1, 3))
      cpars_fleet$V[, a1, ] <- vulPT
    }

  } else if (stage == "return") {

    if (all(!Harvest@vulT)) {
      cpars_fleet$V[, n_age, ] <- 1
    } else {
      vulT <- array(Harvest@vulT, c(Bio@maxage, SOM@nsim, nyears + proyears)) %>%
        aperm(c(2, 1, 3))
      cpars_fleet$V[, a2, ] <- vulT
    }

  } else {
    cpars_fleet$V[, n_age, ] <- 1
  }

  if (Harvest@MSF) {
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


calc_survival <- function(Mjuv, p_mature, maxage = length(Mjuv)) {
  surv <- rep(NA_real_, maxage)
  surv[1] <- 1
  for (a in 2:maxage) {
    surv[a] <- surv[a-1] * exp(-Mjuv[a-1]) * (1 - p_mature[a-1])
  }
  return(surv)
}

make_Stock_objects <- function(SOM, s = 1) {
  do_hatchery <- SOM@Hatchery[[s]]@n_yearling > 0 || SOM@Hatchery[[s]]@n_subyearling > 0
  has_strays <- any(SOM@stray[-s, s] > 0)

  Stocks <- lapply(1:SOM@Bio[[s]]@n_g, function(g) {
    S <- list()
    S[[1]] <- make_Stock(SOM, s, g, NOS = TRUE, stage = "immature")   # NOS_juv
    S[[2]] <- make_Stock(SOM, s, g, NOS = TRUE, stage = "return")     # NOS_recruitment
    S[[3]] <- make_Stock(SOM, s, g, NOS = TRUE, stage = "escapement") # NOS_escapement
    return(S)
  })
  Stocks <- do.call(c, Stocks)
  np <- length(Stocks)

  if (do_hatchery || has_strays) {
    Stocks[[np + 1]] <- make_Stock(SOM, s, NOS = FALSE, stage = "immature")   # HOS_juv
    Stocks[[np + 2]] <- make_Stock(SOM, s, NOS = FALSE, stage = "return")     # HOS_recruitment
    Stocks[[np + 3]] <- make_Stock(SOM, s, NOS = FALSE, stage = "escapement") # HOS_escapement
  }
  return(Stocks)
}

make_Fleet_objects <- function(SOM, s = 1) {
  do_hatchery = SOM@Hatchery[[s]]@n_yearling > 0 || SOM@Hatchery[[s]]@n_subyearling > 0
  has_strays <- any(SOM@stray[-s, s] > 0)

  Fleets <- lapply(1:SOM@Bio[[s]]@n_g, function(g) {
    FF <- list()
    FF[[1]] <- make_Fleet(SOM, s, NOS = TRUE, stage = "immature")   # NOS_juv
    FF[[2]] <- make_Fleet(SOM, s, NOS = TRUE, stage = "return")     # NOS_recruitment
    FF[[3]] <- make_Fleet(SOM, s, NOS = TRUE, stage = "escapement") # NOS_escapement
    return(FF)
  })
  Fleets <- do.call(c, Fleets)
  np <- length(Fleets)

  if (do_hatchery || has_strays) {
    Fleets[[np + 1]] <- make_Fleet(SOM, s, NOS = FALSE, stage = "immature")   # HOS_juv
    Fleets[[np + 2]] <- make_Fleet(SOM, s, NOS = FALSE, stage = "return")     # HOS_recruitment
    Fleets[[np + 3]] <- make_Fleet(SOM, s, NOS = FALSE, stage = "escapement") # HOS_escapement
  }
  return(Fleets)
}

make_stock_index <- function(SOM, check = FALSE) {
  if (check) SOM <- check_SOM(SOM)
  ns <- length(SOM@Bio)
  n_g <- sapply(1:ns, function(s) SOM@Bio[[s]]@n_g)
  np_s <- sapply(1:ns, function(s) {
    do_hatchery <- SOM@Hatchery[[s]]@n_yearling > 0 || SOM@Hatchery[[s]]@n_subyearling > 0
    has_strays <- any(SOM@stray[-s, s] > 0)

    3 * SOM@Bio[[s]]@n_g + ifelse(do_hatchery || has_strays, 3, 0)
  })
  .make_stock_index(ns, n_g, np_s)
}

.make_stock_index <- function(ns, n_g, np_s) {

  dat <- lapply(1:ns, function(s) {

    g <- rep(1:n_g[s], each = 3)
    n_natural <- n_g[s] * 3
    n_hatchery <- np_s[s] - n_natural

    df_natural <- data.frame(
      s = s,
      g = g,
      origin = "natural",
      stage = rep(c("juvenile", "recruitment", "escapement"), times = n_g[s])
    )

    if (n_hatchery > 0) {
      df_hatchery <- data.frame(
        s = s,
        g = NA,
        origin = "hatchery",
        stage = rep(c("juvenile", "recruitment", "escapement"), times = n_hatchery/3)
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

calc_Perr_main <- function(Njuv, Fjuv, Fterm, p_mature, nyears, HistSpawner, fec, p_female, gamma = 1, NOS = TRUE,
                           SRRfun, SRRpars) {

  Recruit <- Njuv[, 1:nyears] * exp(-Fjuv) * p_mature[, 1:nyears]
  SpawnerOM <- Recruit * exp(-Fterm) # Escapement

  if (!length(HistSpawner)) HistSpawner <- SpawnerOM

  EggOM <- colSums(SpawnerOM * gamma * fec * p_female)
  EggHist <- colSums(HistSpawner * gamma * fec * p_female)

  fec_ratio <- EggHist/EggOM
  fec_ratio[is.na(fec_ratio)] <- 1 # No spawners

  if (NOS) {
    Smolt_pred <- SRRfun(EggHist, SRRpars)
    Smolt_actual <- Njuv[1, ]
    dev <- Smolt_actual/c(1, Smolt_pred) # Length nyears + 1
  } else {
    dev <- Njuv[1, ]
  }
  dev[is.na(dev)] <- 0
  list(dev = dev, fec_ratio = fec_ratio)
}
