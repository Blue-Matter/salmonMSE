
#' @rdname salmonMSE-int
#' @param NOS Logical, whether the Stock or Fleet object corresponds to natural origin or hatchery origin fish
#' @param stage Character indicating the corresponding salmon life stage of the Stock or Fleet object
#' @return
#' `make_Stock`: List containing a \linkS4class{Stock} object and accompanying custom parameters list
#' @export
make_Stock <- function(SOM, NOS = TRUE, stage = c("immature", "return", "escapement"), start = list()) {
  stage <- match.arg(stage)

  #### Stock object ----
  Stock <- new("Stock")
  cpars_bio <- list()

  Stock@Name <- paste(ifelse(NOS, "NOS:", "HOS:"), stage)

  # Fish spawn at the beginning of the terminal age + 1
  Stock@maxage <- SOM@maxage
  n_age <- Stock@maxage + 1

  # Survival to be modeled in cpars
  Stock@M <- Stock@Msd <- c(0, 0)

  if (NOS) {
    # 1 = Beverton-Holt, 2 = Ricker. Ricker is trigger for model with partial maturity
    Stock@SRrel <- ifelse(SOM@SRrel == "BH", 1, 3)

    if (SOM@SRrel == "BH") {
      SRRpars_hist <- sapply(1:SOM@nsim, function(x) {
        calc_SRRpars(SOM@prod_smolt[x], SOM@capacity_smolt[x], SOM@fec, SOM@p_female, SOM@SRrel)
      })

      phi0 <- SOM@SAR_NOS * SOM@p_female * SOM@fec

      h <- MSEtool::hconv(SRRpars_hist[1, ], phi0, SR = ifelse(SOM@SRrel == "BH", 1, 2))
      R0 <- MSEtool::R0conv(SRRpars_hist[1, ], SRRpars_hist[2, ], phi0, SR = ifelse(SOM@SRrel == "BH", 1, 2))

      Stock@h <- range(h)
      Stock@R0 <- mean(R0)
      cpars_bio$hs <- h
      cpars_bio$R0 <- R0
    } else {

      # Need custom SRR for OM with partial maturity at age
      h <- R0 <- rep(1, SOM@nsim)

      Stock@SRrel <- 3

      Stock@h <- rep(0.99, 2)
      Stock@R0 <- 1

      SRRfun <- function(SB, SRRpars) SRRpars$a * SB * exp(-SRRpars$b * SB)
      SRRpars <- data.frame(a = SOM@a, b = 1/SOM@Smax)
      relRfun <- function(SSBpR, SRRpars) log(SRRpars$a * SSBpR)/SRRpars$b/SSBpR
      SPRcrashfun <- function(SSBpR0, SRRpars) SRRpars$a
      cpars_bio$SRR <- list(SRRfun = SRRfun, SRRpars = SRRpars, relRfun = relRfun, SPRcrashfun = SPRcrashfun)
    }

  } else {
    # Custom SRR
    Stock@SRrel <- 3

    Stock@h <- rep(0.99, 2)
    Stock@R0 <- 1

    SRRfun <- function(SB, SRRpars) ifelse(sum(SB), 1, 0)
    SRRpars <- data.frame(x = 1:SOM@nsim)
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
  Stock@a <- Stock@b <- 0

  # Ignore area dynamics
  Stock@Size_area_1 <- Stock@Frac_area_1 <- Stock@Prob_staying <- c(0.5, 0.5)

  # Proportion of releases that die (only if there is discarding from mark-selective fishing)
  if (stage == "immature") {
    Stock@Fdisc <- rep(SOM@release_mort[1], 2)
  } else {
    Stock@Fdisc <- rep(SOM@release_mort[2], 2)
  }

  # Custom pars
  # Catch works on the basis of biomass so we need to set weight at age = 1
  cpars_bio$Wt_age <- array(1, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  # Spawning at the last age class, solely for necessary but irrelevant per recruit calculations
  # Maturity that activates return phase modeled by Herm feature
  # Unfished spawners per recruit must be identical between immature and mature NOS population
  cpars_bio$Mat_age <- array(0, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  if (stage == "escapement") { # All escapement are mature and will spawn
    cpars_bio$Mat_age[] <- 1
  } else {
    cpars_bio$Mat_age[, n_age, ] <- 1 # No spawning, n_age should be an empty age class
  }

  cpars_bio$M_ageArray <- array(.Machine$double.eps, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))
  # All ocean survival occurs in the age class prior to maturity
  # There is a requirement for spawners per recruit spool-up.
  # Should be okay so long as all fish immature prior to mat_age
  if (stage == "immature") {
    if (NOS) {
      cpars_bio$M_ageArray[, 2:n_age - 1, ] <- SOM@Mocean_NOS
    } else {
      cpars_bio$M_ageArray[, 2:n_age - 1, ] <- SOM@Mocean_HOS
    }
  } else if (stage == "return") {
    cpars_bio$M_ageArray[, n_age, ] <- 0.1 # Return does not experience M, n_age should be an empty age class
  } else if (stage == "escapement") {
    cpars_bio$M_ageArray[] <- 1e8  # Escapement either removed (broodtake, etc.) or spawns at beginning of year and then dies
  }

  # Unfished spawners per recruit must be identical between immature and mature NOS populations (for BH..?)
  cpars_bio$Fec_age <- array(0, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  if (stage == "escapement") {
    fec_array <- array(SOM@fec, c(SOM@maxage, SOM@nsim, SOM@nyears + SOM@proyears)) %>%
      aperm(c(2, 1, 3))
    if (NOS) {
      cpars_bio$Fec_age[, -1, ] <- SOM@p_female * fec_array
    } else {
      cpars_bio$Fec_age[, -1, ] <- SOM@p_female * SOM@gamma * fec_array
    }
  } else {
    cpars_bio$Fec_age[, n_age, ] <- 1
  }

  #cpars_bio$spawn_time_frac <- rep(0, SOM@nsim) # Spawn at the beginning of the year in the age of maturity + 1

  # No plus group
  cpars_bio$plusgroup <- 0L

  # Generate recruitment devs from Historical object. Specify future rec dev = 1, updated by hatchery Rel
  Perr_y <- matrix(0, SOM@nsim, SOM@maxage + SOM@nyears + SOM@proyears)

  if (length(SOM@HistN)) {
    if (stage == "immature") {
      pind <- ifelse(NOS, 1, 2)

      NPR_unfished <- matrix(NA_real_, SOM@nsim, SOM@maxage)
      NPR_unfished[, 1] <- 1
      for (a in 2:SOM@maxage) {
        NPR_unfished[, a] <- NPR_unfished[, a-1] * exp(-cpars_bio$M_ageArray[, a-1, 1])
      }

      # With custom SRR, we set R0 = 1!
      Ninit <- SOM@HistN[, , 1, pind]
      Perr_init <- Ninit/NPR_unfished
      Perr_y[, seq(SOM@maxage + 1, 2)] <- Perr_init
      Perr_y[, 1] <- 0

      # Main ts
      if (NOS) {
        Perr_y[, SOM@maxage + 2:SOM@nyears] <- 1
      } else {
        Perr_y[, SOM@maxage + 2:SOM@nyears] <- SOM@HistN[, 1, -1, pind]
      }

      # Projected rec dev
      Perr_y[, SOM@maxage + SOM@nyears + 1:SOM@proyears] <- 1
    }
  }

  cpars_bio$Perr_y <- Perr_y

  return(list(Stock = Stock, cpars_bio = cpars_bio))
}

#' @rdname salmonMSE-int
#' @param NOS Logical, whether the Stock or Fleet object corresponds to natural origin or hatchery origin fish
#' @param stage Character indicating the corresponding salmon life stage of the Stock or Fleet object
#' @return
#' `make_Stock`: List containing a \linkS4class{Fleet} object and accompanying custom parameters list
#' @export
make_Fleet <- function(SOM, NOS = TRUE, stage = c("immature", "return", "escapement")) {
  stage <- match.arg(stage)

  Fleet <- new("Fleet")
  Fleet@Name <- switch(
    stage,
    "immature" = "Preterminal fishery",
    "return" = "Terminal fishery",
    "escapement" = "Escapement placeholder"
  )
  Fleet@nyears <- SOM@nyears
  Fleet@CurrentYr <- Sys.Date() %>% format("%Y") %>% as.numeric()
  Fleet@EffYears <- 1:Fleet@nyears

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

  Find <- ifelse(NOS, 1, 2)
  if (stage == "immature") {
    if (length(SOM@HistFPT)) {
      F_hist <- SOM@HistFPT[, , Find]
    } else {
      F_hist <- matrix(0.1, SOM@nsim, SOM@nyears)
    }
  } else if (stage == "return") {

    if (length(SOM@HistFT)) {
      F_hist <- SOM@HistFT[, , Find]
    } else {
      F_hist <- matrix(0.1, SOM@nsim, SOM@nyears)
    }

  } else {
    F_hist <- matrix(0, SOM@nsim, SOM@nyears)
  }
  cpars_fleet$Find <- F_hist

  maxage <- SOM@maxage
  n_age <- maxage + 1
  cpars_fleet$V <- array(0, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  # Need to define selectivity in all populations (life stages) in order to remove excessive messages about selectivity < 0.01.
  if (stage == "immature") {

    if (all(!SOM@vulPT)) {
      cpars_fleet$V[, n_age, ] <- 1
    } else {
      cpars_fleet$V[, -n_age, ] <- SOM@vulPT %>%
        array(c(SOM@maxage, SOM@nsim, SOM@nyears + SOM@proyears)) %>%
        aperm(c(2, 1, 3))
    }

  } else if (stage == "return") {

    if (all(!SOM@vulT)) {
      cpars_fleet$V[, n_age, ] <- 1
    } else {
      cpars_fleet$V[, -n_age, ] <- SOM@vulT %>%
        array(c(SOM@maxage, SOM@nsim, SOM@nyears + SOM@proyears)) %>%
        aperm(c(2, 1, 3))
    }

  } else {
    cpars_fleet$V[, n_age - 1, ] <- 1
  }

  cpars_fleet$retA <- cpars_fleet$V
  if (SOM@m > 0) {
    if (NOS) {
      cpars_fleet$retA[] <- 1e-8
    } else {
      cpars_fleet$retA <- cpars_fleet$V * SOM@m
    }
  }

  return(list(Fleet = Fleet, cpars_fleet = cpars_fleet))
}

