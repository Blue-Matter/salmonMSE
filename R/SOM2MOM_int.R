

make_Stock <- function(SOM, NOS = TRUE, mature = TRUE) {

  #### Stock object ----
  Stock <- new("Stock")
  cpars_bio <- list()

  if (NOS) {
    if (mature) {
      Stock@Name <- "Adult NOS"
    } else {
      Stock@Name <- "Juvenile NOS"
    }
  } else {
    if (mature) {
      Stock@Name <- "Adult HOS"
    } else {
      Stock@Name <- "Juvenile HOS"
    }
  }

  # Two sub-annual steps needed for return year (Harvest in first half, then broodtake before spawning)
  Stock@maxage <- SOM@maxage
  n_age <- Stock@maxage + 1

  # Survival to be modeled in cpars
  Stock@M <- Stock@Msd <- c(0, 0)

  Stock@h <- rep(0.25, 2)
  Stock@R0 <- 1

  if (NOS) {
    # Beverton-Holt SRR
    Stock@SRrel <- 1

    SRRpars_hist <- sapply(1:SOM@nsim, function(x) {
      .AHA_SRRpars(SOM@prod_smolt[x], SOM@capacity_smolt[x], SOM@fec, SOM@p_female)
    })

    phi0 <- SOM@SAR * SOM@p_female * SOM@fec

    h <- MSEtool::hconv(SRRpars_hist["alpha", ], phi0)
    R0 <- MSEtool::R0conv(SRRpars_hist["alpha", ], SRRpars_hist["beta", ], phi0)

    Stock@h <- range(h)
    Stock@R0 <- mean(R0)
    cpars_bio$hs <- h
    cpars_bio$R0 <- R0
  } else {
    # Custom SRR
    Stock@SRrel <- 3
    SRRfun <- function(SB, SRRpars) ifelse(sum(SB), 1, 0)
    SRRpars <- data.frame(x = 1:SOM@nsim)
    relRfun <- function(SSBpR, SRRpars) return(1)
    SPRcrashfun <- function(SSBpR0, SRRpars) return(0)
    cpars_bio$SRR <- list(SRRfun = SRRfun, SRRpars = SRRpars, relRfun = relRfun, SPRcrashfun = SPRcrashfun)
  }

  # Ignore depletion setup (use cpars$qs = 1, and Eff)
  Stock@D <- c(0, 0)

  # Rec devs to be specified in MICE rel
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
  Stock@Fdisc <- c(0, 0)

  # Custom pars
  # Catch works on the basis of biomass so we need to set weight at age = 1
  cpars_bio$Wt_age <- array(1, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  cpars_bio$Mat_age <-
    cpars_bio$Fec_age <-
    cpars_bio$M_ageArray <- array(0, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  # All ocean survival occurs in the age class prior to maturity
  # There is a requirement for spawners per recruit spool-up.
  # Should be okay so long as all fish immature prior to mat_age
  cpars_bio$M_ageArray[] <- .Machine$double.eps
  cpars_bio$M_ageArray[, n_age - 2, ] <- -log(SOM@SAR)

  # No plus group, terminal spawning at the last age class
  cpars_bio$plusgroup <- 0L

  # Spawning at the last age class (maturity that activates return phase modeled by Herm feature)
  # Unfished spawners per recruit must be identical between immature and mature NOS populations
  if (NOS) {
    cpars_bio$Fec_age[, n_age, ] <- SOM@p_female * SOM@fec # Last age class
  } else {
    cpars_bio$Fec_age[, n_age, ] <- SOM@p_female * SOM@gamma * SOM@fec # Last age class
  }
  cpars_bio$Mat_age[, n_age, ] <- 1
  #cpars_bio$spawn_time_frac <- rep(0, SOM@nsim) # Spawn at the beginning of the year in the age of maturity + 1

  # Generate recruitment every generation
  Perr_y <- matrix(0, SOM@nsim, SOM@nyears + SOM@proyears)

  if (NOS && !mature) {
    # Hatchery population has no recruitment unless Perr_y is specified by MICE Rel
    # Only works for generation by generation
    age_mat <- which(SOM@p_mature > 0)[1]
    Perr_y[, seq(1, SOM@nyears + SOM@proyears, by = age_mat)] <- 1
  }

  cpars_bio$Perr_y <- cbind(matrix(0, SOM@nsim, n_age - 1), Perr_y)

  return(list(Stock = Stock, cpars_bio = cpars_bio))
}


make_Fleet <- function(SOM, NOS = TRUE, mature = TRUE) {

  Fleet <- new("Fleet")
  Fleet@Name <- "Fishery"
  Fleet@nyears <- SOM@nyears
  Fleet@CurrentYr <- Sys.Date() %>% format("%Y") %>% as.numeric()
  Fleet@EffYears <- 1:Fleet@nyears
  Fleet@EffLower <- c(0, 0.1)
  Fleet@EffUpper <- c(0, 0.1)

  # No change in fishing efficiency, projected F proportional to effort
  Fleet@Esd <- Fleet@qinc <- Fleet@qcv <- c(0, 0)

  # Fishery selectivity
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

  maxage <- SOM@maxage
  n_age <- maxage + 1
  cpars_fleet$V <- array(0, c(SOM@nsim, n_age, SOM@nyears + SOM@proyears))

  # Need to select for all populations in order to remove excessive messages about selectivity < 0.01.
  # Should be fine if all fish are mature before fishery
  cpars_fleet$V[, n_age - 1, ] <- 1

  # Future feature for mark rate affecting harvest rate or retention of hatchery origin vs natural origin fish
  #if (NOS) {
  #  cpars_fleet$V[, 2 * age_mat, ] <- Harvest$m
  #} else {
  #  cpars_fleet$V[, 2 * age_mat, ] <- 1
  #}
  return(list(Fleet = Fleet, cpars_fleet = cpars_fleet))
}

