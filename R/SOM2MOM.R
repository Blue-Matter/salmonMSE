
#' @export
#' @rdname salmonMSE
SOM2MOM <- function(SOM) {

  # Check SOM here
  SOM <- check_SOM(SOM)

  MOM <- suppressMessages(new("MOM"))
  MOM@Name <- SOM@Name
  MOM@proyears <- SOM@proyears
  MOM@nsim <- SOM@nsim
  MOM@interval <- SOM@proyears + 1
  MOM@pstar <- 0.5
  MOM@maxF <- 4.61 #-log(0.01) # Harvest rate can't exceed 0.99
  MOM@reps <- 1

  # Stock objects
  do_hatchery = SOM@n_yearling > 0 || SOM@n_subyearling > 0

  ns <- 1
  Stocks <- list()
  Stocks[[1]] <- make_Stock(SOM, NOS = TRUE, mature = FALSE)    # NOS_juv
  Stocks[[2]] <- make_Stock(SOM, NOS = TRUE, mature = TRUE)     # NOS_adult

  if (do_hatchery) {
    Stocks[[3]] <- make_Stock(SOM, NOS = FALSE, mature = FALSE) # HOS_juv
    Stocks[[4]] <- make_Stock(SOM, NOS = FALSE, mature = TRUE)  # HOS_adult
  }
  np <- length(Stocks)

  # Single fleet object ----
  # In AHA, there is no harvest on immature fish
  nf <- 1

  Fleets <- list()
  Fleets[[1]] <- make_Fleet(SOM, NOS = TRUE, mature = FALSE)   # NOS_juv
  Fleets[[2]] <- make_Fleet(SOM, NOS = TRUE, mature = TRUE)    # NOS_adult

  if (do_hatchery) {
    Fleets[[3]] <- make_Fleet(SOM, NOS = FALSE, mature = FALSE) # HOS_juv
    Fleets[[4]] <- make_Fleet(SOM, NOS = FALSE, mature = TRUE)  # HOS_adult
  }

  MOM@Stocks <- lapply(Stocks, getElement, "Stock")
  MOM@Fleets <- lapply(1:np, function(p) list(Fleets[[p]]$Fleet))
  MOM@Obs <- lapply(1:np, function(...) {
    lapply(1:nf, function(...) MSEtool::Generic_Obs)
  })
  MOM@Imps <- lapply(1:np, function(...) {
    lapply(1:nf, function(...) MSEtool::Perfect_Imp)
  })

  MOM@cpars <- lapply(1:np, function(p) {
    lapply(1:nf, function(...) {
      c(Stocks[[p]]$cpars_bio, Fleets[[p]]$cpars_fleet)
    })
  })

  # Combine all fisheries into a single fleet
  MOM@CatchFrac <- MOM@Allocation <- lapply(1:np, function(...) matrix(1, SOM@nsim, nf))

  # Ignore Allocation, nf = 1
  # Ignore Efactor, nf = 1
  # Ignore Complexes, ns = 1

  # Specify smolt natural production from NOS population
  if (do_hatchery) {
    SSBfrom <- matrix(0, 4, 4)
    SSBfrom[3, 2] <- 1 # Triggers SRRfun so that there is hatchery production only when there are mature fish
  } else {
    SSBfrom <- matrix(0, 2, 2)
  }
  SSBfrom[1, 2] <- 1 # Recruitment of juv NOS comes from adult NOS population

  # Move juveniles to adult
  Herm_agemat <- matrix(c(SOM@p_mature, 1), SOM@nsim, SOM@maxage + 1, byrow = TRUE)
  Herm <- list()
  Herm[[1]] <- Herm_agemat
  names(Herm) <- c("H_2_1")

  if (do_hatchery) {
    Herm[[2]] <- Herm_agemat
    names(Herm)[2] <- c("H_4_3")
  }

  MOM@SexPars <- list(SSBfrom = SSBfrom, Herm = Herm, share_par = FALSE)

  # Rel
  Rel <- list()

  # The predicted smolt production (juvenile NOS) is from the adult escapement and adult escapement.
  # Combines hatchery and habitat dynamics
  # Hatchery: Presence/absence of HOS
  # Habitat: (A) Update Perr_y of juvenile NOS from adult escapement.
  # Calculates spawners (after broodtake and weir removal) then recruiement deviation based on new productivity and stock parameter, relative to the historical one

  # Not needed if no habitat improvement and no hatchery! SRR pars in the OM should suffice in that case
  # (the escapement is the spawning output)

  SRRpars_hist <- sapply(1:SOM@nsim, function(x) {
    .AHA_SRRpars(SOM@prod_smolt[x], SOM@capacity_smolt[x], SOM@fec, SOM@p_female)
  })

  SRRpars_proj <- sapply(1:SOM@nsim, function(x) {
    .AHA_SRRpars(SOM@prod_smolt[x] * SOM@prod_smolt_improve,
                 SOM@capacity_smolt[x] * SOM@capacity_smolt_improve, SOM@fec, SOM@p_female)
  })

  habitat_change <- any(abs(SRRpars_hist - SRRpars_proj) < 1e-8)

  if (do_hatchery) {
    # Determine the local brood from the number of yearling and subyearling releases, their survival from egg stage,
    # and fecundity of broodtake (identical between natural and hatchery escapement)

    # This is a management action, cannot be stochastic
    brood_import <- 0
    brood_local <- SOM@n_yearling/(SOM@s_prespawn * SOM@p_female * SOM@fec_brood * SOM@s_egg_smolt) +
      SOM@n_subyearling/(SOM@s_prespawn * SOM@p_female * SOM@fec_brood * SOM@s_egg_subyearling) -
      brood_import
    brood_local <- round(brood_local, 3)

    # Survival of eggs in the hatchery
    p_yearling <- SOM@n_yearling/(SOM@n_yearling + SOM@n_subyearling)
    s_egg_hatchery <- SOM@s_egg_subyearling * (1 - p_yearling) + SOM@s_egg_smolt * p_yearling

  } else {
    brood_local <- 0
    s_egg_hatchery <- NA
  }

  if (do_hatchery || habitat_change) {
    fitness_args <- list()

    if (do_hatchery && SOM@fitness_type == "Ford") {

      fitness_env$Ford <- data.frame()

      fitness_args <- local({
        omega <- sqrt(SOM@fitness_variance) * SOM@selection_strength
        omega2 <- omega * omega

        list(
          rel_loss = SOM@rel_loss,
          omega2 = omega * omega,
          A = 1 - SOM@heritability * SOM@fitness_variance/(omega2 + SOM@fitness_variance),
          fitness_variance = SOM@fitness_variance,
          fitness_floor = SOM@fitness_floor,
          heritability = SOM@heritability,
          theta = SOM@theta,
          pbar_start = SOM@pbar_start
        )
      })

      # Need an age-specific MICE rel for fitness function during the marine phase

    }

    # Natural smolt production from NOS and HOS escapement
    Rel[[1]] <- makeRel_smolt(
      p_r = 1, p_natural = 2, p_hatchery = 4, output = "natural",
      ptarget_NOB = SOM@ptarget_NOB, pmax_NOB = SOM@ptarget_NOB,
      brood_local = brood_local, fec_brood = SOM@fec_brood, s_egg = s_egg_hatchery,
      premove_HOS = SOM@premove_HOS, s_prespawn = SOM@s_prespawn, # Broodtake & hatchery production
      p_female = SOM@p_female, fec = SOM@fec, gamma = SOM@gamma, # Spawning (natural production)
      SRRpars_hist, SRRpars_proj, fitness_type = SOM@fitness_type, # Spawning (natural production)
      fitness_args = fitness_args
    )

  }

  if (do_hatchery) {
    # Hatchery smolt releases from NOS and HOS escapement
    Rel[[2]] <- makeRel_smolt(
      p_r = 3, p_natural = 2, p_hatchery = 4, output = "hatchery",
      ptarget_NOB = SOM@ptarget_NOB, pmax_NOB = SOM@ptarget_NOB,
      brood_local = brood_local, fec_brood = SOM@fec_brood, s_egg = s_egg_hatchery,
      premove_HOS = SOM@premove_HOS, s_prespawn = SOM@s_prespawn,
      p_female = SOM@p_female, fec = SOM@fec, gamma = SOM@gamma
    )
  }

  MOM@Rel <- Rel
  #if (length(Rel)) MOM@cpars$control <- list(HistRel = FALSE)

  # Run simulation
  #multiHist <- SimulateMOM(MOM, parallel = FALSE)
  #MMSE <- ProjectMOM(multiHist, MPs = c("Harvest_MP", "NoHarvest_MMP"), dropHist = TRUE, parallel = FALSE)
  #
  #MMSE@multiHist <- multiHist
  #return(MMSE)
  return(MOM)
}

check_SOM <- function(SOM) {

  var_len1 <- c("nyears", "proyears", "seed", "nsim", "maxage", "fec", "p_female",
                "capacity_smolt_improve", "prod_smolt_improve",
                "n_yearling", "n_subyearling",
                "pmax_NOB", "ptarget_NOB", "premove_HOS", "fec_brood",
                "s_prespawn", "s_egg_smolt", "s_egg_subyearling", "gamma",
                "u", "m",
                "fitness_type", "fitness_variance", "selection_strength", "heritability", "fitness_floor")
  lapply(var_len1, function(i) {
    v <- slot(SOM, i)
    if (length(v) != 1) stop("Slot ", i, " must be a single numeric")
  })

  var_maxage <- "p_mature"
  lapply(var_maxage, function(i) {
    v <- slot(SOM, i)
    if (length(v) != SOM@maxage) stop("Slot ", i, " must be length ", SOM@maxage)
  })

  var_stochastic <- c("capacity_smolt", "prod_smolt", "SAR")
  for(i in var_stochastic) {
    if (length(slot(SOM, i)) == 1) slot(SOM, i) <- rep(slot(SOM, i), SOM@nsim)
    if (length(slot(SOM, i)) != SOM@nsim) stop("Slot ", i, " must be length ", SOM@nsim)
  }

  slot(SOM, "fitness_type") <- match.arg(slot(SOM, "fitness_type"), choices = c("Ford", "none"))

  return(SOM)
}
