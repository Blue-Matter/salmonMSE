
#' @name salmonMSE-int
#' @title Internal salmonMSE functions for converting operating model inputs and outputs
#'
#' @description
#' - [SOM2MOM()] converts a salmon operating model (\linkS4class{SOM}) to a multi-stock operating model (\linkS4class{MOM})
#' - [MMSE2SMSE()] converts the openMSE output, along with additional state variables recorded in [salmonMSE_env], into a salmon MSE object (\linkS4class{SMSE})
#' - [make_Stock()] creates the \linkS4class{Stock} object (openMSE) corresponding to salmon life stage
#' - [make_Fleet()] creates the \linkS4class{Fleet} object (openMSE) corresponding to the fishery that interacts with the various salmon life stages
#'
#' [salmonMSE()] is the wrapper function that coordinates the simulation and the output.
#'
#' @param SOM An object of class \linkS4class{SOM}
#' @param start An optional named list to specify the natural smolt production ('Natural') and smolt releases ('Hatchery') at the start of the
#' simulation
#' @export
#' @examples
#' \dontrun{
#' # One thousand natural and hatchery smolts in the first year
#' SOM <- SOM2MOM(SOM, start = list(Natural = 1e3, Hatchery = 1e3))
#' }
#' @return
#' `SOM2MOM`: \linkS4class{MOM} object
SOM2MOM <- function(SOM, start = list()) {

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
  Stocks[[1]] <- make_Stock(SOM, NOS = TRUE, stage = "immature", start)   # NOS_juv
  Stocks[[2]] <- make_Stock(SOM, NOS = TRUE, stage = "return")            # NOS_recruitment
  Stocks[[3]] <- make_Stock(SOM, NOS = TRUE, stage = "escapement")        # NOS_escapement

  if (do_hatchery) {
    Stocks[[4]] <- make_Stock(SOM, NOS = FALSE, stage = "immature", start)   # HOS_juv
    Stocks[[5]] <- make_Stock(SOM, NOS = FALSE, stage = "return")            # HOS_recruitment
    Stocks[[6]] <- make_Stock(SOM, NOS = FALSE, stage = "escapement")        # HOS_escapement
  }
  np <- length(Stocks)

  # Single fleet object ----
  # In AHA, there is no harvest on immature fish
  nf <- 1

  Fleets <- list()
  Fleets[[1]] <- make_Fleet(SOM, NOS = TRUE, stage = "immature")   # NOS_juv
  Fleets[[2]] <- make_Fleet(SOM, NOS = TRUE, stage = "return")     # NOS_recruitment
  Fleets[[3]] <- make_Fleet(SOM, NOS = TRUE, stage = "escapement") # NOS_escapement

  if (do_hatchery) {
    Fleets[[4]] <- make_Fleet(SOM, NOS = FALSE, stage = "immature")   # HOS_juv
    Fleets[[5]] <- make_Fleet(SOM, NOS = FALSE, stage = "return")     # HOS_recruitment
    Fleets[[6]] <- make_Fleet(SOM, NOS = FALSE, stage = "escapement") # HOS_escapement
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
  # Integrated hatchery: triggers SRRfun so that there is hatchery production (p = 4) only when there is natural escapement (p = 3)
  SSBfrom <- matrix(0, np, np)
  if (do_hatchery) SSBfrom[4, 3] <- 1
  SSBfrom[1, 3] <- 1 # Juv NOS predicted from NOS escapement

  # Move stocks around
  # First generate the escapement, then move juveniles to adult
  first_mature_age <- sapply(1:(SOM@nyears + SOM@proyears), function(y) {
    sapply(1:SOM@nsim, function(x) {
      which(SOM@p_mature[x, , y] > 0)[1]
    })
  })
  Herm_escapement <- sapply(1:SOM@nsim, function(x) {
    sapply(1:(SOM@nyears + SOM@proyears), function(y) ifelse(0:SOM@maxage >= first_mature_age[x, y], 1, 0))
  }, simplify = "array") %>%
    aperm(c(3, 1, 2))

  Herm_mature <- sapply(1:SOM@nsim, function(x) {
    sapply(1:(SOM@nyears + SOM@proyears), function(y) c(SOM@p_mature[x, , y], 1))
  }, simplify = "array") %>%
    aperm(c(3, 1, 2))

  # For NOS
  Herm <- list(Herm_escapement, Herm_mature)
  names(Herm) <- c("H_3_2", "H_2_1")

  # HOS
  if (do_hatchery) {
    Herm <- c(Herm, Herm)
    names(Herm)[3:4] <- c("H_6_5", "H_5_4")
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
  SRRpars_hist <- Stocks[[1]]$cpars_bio$SRR$SRRpars
  SRRpars_proj <- local({
    if (SOM@SRrel == "BH") {
      pars <- lapply(1:SOM@nsim, function(x) {
        calc_SRRpars(SOM@prod_smolt[x] * SOM@prod_smolt_improve,
                     SOM@capacity_smolt[x] * SOM@capacity_smolt_improve, sum(SOM@fec), SOM@p_female)
      })
      a <- vapply(pars, getElement, numeric(1), 1)
      b <- vapply(pars, getElement, numeric(1), 2)
    } else {
      a <- SOM@a * SOM@prod_smolt_improve
      b <- 1/SOM@Smax/SOM@capacity_smolt_improve
    }
    data.frame(a = a, b = b, phi0 = phi0, SPRcrash = 1/a/phi0, SRrel = SOM@SRrel)
  })

  #habitat_change <- any(abs(SRRpars_hist[, 1:2] - SRRpars_proj[, 1:2]) > 1e-8)
  habitat_change <- SOM@prod_smolt_improve != 1 || SOM@capacity_smolt_improve != 1

  if (do_hatchery) {
    # Determine the total number of eggs needed from the number of yearling and subyearling releases, their survival from egg stage
    # and fecundity of broodtake (identical between natural and hatchery escapement)
    # This is a management action, cannot be stochastic
    # No imports
    egg_local <- SOM@n_yearling/SOM@s_egg_smolt + SOM@n_subyearling/SOM@s_egg_subyearling

    # Survival of eggs in the hatchery
    p_yearling <- SOM@n_yearling/(SOM@n_yearling + SOM@n_subyearling)
    s_egg_hatchery <- SOM@s_egg_subyearling * (1 - p_yearling) + SOM@s_egg_smolt * p_yearling

  } else {
    egg_local <- 0
    s_egg_hatchery <- NA
  }

  if (do_hatchery || habitat_change) {
    fitness_args <- list()

    if (do_hatchery && SOM@fitness_type == "Ford") {
      fitness_args <- local({
        omega <- sqrt(SOM@fitness_variance) * SOM@selection_strength
        omega2 <- omega * omega

        list(
          rel_loss = SOM@rel_loss,
          omega2 = omega2,
          A = 1 - SOM@heritability * SOM@fitness_variance/(omega2 + SOM@fitness_variance),
          fitness_variance = SOM@fitness_variance,
          fitness_floor = SOM@fitness_floor,
          heritability = SOM@heritability,
          theta = SOM@theta,
          pbar_start = SOM@pbar_start
        )
      })
    }

    # Natural smolt production from NOS and HOS escapement and habitat
    Rel[[1]] <- makeRel_smolt(
      p_smolt = 1, p_natural = 3, p_hatchery = 6, # NOTE: what to do if no hatchery?? p_hatchery needs to be undefined
      output = "natural",
      ptarget_NOB = SOM@ptarget_NOB, pmax_NOB = SOM@pmax_NOB,
      egg_local = egg_local, fec_brood = SOM@fec_brood, s_egg = s_egg_hatchery,
      phatchery = SOM@phatchery, premove_HOS = SOM@premove_HOS, s_prespawn = SOM@s_prespawn, # Broodtake & hatchery production
      p_female = SOM@p_female, fec = SOM@fec, gamma = SOM@gamma, # Spawning (natural production)
      SRRpars_hist, SRRpars_proj, SRrel = SOM@SRrel, fitness_type = SOM@fitness_type, # Spawning (natural production)
      fitness_args = fitness_args
    )

    # Marine survival of natural origin fish
    if (SOM@fitness_type != "none") {
      Rel[[2]] <- makeRel_SAR(
        p_smolt = 1, p_naturalsmolt = 1, fitness_type = SOM@fitness_type, SAR = SOM@SAR_NOS, rel_loss = SOM@rel_loss[3],
        age_mat = which(SOM@p_mature == 1)[1]
      )
    }
  }

  if (do_hatchery) {
    nRel <- length(Rel)
    # Hatchery smolt releases from NOS and HOS escapement
    Rel[[nRel + 1]] <- makeRel_smolt(
      p_smolt = 4, p_natural = 3, p_hatchery = 6, output = "hatchery",
      ptarget_NOB = SOM@ptarget_NOB, pmax_NOB = SOM@ptarget_NOB,
      brood_local = brood_local, fec_brood = SOM@fec_brood, s_egg = s_egg_hatchery,
      phatchery = SOM@phatchery, premove_HOS = SOM@premove_HOS, s_prespawn = SOM@s_prespawn,
      p_female = SOM@p_female, fec = SOM@fec, gamma = SOM@gamma
    )
  }

  MOM@Rel <- Rel
  #if (length(Rel)) MOM@cpars$control <- list(HistRel = FALSE)

  return(MOM)
}

check_SOM <- function(SOM) {

  slot(SOM, "fitness_type") <- match.arg(slot(SOM, "fitness_type"), choices = c("Ford", "none"))
  slot(SOM, "SRrel") <- match.arg(slot(SOM, "SRrel"), choices = c("BH", "Ricker"))

  # Length 1
  var_len1 <- c("nyears", "proyears", "seed", "nsim", "maxage", "p_female",
                "capacity_smolt_improve", "prod_smolt_improve",
                "n_yearling", "n_subyearling",
                "pmax_NOB", "ptarget_NOB", "phatchery", "premove_HOS",
                "s_prespawn", "s_egg_smolt", "s_egg_subyearling", "gamma",
                "u_preterminal", "u_terminal", "m",
                "fitness_type", "fitness_variance", "selection_strength", "heritability", "fitness_floor")
  lapply(var_len1, function(i) {
    v <- slot(SOM, i)
    if (length(v) != 1) stop("Slot ", i, " must be a single numeric")
  })

  # Length 2
  var_len2 <- "release_mort"
  lapply(var_len2, function(i) {
    v <- slot(SOM, i)
    if (length(v) != 2) stop("Slot ", i, " must be a vector of length 2")
  })

  # Length maxage
  var_maxage <- c("fec", "vulPT", "vulT", "fec_brood")
  lapply(var_maxage, function(i) {
    v <- slot(SOM, i)
    if (length(v) != SOM@maxage) stop("Slot ", i, " must be length ", SOM@maxage)
  })

  # Length maxage to full array
  var_dyn <- c("p_mature", "Mocean_NOS", "Mocean_HOS")
  do_hatchery <- SOM@n_subyearling > 0 || SOM@n_yearling > 0
  for(i in var_dyn) {
    if (i != "Mocean_HOS" || do_hatchery) {
      if (is.array(slot(SOM, i))) {

        dim_check <- all(dim(slot(SOM, i)) == c(SOM@nsim, SOM@maxage, SOM@nyears + SOM@proyears))
        if (!dim_check) {
          stop("Slot ", i, " must be an array of dimension ",
               paste(c(SOM@nsim, SOM@maxage, SOM@nyears + SOM@proyears), collapse = ", "))
        }

      } else {
        if (length(slot(SOM, i)) == 1) slot(SOM, i) <- rep(slot(SOM, i), SOM@maxage)
        if (length(slot(SOM, i)) != SOM@maxage) stop("Slot ", i, " must be length ", SOM@maxage)

        dyn_array <- slot(SOM, i) %>%
          array(c(SOM@maxage, SOM@nsim, SOM@nyears + SOM@proyears)) %>%
          aperm(c(2, 1, 3))
        slot(SOM, i) <- dyn_array
      }
    }
  }

  # Length nsim
  if (SOM@SRrel == "BH") {
    var_stochastic <- c("capacity_smolt", "prod_smolt")
  } else {
    var_stochastic <- c("a", "Smax")
  }
  for(i in var_stochastic) {
    if (length(slot(SOM, i)) == 1) slot(SOM, i) <- rep(slot(SOM, i), SOM@nsim)
    if (length(slot(SOM, i)) != SOM@nsim) stop("Slot ", i, " must be length ", SOM@nsim)
  }

  # Various hist objects
  # Length nyears to `nsim, nyears` matrix or `nsim, nyears, 2` array or `nsim, maxage, nyears` array or `nsim, maxage, nyears, 2`
  var_hist <- c("HistSmolt", "HistYearling", "HistSpawner", "HistFPT", "HistFT", "HistN")
  check_hist <- sapply(var_hist, function(i) length(slot(SOM, i)) > 0)
  if (all(check_hist)) {
    for(i in var_hist) {

      dim_i <- dim(slot(SOM, i))
      if (is.null(dim_i)) {

        if (i == "HistSpawner") {
          stop("Slot ", i, " must be an array of dimension ",
               paste(c(SOM@nsim, SOM@maxage, SOM@nyears), collapse = ", "))
        } else if (i == "HistN") {
          stop("Slot ", i, " must be an array of dimension ",
               paste(c(SOM@nsim, SOM@maxage, SOM@nyears, 2), collapse = ", "))
        } else {
          if (length(slot(SOM, i)) != SOM@nyears) {
            stop("Slot ", i, " must be length ", SOM@nyears)
          }

          if (i == "HistSmolt") {
            slot(SOM, i) <- matrix(slot(SOM, i), SOM@nsim, SOM@nyears, byrow = TRUE)
          } else if (i %in% c("HistFPT", "HistFT")) {
            slot(SOM, i) <- array(slot(SOM, i), c(SOM@nyears, SOM@nsim, 2)) %>%
              aperm(c(2, 1, 3))
          }
        }

      }

      dim_i <- dim(slot(SOM, i)) # Re-check
      if (i == "HistSmolt") {
        dim_check <- length(dim_i) == 2 && all(dim_i == c(SOM@nsim, SOM@nyears))

        if (!dim_check) {
          stop("Slot ", i, " must be a matrix of dimension ",
               paste(c(SOM@nsim, SOM@nyears), collapse = ", "))
        }
      } else if (i == "HistSpawner") {
        dim_check <- length(dim_i) == 3 && all(dim_i == c(SOM@nsim, SOM@maxage, SOM@nyears))

        if (!dim_check) {
          stop("Slot ", i, " must be an array of dimension ",
               paste(c(SOM@nsim, SOM@maxage, SOM@nyears), collapse = ", "))
        }
      } else if (i == "HistN") {
        dim_check <- length(dim_i) == 4 && all(dim_i == c(SOM@nsim, SOM@maxage, SOM@nyears, 2))

        if (!dim_check) {
          stop("Slot ", i, " must be an array of dimension ",
               paste(c(SOM@nsim, SOM@maxage, SOM@nyears, 2), collapse = ", "))
        }
      } else if (i %in% c("HistFPT", "HistFT")) {
        dim_check <- length(dim_i) == 3 && all(dim_i == c(SOM@nsim, SOM@nyears, 2))

        if (!dim_check) {
          stop("Slot ", i, " must be an array of dimension ",
               paste(c(SOM@nsim, SOM@nyears, 2), collapse = ", "))
        }
      }
    }

  } else if (any(check_hist)) {
    stop("All historical slots should be filled")
  }

  return(SOM)
}
