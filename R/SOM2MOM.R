
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
#' @export
#' @return
#' `SOM2MOM`: \linkS4class{MOM} object
SOM2MOM <- function(SOM) {
  SOM <- check_SOM(SOM)

  MOM <- suppressMessages(new("MOM"))
  MOM@Name <- SOM@Name
  MOM@proyears <- 2 * SOM@proyears
  MOM@nsim <- SOM@nsim
  MOM@interval <- 1
  MOM@pstar <- 0.5
  MOM@maxF <- 4.61 #-log(0.01) # Harvest rate can't exceed 0.99
  MOM@reps <- 1

  # Stock objects
  do_hatchery = SOM@n_yearling > 0 || SOM@n_subyearling > 0

  ns <- 1
  Stocks <- list()
  Stocks[[1]] <- make_Stock(SOM, NOS = TRUE, stage = "immature")   # NOS_juv
  Stocks[[2]] <- make_Stock(SOM, NOS = TRUE, stage = "return")            # NOS_recruitment
  Stocks[[3]] <- make_Stock(SOM, NOS = TRUE, stage = "escapement")        # NOS_escapement

  if (do_hatchery) {
    Stocks[[4]] <- make_Stock(SOM, NOS = FALSE, stage = "immature")   # HOS_juv
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
  nage <- 2 * SOM@maxage + 1
  nyears <- 2 * SOM@nyears
  proyears <- 2 * SOM@proyears
  Herm_escapement <- Herm_mature <- array(0, c(SOM@nsim, nage, nyears + proyears))

  # Maturation occurs for even age classes (age class 2, 4, 6) at the beginning of even time steps.
  Herm_mature[, 2 * seq(1, SOM@maxage), seq(2, nyears + proyears, 2)] <- SOM@p_mature

  # In reality, escapement occurs for even age classes (age class 2, 4, 6) at the end of even time steps.
  # In openMSE, we do escapement for odd age classes at the beginning of the subsequent odd time steps
  Herm_escapement[, 2 * seq(1, SOM@maxage) + 1, seq(1, nyears + proyears, 2)] <- sapply(1:SOM@nsim, function(x) {
    sapply(1:(SOM@nyears + SOM@proyears), function(y) ifelse(1:SOM@maxage >= first_mature_age[x, y], 1, 0))
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
  SRRpars <- Stocks[[1]]$cpars_bio$SRR$SRRpars

  habitat_change <- SOM@kappa_improve != 1 || SOM@capacity_smolt_improve != 1

  if (do_hatchery) {
    # Determine the total number of eggs needed from the number of yearling and subyearling releases, their survival from egg stage
    # and fecundity of broodtake (identical between natural and hatchery escapement)
    # This is a management action, cannot be stochastic
    # No imports
    egg_yearling <- ifelse(SOM@n_yearling > 0, SOM@n_yearling/SOM@s_egg_smolt, 0)
    egg_subyearling <- ifelse(SOM@n_subyearling > 0, SOM@n_subyearling/SOM@s_egg_subyearling, 0)
    egg_local <- egg_yearling + egg_subyearling

    p_yearling <- SOM@n_yearling/(SOM@n_yearling + SOM@n_subyearling)

  } else {
    egg_local <- 0
    p_yearling <- NA
  }

  if (do_hatchery || habitat_change) {
    fitness_args <- list()

    if (do_hatchery && any(SOM@fitness_type == "Ford")) {
      fitness_args <- local({
        omega <- sqrt(SOM@fitness_variance) * SOM@selection_strength
        omega2 <- omega * omega

        list(
          fitness_type = SOM@fitness_type,
          rel_loss = SOM@rel_loss,
          omega2 = omega2,
          fitness_variance = SOM@fitness_variance,
          fitness_floor = SOM@fitness_floor,
          heritability = SOM@heritability,
          theta = SOM@theta,
          zbar_start = SOM@zbar_start
        )
      })
    }

    # Natural smolt production from NOS and HOS escapement and habitat
    Rel[[1]] <- makeRel_smolt(
      p_smolt = 1, p_natural = 3, p_hatchery = ifelse(do_hatchery, 6, NA_real_), output = "natural",
      ptarget_NOB = SOM@ptarget_NOB, pmax_NOB = SOM@pmax_NOB,
      egg_local = egg_local, fec_brood = SOM@fec_brood,
      s_yearling = SOM@s_egg_smolt, s_subyearling = SOM@s_egg_subyearling, p_yearling = p_yearling,
      phatchery = SOM@phatchery, premove_HOS = SOM@premove_HOS, s_prespawn = SOM@s_prespawn,
      p_female = SOM@p_female, fec = SOM@fec, gamma = SOM@gamma, SRRpars = SRRpars,
      fitness_args = fitness_args
    )

    # Marine survival of natural origin fish
    if (SOM@fitness_type[1] != "none") {
      Rel[[2]] <- makeRel_SAR(
        p_smolt = 1, p_naturalsmolt = 1, envir = "natural",
        rel_loss = SOM@rel_loss[3], nyears = 2 * SOM@nyears,
        Mbase = Stocks[[1]]$cpars_bio$M_ageArray[, , 2 * SOM@nyears + seq(1, 2 * SOM@proyears)]
      )
    }

    # Marine survival of hatchery origin fish
    if (SOM@fitness_type[2] != "none") {
      nRel <- length(Rel)
      Rel[[nRel + 1]] <- makeRel_SAR(
        p_smolt = 4, p_naturalsmolt = 1, envir = "hatchery",
        rel_loss = SOM@rel_loss[3], nyears = 2 * SOM@nyears,
        Mbase = Stocks[[4]]$cpars_bio$M_ageArray[, , 2 * SOM@nyears + seq(1, 2 * SOM@proyears)]
      )
    }
  }

  if (do_hatchery) {
    nRel <- length(Rel)
    # Hatchery smolt releases from NOS and HOS escapement
    Rel[[nRel + 1]] <- makeRel_smolt(
      p_smolt = 4, p_natural = 3, p_hatchery = 6, output = "hatchery",
      ptarget_NOB = SOM@ptarget_NOB, pmax_NOB = SOM@ptarget_NOB,
      egg_local = egg_local, fec_brood = SOM@fec_brood,
      s_yearling = SOM@s_egg_smolt, s_subyearling = SOM@s_egg_subyearling, p_yearling = p_yearling,
      phatchery = SOM@phatchery, premove_HOS = SOM@premove_HOS, s_prespawn = SOM@s_prespawn,
      p_female = SOM@p_female, fec = SOM@fec, gamma = SOM@gamma, SRRpars = SRRpars,
      fitness_args = fitness_args
    )
  }

  MOM@Rel <- Rel
  MOM@cpars$control <- list(HistRel = FALSE, HermEq = FALSE)

  return(MOM)
}

check_SOM <- function(SOM) {

  slot(SOM, "SRrel") <- match.arg(slot(SOM, "SRrel"), choices = c("BH", "Ricker"))

  # Length 1
  var_len1 <- c("nyears", "proyears", "seed", "nsim", "maxage", "p_female",
                "capacity_smolt_improve", "kappa_improve",
                "n_yearling", "n_subyearling",
                "pmax_NOB", "ptarget_NOB", "phatchery", "premove_HOS",
                "s_prespawn", "s_egg_smolt", "s_egg_subyearling", "gamma",
                "u_preterminal", "u_terminal", "m",
                "fitness_variance", "selection_strength", "heritability", "fitness_floor")
  lapply(var_len1, function(i) {
    v <- slot(SOM, i)
    if (length(v) != 1) stop("Slot ", i, " must be a single numeric")
  })

  # Length 2
  var_len2 <- c("release_mort", "fitness_type")
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
  var_dyn <- c("p_mature", "Mjuv_NOS", "Mjuv_HOS")
  do_hatchery <- SOM@n_subyearling > 0 || SOM@n_yearling > 0
  for(i in var_dyn) {
    if (i != "Mjuv_HOS" || do_hatchery) {
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
    var_stochastic <- c("capacity_smolt", "kappa", "phi")
  } else {
    var_stochastic <- c("Smax", "kappa", "phi")
  }
  for(i in var_stochastic) {
    if (i == "phi") {
      if (!length(slot(SOM, i))) {

        slot(SOM, i) <- local({
          NPR <- matrix(NA_real_, SOM@nsim, SOM@maxage)
          NPR[, 1] <- 1
          for (a in 2:SOM@maxage) {
            NPR[, a] <- NPR[, a-1] * exp(-SOM@Mjuv_NOS[, a-1, 1]) * (1 - SOM@p_mature[, a-1, 1])
          }
          EPR <- NPR * SOM@p_mature[, , 1]
          colSums(t(EPR) * SOM@p_female * SOM@fec)
        })

      }
    }
    if (length(slot(SOM, i)) == 1) slot(SOM, i) <- rep(slot(SOM, i), SOM@nsim)
    if (length(slot(SOM, i)) != SOM@nsim) stop("Slot ", i, " must be length ", SOM@nsim)
  }

  # Various hist objects


  var_F <- c("HistFPT", "HistFT")
  for(i in var_F) {
    if (!length(slot(SOM, i))) slot(SOM, i) <- rep(0, SOM@nyears)

    dim_i <- dim(slot(SOM, i))
    if (is.null(dim_i)) { # Check if vector
      if (length(slot(SOM, i)) != SOM@nyears) {
        stop("Slot ", i, " must be length ", SOM@nyears)
      }
      slot(SOM, i) <- array(slot(SOM, i), c(SOM@nyears, SOM@nsim, 2)) %>%
        aperm(c(2, 1, 3))
    }

    dim_i <- dim(slot(SOM, i)) # Re-check with full dimensions
    dim_check <- length(dim_i) == 3 && all(dim_i == c(SOM@nsim, SOM@nyears, 2))

    if (!dim_check) {
      stop("Slot ", i, " must be an array of dimension ",
           paste(c(SOM@nsim, SOM@nyears, 2), collapse = ", "))
    }
  }

  # HistSpawner check
  if (length(slot(SOM, "HistSpawner"))) {
    dim_i <- dim(slot(SOM, "HistSpawner"))
    if (is.null(dim_i)) {
      stop("Slot ", i, " must be an array of dimension ",
           paste(c(SOM@nsim, SOM@maxage, SOM@nyears, 2), collapse = ", "))
    } else {
      dim_check <- length(dim_i) == 4 && all(dim_i == c(SOM@nsim, SOM@maxage, SOM@nyears, 2))
      if (!dim_check) {
        stop("Slot ", i, " must be an array of dimension ",
             paste(c(SOM@nsim, SOM@maxage, SOM@nyears), collapse = ", "))
      }
    }
  }

  # HistN check
  if (!length(slot(SOM, "HistN"))) {
    HistN <- array(0, c(SOM@nsim, SOM@maxage, SOM@nyears, 2))
    HistN[, 1, , ] <- 1000
    for (y in 2:SOM@nyears) {
      ZNOS <- SOM@Mjuv_NOS[, 2:SOM@maxage - 1, y-1] + SOM@HistFPT[, 2:SOM@maxage - 1, y-1, 1]
      HistN[, 2:SOM@maxage, y, 1] <- HistN[, 2:SOM@maxage - 1, y-1, 1] * exp(-ZNOS)

      ZHOS <- SOM@Mjuv_HOS[, 2:SOM@maxage - 1, y-1] + SOM@HistFPT[, 2:SOM@maxage - 1, y-1, 2]
      HistN[, 2:SOM@maxage, y, 2] <- HistN[, 2:SOM@maxage - 1, y-1, 2] * exp(-ZHOS)
    }
    slot(SOM, "HistN") <- HistN
  }
  dim_i <- dim(slot(SOM, "HistN"))
  if (is.null(dim_i)) {
    stop("Slot ", i, " must be an array of dimension ",
         paste(c(SOM@nsim, SOM@maxage, SOM@nyears, 2), collapse = ", "))
  } else {
    dim_check <- length(dim_i) == 4 && all(dim_i == c(SOM@nsim, SOM@maxage, SOM@nyears, 2))
    if (!dim_check) {
      stop("Slot ", i, " must be an array of dimension ",
           paste(c(SOM@nsim, SOM@maxage, SOM@nyears), collapse = ", "))
    }
  }
  if (is.null(dim_i)) {
    stop("Slot ", i, " must be an array of dimension ",
         paste(c(SOM@nsim, SOM@maxage, SOM@nyears, 2), collapse = ", "))
  } else {
    dim_check <- length(dim_i) == 4 && all(dim_i == c(SOM@nsim, SOM@maxage, SOM@nyears, 2))
    if (!dim_check) {
      stop("Slot ", i, " must be an array of dimension ",
           paste(c(SOM@nsim, SOM@maxage, SOM@nyears), collapse = ", "))
    }
  }

  return(SOM)
}
