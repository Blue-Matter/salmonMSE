
#' @name salmonMSE-int
#' @title Internal salmonMSE functions for converting operating model inputs and outputs
#'
#' @description
#' - [SOM2MOM()] converts a salmon operating model (\linkS4class{SOM}) to a multi-stock operating model (\linkS4class{MOM})
#' - [make_Stock()] creates the \linkS4class{Stock} object (openMSE) corresponding to salmon life stage
#' - [make_Fleet()] creates the \linkS4class{Fleet} object (openMSE) corresponding to the fishery that interacts with the various salmon life stages
#' - [multiHist2SHist()] converts the openMSE historical reconstruction into a salmon Hist object (\linkS4class{SHist})
#' - [MMSE2SMSE()] converts the openMSE projection output, along with additional state variables recorded in [salmonMSE_env], into a salmon MSE object (\linkS4class{SMSE})
#'
#' [salmonMSE()] is the wrapper function that coordinates the simulation and the output.
#'
#' @param SOM An object of class \linkS4class{SOM}
#' @param check Logical, whether to check the `SOM` object using [check_SOM()]
#' @export
#' @return
#' `SOM2MOM`: \linkS4class{MOM} object
SOM2MOM <- function(SOM, check = TRUE) {
  if (check) SOM <- check_SOM(SOM)

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
  Stocks <- make_Stock_objects(SOM)
  np <- length(Stocks)

  # Single fleet object ----
  nf <- 1
  Fleets <- make_Fleet_objects(SOM)

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

  if (do_hatchery || habitat_change || SOM@s_enroute < 1) {

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
    }

    fitness_args <- list()
    hatchery_args <- list(
      egg_local = egg_local
    )

    if (do_hatchery) {

      hatchery_args <- list(
        ptarget_NOB = SOM@ptarget_NOB,
        pmax_esc = SOM@pmax_esc,
        pmax_NOB = SOM@pmax_NOB,
        fec_brood = SOM@fec_brood,
        egg_local = egg_local,
        p_female = SOM@p_female,
        s_yearling = SOM@s_egg_smolt,
        s_subyearling = SOM@s_egg_subyearling,
        p_yearling = p_yearling,
        phatchery = SOM@phatchery,
        premove_HOS = SOM@premove_HOS,
        s_prespawn = SOM@s_prespawn,
        gamma = SOM@gamma,
        m = SOM@m
      )

      if (any(SOM@fitness_type == "Ford")) {
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
            #zbar_start = SOM@zbar_start # Now assigned by salmonMSE::salmonMSE()
            theta = SOM@theta
          )
        })
      }
    }

    # Natural smolt production from NOS and HOS escapement and habitat
    Rel[[1]] <- makeRel_smolt(
      p_smolt = 1, p_naturalsmolt = 1, p_natural = 3, p_hatchery = ifelse(do_hatchery, 6, NA_real_), output = "natural",
      s_enroute = SOM@s_enroute, p_female = SOM@p_female, fec = SOM@fec, SRRpars = SRRpars,
      hatchery_args = hatchery_args, fitness_args = fitness_args
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
      p_smolt = 4, p_naturalsmolt = 1, p_natural = 3, p_hatchery = 6, output = "hatchery",
      s_enroute = SOM@s_enroute, p_female = SOM@p_female, fec = SOM@fec, SRRpars = SRRpars,
      hatchery_args = hatchery_args, fitness_args = fitness_args
    )
  }

  MOM@Rel <- Rel
  MOM@cpars$control <- list(HistRel = FALSE, HermEq = FALSE)

  return(MOM)
}


#' Check inputs to SOM object
#'
#' Ensures that the slots in the [salmonMSE::SOM-class] object have the correct dimensions. Function will
#' update some slots to their full dimensions.
#'
#' @param SOM [salmonMSE::SOM-class] object
#' @returns Updated [salmonMSE::SOM-class] object with full dimensions in various slots
#' @export
check_SOM <- function(SOM) {

  ### Base slots ----
  if (!length(SOM@nsim)) stop("Need SOM@nsim")
  if (!length(SOM@nyears)) stop("Need SOM@nyears")
  if (!length(SOM@proyears)) stop("Need SOM@proyears")
  if (!length(SOM@seed)) stop("Need SOM@seed")
  nsim <- SOM@nsim

  ### Check Bio ----
  Bio <- SOM@Bio

  if (!length(Bio@maxage)) stop("Need SOM@Bio@maxage")
  maxage <- Bio@maxage
  years <- SOM@nyears + SOM@proyears

  Bio <- check_maxage2array(Bio, "p_mature", maxage, nsim, years)
  Bio <- check_numeric(Bio, "p_female", default = 0.5)
  Bio <- check_numeric(Bio, "s_enroute", default = 1)
  Bio <- check_maxage(Bio, "fec", maxage)

  # SRR pars
  slot(Bio, "SRrel") <- match.arg(slot(Bio, "SRrel"), choices = c("BH", "Ricker"))
  Bio <- check_numeric2nsim(Bio, "kappa", nsim)

  if (SOM@Bio@SRrel == "BH") {
    Bio <- check_numeric2nsim(Bio, "capacity_smolt", nsim)
  } else {
    Bio <- check_numeric2nsim(Bio, "Smax", nsim)
  }

  # Mjuv_NOS
  Bio <- check_maxage2array(Bio, "Mjuv_NOS", maxage, nsim, years)

  # phi
  if (!length(Bio@phi)) {
    Bio@phi <- local({
      surv_juv <- sapply(1:SOM@nsim, function(x) {
        calc_survival(Bio@Mjuv_NOS[x, , 1], Bio@p_mature[x, , 1])
      }) # maxage x nsim
      EPR <- t(surv_juv) * Bio@p_mature[, , 1]
      colSums(t(EPR) * Bio@p_female * Bio@fec)
    })
  }
  Bio <- check_numeric2nsim(Bio, "phi", nsim)

  # Habitat
  Habitat <- SOM@Habitat
  Habitat <- check_numeric(Habitat, "capacity_smolt_improve", default = 1)
  Habitat <- check_numeric(Habitat, "kappa_improve", default = 1)

  # Hatchery
  Hatchery <- SOM@Hatchery
  Hatchery <- check_numeric(Hatchery, "n_yearling", default = 0)
  Hatchery <- check_numeric(Hatchery, "n_subyearling", default = 0)

  do_hatchery <- Hatchery@n_yearling > 0 || Hatchery@n_subyearling > 0

  if (do_hatchery) {
    Hatchery <- check_numeric(Hatchery, "s_prespawn", default = 1)
    Hatchery <- check_numeric(Hatchery, "s_egg_smolt", default = 1)
    Hatchery <- check_numeric(Hatchery, "s_egg_subyearling", default = 1)

    if (!length(Hatchery@Mjuv_HOS)) Hatchery@Mjuv_HOS <- Bio@Mjuv_NOS
    Hatchery <- check_maxage2array(Hatchery, "Mjuv_HOS")

    Hatchery <- check_numeric(Hatchery, "gamma", default = 1)
    Hatchery <- check_numeric(Hatchery, "m", default = 0)
    Hatchery <- check_numeric(Hatchery, "pmax_esc", default = 0.75)
    Hatchery <- check_numeric(Hatchery, "pmax_NOB", default = 1)
    Hatchery <- check_numeric(Hatchery, "ptarget_NOB", default = 0.9)
    Hatchery <- check_numeric(Hatchery, "phatchery", default = 0.8)
    Hatchery <- check_numeric(Hatchery, "premove_HOS", default = 0)

    if (!length(Hatchery@fec_brood)) Hatchery@fec_brood <- Bio@fec
    Hatchery <- check_maxage(Hatchery, "fec_brood", maxage)

    Hatchery <- check_numeric(Hatchery, "fitness_type", size = 2)
    Hatchery <- check_numeric(Hatchery, "theta", size = 2)
    Hatchery <- check_numeric(Hatchery, "rel_loss", size = 3)

    if (!is.array(Hatchery@zbar_start)) {
      Hatchery <- check_numeric(Hatchery, "zbar_start", size = 2)
      Hatchery@zbar_start <- array(Hatchery@zbar_start, c(2, nsim, maxage)) %>%
        aperm(c(2, 3, 1))
    } else {
      Hatchery <- check_array(Hatchery, "zbar_start", c(nsim, maxage, 2))
    }

    Hatchery <- check_numeric(Hatchery, "fitness_variance")
    Hatchery <- check_numeric(Hatchery, "selection_strength")
    Hatchery <- check_numeric(Hatchery, "heritability")
    Hatchery <- check_numeric(Hatchery, "fitness_floor", default = 0.5)
  }

  # Harvest
  Harvest <- SOM@Harvest
  Harvest <- check_numeric(Harvest, "u_preterminal")
  Harvest <- check_numeric(Harvest, "u_terminal")
  Harvest <- check_numeric(Harvest, "MSF")
  if (Harvest@MSF) Harvest <- check_numeric(Harvest, "release_mortality", size = 2)
  if (Harvest@u_preterminal > 0) {
    Harvest <- check_numeric(Harvest, "vulPT", size = maxage)
  } else if (!length(Harvest@vulPT)) {
    Harvest@vulPT <- rep(0, maxage)
  }
  if (Harvest@u_terminal > 0) {
    Harvest <- check_numeric(Harvest, "vulT", size = maxage)
  } else if (!length(Harvest@vulT)) {
    Harvest@vulT <- rep(0, maxage)
  }

  # Historical
  Historical <- SOM@Historical

  var_F <- c("HistFPT", "HistFT")
  for (i in var_F) {
    if (!length(slot(Historical, i))) slot(Historical, i) <- rep(0, SOM@nyears)
    if (!is.array(slot(Historical, i))) {
      Historical <- check_numeric(Historical, i, SOM@nyears)
      slot(Historical, i) <- slot(Historical, i) %>% array(c(SOM@nyears, nsim, 2)) %>% aperm(c(2, 1, 3))
    }
    Historical <- check_array(Historical, i, dims = c(nsim, SOM@nyears, 2))
  }

  if (length(Historical@HistSpawner)) {
    Historical <- check_array(Historical, "HistSpawner", dims = c(nsim, maxage, SOM@nyears, 2))
  }

  if (!length(Historical@HistN)) {

    HistN <- array(0, c(nsim, maxage, SOM@nyears + 1, 2))
    HistN[, 1, , ] <- 1000
    for (y in seq(2, SOM@nyears + 1)) {
      ZNOS <- Bio@Mjuv_NOS[, 2:maxage - 1, y-1] + Historical@HistFPT[, 2:maxage - 1, y-1, 1]
      HistN[, 2:maxage, y, 1] <- HistN[, 2:maxage - 1, y-1, 1] * exp(-ZNOS)

      ZHOS <- Hatchery@Mjuv_HOS[, 2:maxage - 1, y-1] + Historical@HistFPT[, 2:maxage - 1, y-1, 2]
      HistN[, 2:maxage, y, 2] <- HistN[, 2:maxage - 1, y-1, 2] * exp(-ZHOS)
    }
    Historical@HistN <- HistN

  }
  Historical <- check_array(Historical, "HistN", dims = c(nsim, maxage, SOM@nyears+1, 2))

  SOM@Bio <- Bio
  SOM@Habitat <- Habitat
  SOM@Hatchery <- Hatchery
  SOM@Harvest <- Harvest
  SOM@Historical <- Historical

  return(SOM)
}

check_numeric <- function(object, name, size = 1, default) {
  object_name <- as.character(substitute(object))

  if (!length(slot(object, name))) {
    if (!missing(default)) {
      slot(object, name) <- default
    } else {
      stop(paste0("Need ", object_name, "@", name, " (numeric)"))
    }
  } else if (length(slot(object, name)) != size) {
    stop(paste0("Length of ", object_name, "@", name, " needs to be ", size))
  }

  return(invisible(object))
}

check_maxage <- function(object, name, maxage) {
  object_name <- as.character(substitute(object))

  if (!length(slot(object, name)) || length(slot(object, name)) != maxage) {
    if (!missing(default)) {
      slot(object, name) <- default
    } else {
      stop(paste0("Need ", object_name, "@", name, " (maxage vector)"))
    }
  }
  return(invisible(object))
}

check_maxage2array <- function(object, name, maxage, nsim, years) {
  object_name <- as.character(substitute(object))

  if (!is.array(slot(object, name))) {
    if (length(slot(object, name) != maxage)) {
      stop(paste0(object_name, "@", name, " needs to be length maxage (", maxage, ")"))
    }
    slot(object, name) <- array(slot(object, name), c(maxage, nsim, years)) %>%
      aperm(c(2, 1, 3))
  }

  dim_i <- dim(slot(object, name))
  dim_check <- length(dim_i) == 3 && all(dim(slot(object, name) == c(nsim, maxage, years)))
  if (!dim_check) {
    stop(
      paste0(object_name, "@", name, " must be an array with dimension ",
             paste(c(nsim, maxage, years), collapse = ", "))
    )
  }
  return(invisible(object))
}

check_array <- function(object, name, dims) {
  object_name <- as.character(substitute(object))

  dim_i <- dim(slot(object, name))

  dim_check <- length(dim_i) == length(dims) && all(dim(slot(object, name) == dims))
  if (!dim_check) {
    stop(
      paste0(object_name, "@", name, " must be an array with dimension ",
             paste(dims, collapse = ", "))
    )
  }
  return(invisible(object))
}

check_numeric2nsim <- function(object, name, nsim) {
  object_name <- as.character(substitute(object))

  if (!length(slot(object, name))) {
    stop(paste0("Need ", object_name, "@", name))
  }

  if (length(slot(object, name)) == 1) slot(object, name) <- rep(slot(object, name), nsim)

  dim_check <- all(dim(slot(object, name) == c(nsim, maxage, years)))
  if (length(slot(object, name)) != nsim) {
    stop(paste0(object_name, "@", name, " needs to length nsim (", nsim, ")"))
  }
  return(invisible(object))
}
