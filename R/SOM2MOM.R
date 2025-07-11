
#' @name salmonMSE-int
#' @title Internal salmonMSE functions for converting operating model inputs and outputs
#'
#' @description
#' - [SOM2MOM()] converts a salmon operating model (\linkS4class{SOM}) to a multi-stock operating model ([MSEtool::MOM-class])
#' - [make_Stock()] creates the [MSEtool::Stock-class] object (openMSE) corresponding to salmon life stage
#' - [make_Fleet()] creates the [MSEtool::Fleet-class] object (openMSE) corresponding to the fishery that interacts with the various salmon life stages
#' - [multiHist2SHist()] converts the openMSE historical reconstruction into a salmon Hist object (\linkS4class{SHist})
#' - [MMSE2SMSE()] converts the openMSE projection output, along with additional state variables recorded in [salmonMSE_env], into a salmon MSE object (\linkS4class{SMSE})
#' - [make_Harvest_MMP()] creates a multi-stock management procedure for the harvest component of the operating model by
#' specifying exploitation rates through updating the formal arguments for [Harvest_MMP()]
#'
#' [salmonMSE()] is the wrapper function that coordinates the simulation and the output.
#'
#' @param SOM An object of class \linkS4class{SOM}
#' @param check Logical, whether to check the `SOM` object using [check_SOM()]
#' @export
#' @return
#' `SOM2MOM`: [MSEtool::MOM-class] object
SOM2MOM <- function(SOM, check = TRUE) {
  if (check) SOM <- check_SOM(SOM)

  MOM <- suppressMessages(new("MOM"))
  MOM@Name <- SOM@Name
  MOM@proyears <- 2 * SOM@proyears
  MOM@nsim <- SOM@nsim
  MOM@interval <- 1
  MOM@pstar <- 0.5
  MOM@maxF <- 20
  MOM@reps <- 1

  # Stock objects
  ns <- length(SOM@Bio)
  Stocks_s <- lapply(1:ns, function(s) make_Stock_objects(SOM, s = s))
  Stocks <- do.call(c, Stocks_s)
  #ns <- length(Stocks_s)
  np_s <- sapply(Stocks_s, length)
  np <- length(Stocks)
  n_g <- sapply(1:ns, function(s) SOM@Bio[[s]]@n_g)
  n_r <- sapply(1:ns, function(s) SOM@Hatchery[[s]]@n_r)

  do_hatchery <- vapply(SOM@Hatchery, function(Hatchery) sum(Hatchery@n_yearling, Hatchery@n_subyearling) > 0, logical(1))
  has_strays <- vapply(1:ns, function(s) any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external), logical(1))

  # Fleet objects (one per openMSE population/life stage) ----
  nf <- 1
  Fleets_s <- lapply(1:ns, function(s) make_Fleet_objects(SOM, s = s))
  Fleets <- do.call(c, Fleets_s)

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
  # Ignore Complexes for now

  # Specify smolt natural production from NOS population
  pindex <- .make_stock_index(ns, n_g, n_r, np_s)
  SSBfrom <- matrix(0, np, np)
  for (s in 1:ns) {
    for (g in 1:n_g[s]) {
      # Juv NOS predicted from NOS escapement of the corresponding LHG. Note p_LHG must be full rank
      p_nat_juv <- pindex$p[pindex$s == s & pindex$g == g & pindex$origin == "natural" & pindex$stage == "juvenile"]
      p_nat_esc <- pindex$p[pindex$s == s & pindex$g == g & pindex$origin == "natural" & pindex$stage == "escapement"]
      SSBfrom[p_nat_juv, p_nat_esc] <- 1

      # Integrated hatchery: triggers SRRfun so that there is hatchery production only when there is natural escapement for all LHG
      if (do_hatchery[s]) {
        for (r in 1:n_r[s]) {
          p_hat_juv <- pindex$p[pindex$s == s & pindex$r == r & pindex$origin == "hatchery" & pindex$stage == "juvenile"]
          SSBfrom[p_hat_juv, p_nat_esc] <- 1
        }
      }
    }
  }

  # Move stocks around
  # First generate the escapement, then move juveniles to adult
  maxage_s <- vapply(SOM@Bio, slot, numeric(1), "maxage")
  nage <- 2 * maxage_s[1] + 1
  nyears <- 2 * SOM@nyears
  proyears <- 2 * SOM@proyears

  Herm_mature <- lapply(1:ns, function(s) {

    # Maturation occurs for even age classes (age class 2, 4, 6) at the beginning of even time steps.
    prop_mature_NOS <- prop_escapement_NOS <- array(0, c(SOM@nsim, nage, nyears + proyears))
    prop_mature_NOS[, 2 * seq(1, maxage_s[s]), seq(2, nyears + proyears, 2)] <- SOM@Bio[[s]]@p_mature

    Herm_mature_NOS <- sapply(1:(nyears + proyears), function(y) {
      sapply(1:SOM@nsim, function(x) solve_Herm(prop_mature_NOS[x, , y]))
    }, simplify = "array") %>%
      aperm(c(2, 1, 3))

    # In reality, escapement occurs for even age classes (age class 2, 4, 6) at the end of even time steps.
    # In openMSE, we do escapement for odd age classes at the beginning of the subsequent odd time steps
    first_mature_age <- sapply(seq(1, SOM@nyears + SOM@proyears), function(y) {
      sapply(1:SOM@nsim, function(x) which(SOM@Bio[[s]]@p_mature[x, , y] > 0)[1])
    })

    prop_escapement_NOS[, 2 * seq(1, maxage_s[s]) + 1, seq(1, nyears + proyears, 2)] <- sapply(1:SOM@nsim, function(x) {
      sapply(1:(SOM@nyears + SOM@proyears), function(y) ifelse(1:maxage_s[s] >= first_mature_age[x, y], 1, 0))
    }, simplify = "array") %>%
      aperm(c(3, 1, 2))

    Herm_escapement_NOS <- sapply(1:(nyears + proyears), function(y) {
      sapply(1:SOM@nsim, function(x) solve_Herm(prop_escapement_NOS[x, , y]))
    }, simplify = "array") %>%
      aperm(c(2, 1, 3))

    # For NOS
    Herm <- list()
    Herm_names <- character(0)
    for (g in 1:n_g[s]) {
      p_nat_esc <- pindex$p[pindex$s == s & pindex$g == g & pindex$origin == "natural" & pindex$stage == "escapement"]
      p_nat_rec <- pindex$p[pindex$s == s & pindex$g == g & pindex$origin == "natural" & pindex$stage == "recruitment"]
      p_nat_juv <- pindex$p[pindex$s == s & pindex$g == g & pindex$origin == "natural" & pindex$stage == "juvenile"]

      Herm <- c(Herm, list(Herm_escapement_NOS, Herm_mature_NOS))
      Herm_names <- c(Herm_names, paste0("H_", p_nat_esc, "_", p_nat_rec), paste0("H_", p_nat_rec, "_", p_nat_juv))
    }

    # HOS
    if (do_hatchery[s] || has_strays[s]) {
      for (r in 1:n_r[s]) {

        # Maturation occurs for even age classes (age class 2, 4, 6) at the beginning of even time steps.
        prop_mature_HOS <- prop_escapement_HOS <- array(0, c(SOM@nsim, nage, nyears + proyears))
        prop_mature_HOS[, 2 * seq(1, maxage_s[s]), seq(2, nyears + proyears, 2)] <- SOM@Hatchery[[s]]@p_mature_HOS[, , , r]

        Herm_mature_HOS <- sapply(1:(nyears + proyears), function(y) {
          sapply(1:SOM@nsim, function(x) solve_Herm(prop_mature_HOS[x, , y]))
        }, simplify = "array") %>%
          aperm(c(2, 1, 3))

        # In reality, escapement occurs for even age classes (age class 2, 4, 6) at the end of even time steps.
        # In openMSE, we do escapement for odd age classes at the beginning of the subsequent odd time steps
        first_mature_age <- sapply(seq(1, SOM@nyears + SOM@proyears), function(y) {
          sapply(1:SOM@nsim, function(x) which(SOM@Hatchery[[s]]@p_mature_HOS[x, , y, r] > 0)[1])
        })

        prop_escapement_HOS[, 2 * seq(1, maxage_s[s]) + 1, seq(1, nyears + proyears, 2)] <- sapply(1:SOM@nsim, function(x) {
          sapply(1:(SOM@nyears + SOM@proyears), function(y) ifelse(1:maxage_s[s] >= first_mature_age[x, y], 1, 0))
        }, simplify = "array") %>%
          aperm(c(3, 1, 2))

        Herm_escapement_HOS <- sapply(1:(nyears + proyears), function(y) {
          sapply(1:SOM@nsim, function(x) solve_Herm(prop_escapement_HOS[x, , y]))
        }, simplify = "array") %>%
          aperm(c(2, 1, 3))

        p_hat_esc <- pindex$p[pindex$s == s & pindex$r == r & pindex$origin == "hatchery" & pindex$stage == "escapement"]
        p_hat_rec <- pindex$p[pindex$s == s & pindex$r == r & pindex$origin == "hatchery" & pindex$stage == "recruitment"]
        p_hat_juv <- pindex$p[pindex$s == s & pindex$r == r & pindex$origin == "hatchery" & pindex$stage == "juvenile"]

        Herm <- c(Herm, list(Herm_escapement_HOS, Herm_mature_HOS))
        Herm_names <- c(Herm_names, paste0("H_", p_hat_esc, "_", p_hat_rec), paste0("H_", p_hat_rec, "_", p_hat_juv))
      }
    }

    structure(Herm, names = Herm_names)
  })
  Herm_mature <- do.call(c, Herm_mature)

  MOM@SexPars <- list(SSBfrom = SSBfrom, Herm = Herm_mature, share_par = FALSE)

  # Rel
  # The predicted smolt production (juvenile NOS) is from the adult escapement (natural and hatchery)
  # Combines hatchery and habitat dynamics
  # Hatchery: Presence/absence of HOS (from strays or own system)
  # Habitat: Update Perr_y of juvenile NOS from adult escapement from stage specific relationships
  Rel_s <- lapply(1:ns, function(s) {
    Bio <- SOM@Bio[[s]]
    Habitat <- SOM@Habitat[[s]]
    Hatchery <- SOM@Hatchery[[s]]
    S <- Stocks_s[[s]]
    SRRpars <- S[[1]]$cpars_bio$SRR$SRRpars # Assume SRR pars are identical for all LHG

    p_LHG <- Bio@p_LHG

    pind <- pindex[pindex$s == s, ]

    # Not needed if all of the following are true: no hatchery, no strays, s_enroute = 1, n_g = 1
    # SRR pars in the OM should suffice in that case (the escapement is the spawning output)
    do_hatchery_s <- do_hatchery[s]
    has_strays_s <- has_strays[s]
    g_s <- n_g[s]
    r_s <- n_r[s]

    use_smolt_func <- TRUE

    Rel_smolt_g <- list()
    Rel_hatchrel <- list()
    Rel_SAR_NOS_g <- list()
    Rel_SAR_HOS <- list()

    if (use_smolt_func) {

      if (do_hatchery_s) {
        # Determine the total number of eggs needed from the number of yearling and subyearling releases, their survival from egg stage
        # and fecundity of broodtake (identical between natural and hatchery escapement)
        # This is a management action, cannot be stochastic
        # No imports
        egg_yearling <- ifelse(sum(Hatchery@n_yearling) > 0, sum(Hatchery@n_yearling)/Hatchery@s_egg_smolt, 0)
        egg_subyearling <- ifelse(sum(Hatchery@n_subyearling) > 0, sum(Hatchery@n_subyearling)/Hatchery@s_egg_subyearling, 0)
        egg_target <- egg_yearling + egg_subyearling

        p_yearling <- Hatchery@n_yearling/sum(Hatchery@n_yearling, Hatchery@n_subyearling) # Vector by release strategy
        p_subyearling <- Hatchery@n_subyearling/sum(Hatchery@n_yearling, Hatchery@n_subyearling) # Vector by release strategy
      } else {
        egg_target <- p_yearling <- p_subyearling <- 0
      }

      fitness_args <- list()
      hatchery_args <- list(
        egg_target = egg_target
      )
      habitat_args <- list(
        use_habitat = Habitat@use_habitat
      )
      stray_args <- list(
        stray_external = Hatchery@stray_external
      )

      donate_stray <- SOM@stray[s, s] < 1
      receive_stray <- any(SOM@stray[-s, s] > 0)

      if (donate_stray || receive_stray) {
        if (any(n_r > 1)) stop("Straying matrix not yet configured for multiple release strategies.")
        stray_args$stray_matrix <- SOM@stray
        stray_args$p_donor <- pindex[pindex$stage == "escapement" & pindex$origin == "hatchery", c("s", "p")]
      }

      if (habitat_args$use_habitat) {
        habitat_args$Habitat <- Habitat
        habitat_args$nyears_om <- 2 * SOM@nyears
      }

      if (do_hatchery_s || has_strays_s) {

        hatchery_args <- list(
          ptarget_NOB = Hatchery@ptarget_NOB,
          pmax_esc = Hatchery@pmax_esc,
          pmax_NOB = Hatchery@pmax_NOB,
          fec_brood = Hatchery@fec_brood,
          egg_target = egg_target,
          brood_import = Hatchery@brood_import,
          p_female = Bio@p_female,
          s_yearling = Hatchery@s_egg_smolt,
          s_subyearling = Hatchery@s_egg_subyearling,
          p_yearling = p_yearling,
          p_subyearling = p_subyearling,
          phatchery = Hatchery@phatchery,
          premove_HOS = Hatchery@premove_HOS,
          s_prespawn = Hatchery@s_prespawn,
          gamma = Hatchery@gamma,
          m = Hatchery@m
        )

        if (any(Hatchery@fitness_type == "Ford") && (do_hatchery_s || has_strays_s)) {
          fitness_args <- list(
            fitness_type = Hatchery@fitness_type,
            rel_loss = Hatchery@rel_loss,
            phenotype_variance = Hatchery@phenotype_variance,
            fitness_variance = Hatchery@fitness_variance,
            fitness_floor = Hatchery@fitness_floor,
            heritability = Hatchery@heritability,
            #zbar_start = Hatchery@zbar_start # Now assigned by salmonMSE::salmonMSE()
            theta = Hatchery@theta
          )
        }
      }

      p_nat_esc <- pind$p[pind$origin == "natural" & pind$stage == "escapement"] # length ng
      p_hat_esc <- if(do_hatchery_s || has_strays_s) {
        pind$p[pind$origin == "hatchery" & pind$stage == "escapement"]
      } else {
        NULL
      }

      for (g in seq_len(g_s)) {
        # Natural smolt production from NOS and HOS escapement (including strays), habitat, and/or en-route mortality
        p_nat_smolt <- pind$p[pind$origin == "natural" & pind$g == g & pind$stage == "juvenile"]

        Rel_smolt_g[[g]] <- makeRel_smolt(
          p_smolt = p_nat_smolt, s = s, p_natural = p_nat_esc, p_hatchery = p_hat_esc,
          p_stray = setdiff(stray_args$p_donor[, "p"], p_hat_esc),
          output = "natural", s_enroute = Bio@s_enroute, p_female = Bio@p_female, fec = Bio@fec,
          SRRpars = SRRpars,
          hatchery_args = hatchery_args, fitness_args = fitness_args, habitat_args = habitat_args,
          stray_args = stray_args,
          g = g, prop_LHG = p_LHG[g]
        )

        # Marine survival of natural origin fish
        if (do_hatchery_s && Hatchery@fitness_type[1] != "none" && Hatchery@rel_loss[3] > 0) {
          Rel_SAR_NOS_g[[g]] <- makeRel_SAR(
            p_smolt = p_nat_smolt, s = s, envir = "natural",
            rel_loss = Hatchery@rel_loss[3], nyears = 2 * SOM@nyears,
            Mbase = S[[p_nat_smolt]]$cpars_bio$M_ageArray[, , 2 * SOM@nyears + seq(1, 2 * SOM@proyears)]
          )
        }
      }

      if (do_hatchery_s) {
        for (r in seq_len(r_s)) {
          p_hat_smolt <- pind$p[pind$origin == "hatchery" & pind$r == r & pind$stage == "juvenile"]

          # Hatchery smolt releases from NOS and HOS escapement
          Rel_hatchrel[[r]] <- makeRel_smolt(
            p_smolt = p_hat_smolt, s = s, p_natural = p_nat_esc, p_hatchery = p_hat_esc,
            p_stray = setdiff(stray_args$p_donor[, "p"], p_hat_esc),
            output = "hatchery", s_enroute = Bio@s_enroute, p_female = Bio@p_female, fec = Bio@fec,
            SRRpars = SRRpars,
            hatchery_args = hatchery_args, fitness_args = fitness_args, habitat_args = habitat_args,
            stray_args = stray_args,
            r = r
          )

          # Marine survival of hatchery origin fish
          if (Hatchery@fitness_type[2] != "none" && Hatchery@rel_loss[3] > 0) {
            Rel_SAR_HOS[[1]] <- makeRel_SAR(
              p_smolt = p_hat_smolt, s = s, envir = "hatchery",
              rel_loss = Hatchery@rel_loss[3], nyears = 2 * SOM@nyears,
              Mbase = S[[p_hat_smolt]]$cpars_bio$M_ageArray[, , 2 * SOM@nyears + seq(1, 2 * SOM@proyears)]
            )
          }
        }

      }
    }
    Rel <- c(Rel_smolt_g, Rel_SAR_NOS_g, Rel_SAR_HOS, Rel_hatchrel)

    return(Rel)
  })

  MOM@Rel <- do.call(c, Rel_s)
  MOM@cpars$control <- list(HistRel = FALSE, HermEq = FALSE)

  return(MOM)
}

NAor0 <- function(x) !length(x) || is.na(x) || !x

#' Check inputs to SOM object
#'
#' Ensures that the slots in the [salmonMSE::SOM-class] object have the correct dimensions. Function will
#' update some slots to their full dimensions.
#'
#' @param SOM [salmonMSE::SOM-class] object
#' @param silent Logical, whether to report progress in console
#' @returns Updated [salmonMSE::SOM-class] object with full dimensions in various slots
#' @export
check_SOM <- function(SOM, silent = FALSE) {

  ### Base slots ----
  if (!length(SOM@nsim)) stop("Need SOM@nsim")
  if (!length(SOM@nyears)) stop("Need SOM@nyears")
  if (!length(SOM@proyears)) stop("Need SOM@proyears")
  if (!length(SOM@seed)) stop("Need SOM@seed")
  nsim <- SOM@nsim
  years <- SOM@nyears + SOM@proyears

  # Convert sub class objects to lists
  obj <- c("Bio", "Habitat", "Hatchery", "Harvest", "Historical")
  for (i in obj) {
    if (inherits(slot(SOM, i), i)) slot(SOM, i) <- list(slot(SOM, i))
  }

  # Loop over populations
  ns <- length(SOM@Bio)

  # Check straying
  if (!length(SOM@stray)) SOM@stray <- diag(ns)

  for (s in 1:ns) {
    if (!silent && ns > 1) message("Checking parameters for population ", s)

    ### Check Bio ----
    Bio <- SOM@Bio[[s]]

    Bio <- check_numeric(Bio, "maxage")
    maxage <- Bio@maxage

    Bio <- check_numeric(Bio, "n_g", default = 1)
    Bio <- check_numeric(Bio, "p_LHG", size = Bio@n_g, default = rep(1/Bio@n_g, Bio@n_g))

    Bio <- check_maxage2array(Bio, "p_mature", maxage, nsim, years)
    Bio <- check_numeric(Bio, "p_female", default = 0.5)
    Bio <- check_numeric(Bio, "s_enroute", default = 1)
    Bio <- check_maxage(Bio, "fec", maxage)

    ### Check Habitat (initial) ----
    Habitat <- SOM@Habitat[[s]]
    Habitat <- check_numeric(Habitat, "use_habitat", default = FALSE)

    if (length(dim(Bio@Mjuv_NOS)) == 3 && Bio@n_g == 1) {
      Bio <- check_maxage2array(Bio, "Mjuv_NOS", maxage, nsim, years)
      Bio@Mjuv_NOS <- array(Bio@Mjuv_NOS, c(dim(Bio@Mjuv_NOS), Bio@n_g))
    }
    Bio <- check_maxage2garray(Bio, "Mjuv_NOS", maxage, nsim, years, Bio@n_g)

    if (!Habitat@use_habitat) {
      # SRR pars
      slot(Bio, "SRrel") <- match.arg(slot(Bio, "SRrel"), choices = c("BH", "Ricker"))
      Bio <- check_numeric2nsim(Bio, "kappa", nsim)

      if (Bio@SRrel == "BH") {
        Bio <- check_numeric2nsim(Bio, "capacity", nsim)
      } else {
        Bio <- check_numeric2nsim(Bio, "Smax", nsim)
      }

      if (!length(Bio@phi)) {
        Bio@phi <- sapply(1:SOM@nsim, function(x) {
          calc_phi(
            Mjuv = matrix(Bio@Mjuv_NOS[x, , 1, ], Bio@maxage, Bio@n_g),
            p_mature = matrix(Bio@p_mature[x, , 1], Bio@maxage, Bio@n_g),
            p_female = Bio@p_female,
            fec = matrix(Bio@fec, Bio@maxage, Bio@n_g),
            s_enroute = Bio@s_enroute,
            n_g = Bio@n_g,
            p_LHG = Bio@p_LHG
          )
        })
      }
      Bio <- check_numeric2nsim(Bio, "phi", nsim)
    } else {
      Habitat <- check_numeric(Habitat, "egg_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "egg_prod", default = 1)
      Habitat <- check_numeric(Habitat, "egg_capacity", default = Inf)

      Habitat <- check_numeric(Habitat, "fry_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "fry_prod", default = 0.4)
      Habitat <- check_numeric(Habitat, "fry_capacity", default = Inf)
      if (!length(Habitat@fry_sdev)) {
        Habitat@fry_sdev <- matrix(1, nsim, SOM@proyears)
      }
      Habitat <- check_maxage2matrix(Habitat, "fry_sdev", SOM@proyears, nsim)

      Habitat <- check_numeric(Habitat, "smolt_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "smolt_prod", default = 1)
      Habitat <- check_numeric(Habitat, "smolt_capacity", default = Inf)
      if (!length(Habitat@smolt_sdev)) {
        Habitat@smolt_sdev <- matrix(1, nsim, SOM@proyears)
      }
      Habitat <- check_maxage2matrix(Habitat, "smolt_sdev", SOM@proyears, nsim)
    }

    ### Check Hatchery ----
    Hatchery <- SOM@Hatchery[[s]]

    Hatchery <- check_numeric(Hatchery, "n_r", default = 1)
    Hatchery <- check_numeric(Hatchery, "n_yearling", size = Hatchery@n_r, default = rep(0, Hatchery@n_r))
    Hatchery <- check_numeric(Hatchery, "n_subyearling", size = Hatchery@n_r, default = rep(0, Hatchery@n_r))

    do_hatchery <- sum(Hatchery@n_yearling, Hatchery@n_subyearling) > 0

    if (!length(Hatchery@stray_external)) Hatchery@stray_external <- matrix(0, maxage, Hatchery@n_r)
    Hatchery <- check_maxage2matrix(Hatchery, "stray_external", Hatchery@n_r, maxage)

    has_strays <- any(SOM@stray[-s, s] > 0) || sum(Hatchery@stray_external)

    if (do_hatchery || has_strays) {
      Hatchery <- check_numeric(Hatchery, "s_prespawn", default = 1)
      Hatchery <- check_numeric(Hatchery, "s_egg_smolt", default = 1)
      Hatchery <- check_numeric(Hatchery, "s_egg_subyearling", default = 1)
      Hatchery <- check_numeric(Hatchery, "brood_import", size = maxage, default = rep(0, maxage))

      if (length(dim(Hatchery@Mjuv_HOS)) == 3 && Hatchery@n_r == 1) {
        Hatchery <- check_maxage2array(Hatchery, "Mjuv_HOS", maxage, nsim, years)
        Hatchery@Mjuv_HOS <- array(Hatchery@Mjuv_HOS, c(dim(Hatchery@Mjuv_HOS), 1))
      }
      Hatchery <- check_maxage2garray(Hatchery, "Mjuv_HOS", maxage, nsim, years, Hatchery@n_r)

      if (!length(Hatchery@p_mature_HOS)) {
        Hatchery@p_mature_HOS <- array(Bio@p_mature, c(nsim, maxage, years, Hatchery@n_r))
      }
      if (length(dim(Hatchery@p_mature_HOS)) == 3 && Hatchery@n_r == 1) {
        Hatchery <- check_maxage2array(Hatchery, "p_mature_HOS", maxage, nsim, years)
        Hatchery@p_mature_HOS <- array(Hatchery@p_mature_HOS, c(dim(Hatchery@Mjuv_HOS), 1))
      }
      Hatchery <- check_maxage2garray(Hatchery, "p_mature_HOS", maxage, nsim, years, Hatchery@n_r)

      Hatchery <- check_numeric(Hatchery, "m", default = 0)

      Hatchery <- check_numeric(Hatchery, "pmax_esc", default = 0.75)
      Hatchery <- check_numeric(Hatchery, "pmax_NOB", default = 1)
      Hatchery <- check_numeric(Hatchery, "ptarget_NOB", default = 0.9)

      if (!length(Hatchery@fec_brood)) Hatchery@fec_brood <- Bio@fec
      Hatchery <- check_maxage(Hatchery, "fec_brood", maxage)

      Hatchery <- check_numeric(Hatchery, "gamma", default = 1)

      Hatchery <- check_numeric(Hatchery, "phatchery", default = NA)
      Hatchery <- check_numeric(Hatchery, "premove_HOS", default = 0)

      Hatchery <- check_numeric(Hatchery, "fitness_type", size = 2)
      if (any(Hatchery@fitness_type == "Ford")) {
        Hatchery <- check_numeric(Hatchery, "theta", size = 2)
        Hatchery <- check_numeric(Hatchery, "rel_loss", size = 3)

        if (!is.array(Hatchery@zbar_start)) {
          Hatchery <- check_numeric(Hatchery, "zbar_start", size = 2)
          Hatchery@zbar_start <- array(Hatchery@zbar_start, c(2, nsim, maxage)) %>%
            aperm(c(2, 3, 1))
        } else {
          Hatchery <- check_array(Hatchery, "zbar_start", c(nsim, maxage, 2))
        }

        Hatchery <- check_numeric(Hatchery, "fitness_variance", default = 100)
        Hatchery <- check_numeric(Hatchery, "phenotype_variance", default = 10)

        if (length(Hatchery@heritability) == 1) {
          Hatchery@heritability <- rep(Hatchery@heritability, nsim)
        }
        Hatchery <- check_numeric(Hatchery, "heritability", size = nsim, default = rep(0.5, nsim))
        Hatchery <- check_numeric(Hatchery, "fitness_floor", default = 0.5)
      }
    } else {
      Hatchery <- check_numeric(Hatchery, "m", default = 0)
    }

    ### Check for PNI_wild calculation ----
    if (!do_hatchery && sum(Hatchery@stray_external)) {
      if (!length(Hatchery@heritability)) {
        stop("Need hatchery heritability parameter for PNI calculation (no hatchery, presence of external strays)")
      }
      if (!length(Hatchery@phenotype_variance)) {
        stop("Need hatchery phenotype_variance parameter for PNI calculation (no hatchery, presence of external strays)")
      }
      if (!length(Hatchery@fitness_variance)) {
        stop("Need hatchery fitness_variance parameter for PNI calculation (no hatchery, presence of external strays)")
      }
    }

    ### Check Harvest ----
    Harvest <- SOM@Harvest[[s]]

    Harvest <- check_numeric(Harvest, "type_PT", default = "u")
    Harvest <- check_numeric(Harvest, "type_T", default = "u")

    Harvest@type_PT <- match.arg(Harvest@type_PT, choices = c("u", "catch"))
    Harvest@type_T <- match.arg(Harvest@type_T, choices = c("u", "catch"))

    if (Harvest@type_PT == "u") {
      Harvest <- check_numeric(Harvest, "u_preterminal", default = 0)
      Harvest <- check_numeric(Harvest, "K_PT", default = NA_real_)
    } else {
      Harvest <- check_numeric(Harvest, "u_preterminal", default = NA_real_)
      Harvest <- check_numeric(Harvest, "K_PT", default = 0)
    }

    if (Harvest@type_T == "u") {
      Harvest <- check_numeric(Harvest, "u_terminal", default = 0)
      Harvest <- check_numeric(Harvest, "K_T", default = NA_real_)
    } else {
      Harvest <- check_numeric(Harvest, "u_terminal", default = NA_real_)
      Harvest <- check_numeric(Harvest, "K_T", default = 0)
    }

    Harvest <- check_numeric(Harvest, "MSF_PT", default = FALSE)
    Harvest <- check_numeric(Harvest, "MSF_T", default = FALSE)

    Harvest <- check_numeric(Harvest, "release_mort", size = 2, default = c(0, 0))

    if (!NAor0(Harvest@u_preterminal > 0) || !NAor0(Harvest@K_PT)) {
      Harvest <- check_maxage2matrix(Harvest, "vulPT", maxage, nsim)
    } else {
      Harvest@vulPT <- matrix(0, nsim, maxage)
    }
    if (!NAor0(Harvest@u_terminal > 0) || !NAor0(Harvest@K_T)) {
      Harvest <- check_maxage2matrix(Harvest, "vulT", maxage, nsim)
    } else {
      Harvest@vulT <- matrix(0, nsim, maxage)
    }

    ### Check Historical (initial) ----
    Historical <- SOM@Historical[[s]]

    var_F <- c("HistFPT", "HistFT")
    for (i in var_F) {
      if (!length(slot(Historical, i))) slot(Historical, i) <- rep(0, SOM@nyears)
      if (!is.array(slot(Historical, i))) {
        Historical <- check_numeric(Historical, i, SOM@nyears)
        slot(Historical, i) <- slot(Historical, i) %>% array(c(SOM@nyears, nsim, 2)) %>% aperm(c(2, 1, 3))
      }
      Historical <- check_array(Historical, i, dims = c(nsim, SOM@nyears, 2))
    }

    init_NOS <- init_HOS <- NA_real_

    if (!length(Historical@HistSpawner_NOS)) {
      init_NOS <- 1000
    } else if (!is.array(Historical@HistSpawner_NOS) && length(Historical@HistSpawner_NOS) == 1) {
      init_NOS <- Historical@HistSpawner_NOS
    } else {
      if (length(dim(Historical@HistSpawner_NOS)) == 3 && Bio@n_g == 1) {
        Historical <- check_array(Historical, "HistSpawner_NOS", dims = c(nsim, maxage, SOM@nyears))
        Historical@HistSpawner_NOS <- array(Historical@HistSpawner_NOS, c(dim(Historical@HistSpawner_NOS), 1))
      } else {
        Historical <- check_array(Historical, "HistSpawner_NOS", dims = c(nsim, maxage, SOM@nyears, Bio@n_g))
      }
    }

    if (!length(Historical@HistNjuv_NOS)) {
      if (SOM@nyears > maxage) stop("No historical abundance was provided for natural origin fish. Set historical years (SOM@nyears) < maxage")

      HistNjuv_NOS <- array(0, c(nsim, maxage, SOM@nyears + 1, Bio@n_g))
      HistNjuv_NOS[, maxage, SOM@nyears + 1, ] <- local({
        init_Njuv <- init_NOS/(1 - Harvest@u_terminal)/(1 - Harvest@u_preterminal)
        outer(init_Njuv, Bio@p_LHG) # matrix nsim x n_g
      })

      for (a in seq(maxage - 1, 1)) {
        y <- SOM@nyears + a - maxage + 1
        if (y > 0) {
          ZNOS <- matrix(Harvest@vulPT[, a] * Historical@HistFPT[, y, 1], nsim, Bio@n_g) + Bio@Mjuv_NOS[, a, y, ]
          HistNjuv_NOS[, a, y, ] <- HistNjuv_NOS[, a+1, y+1, ] * exp(ZNOS) / (1 - matrix(Bio@p_mature[, a, y], nsim, Bio@n_g))
        }
      }
      Historical@HistNjuv_NOS <- HistNjuv_NOS
    } else if (length(dim(Historical@HistNjuv_NOS)) == 3 && Bio@n_g == 1) {
      Historical <- check_array(Historical, "HistNjuv_NOS", dims = c(nsim, maxage, SOM@nyears+1))
      Historical@HistNjuv_NOS <- array(Historical@HistNjuv_NOS, c(dim(Historical@HistNjuv_NOS), 1))
    } else {
      Historical <- check_array(Historical, "HistNjuv_NOS", dims = c(nsim, maxage, SOM@nyears+1, Bio@n_g))
    }

    if (do_hatchery || has_strays) {

      if (!length(Historical@HistSpawner_HOS)) {
        init_HOS <- ifelse(do_hatchery, 1000, 0)
      } else if (!is.array(Historical@HistSpawner_HOS) && length(Historical@HistSpawner_HOS) == 1) {
        init_HOS <- Historical@HistSpawner_HOS
      } else {
        if (length(dim(Historical@HistSpawner_HOS)) == 3 && Hatchery@n_r == 1) {
          Historical <- check_array(Historical, "HistSpawner_HOS", dims = c(nsim, maxage, SOM@nyears))
          Historical@HistSpawner_HOS <- array(Historical@HistSpawner_HOS, c(dim(Historical@HistSpawner_HOS), 1))
        } else {
          Historical <- check_array(Historical, "HistSpawner_HOS", dims = c(nsim, maxage, SOM@nyears, Hatchery@n_r))
        }
      }

      if (!length(Historical@HistNjuv_HOS)) {
        if (SOM@nyears > maxage) stop("No historical abundance was provided for hatchery origin fish. Set historical years (SOM@nyears) < maxage")
        HistNjuv_HOS <- array(0, c(nsim, maxage, SOM@nyears + 1, Hatchery@n_r))

        HistNjuv_HOS[, maxage, SOM@nyears + 1, ] <- local({

          if (is.na(Hatchery@phatchery)) {
            init_Njuv <- init_HOS/(1 - Harvest@u_terminal)/(1 - Harvest@u_preterminal)/(1 - Hatchery@premove_HOS * Hatchery@m)
          } else {
            init_Njuv <- init_HOS/(1 - Harvest@u_terminal)/(1 - Harvest@u_preterminal)/
              (1 - Hatchery@phatchery)/(1 - Hatchery@premove_HOS * Hatchery@m)
          }
          outer(init_Njuv, 1/Hatchery@n_r) # matrix nsim x n_r
        })

        for (a in seq(maxage - 1, 1)) {
          y <- SOM@nyears + a - maxage + 1
          if (y > 0) {
            ZHOS <- matrix(Harvest@vulPT[, a] * Historical@HistFPT[, y, 2], nsim, Hatchery@n_r) + Hatchery@Mjuv_HOS[, a, y, ]
            HistNjuv_HOS[, a, y, ] <- HistNjuv_HOS[, a+1, y+1, ] * exp(ZHOS) /
              (1 - matrix(Hatchery@p_mature_HOS[, a, y, ], nsim, Hatchery@n_r))
          }
        }
        Historical@HistNjuv_HOS <- HistNjuv_HOS
      } else if (length(dim(Historical@HistNjuv_HOS)) == 3 && Hatchery@n_r == 1) {
        Historical <- check_array(Historical, "HistNjuv_HOS", dims = c(nsim, maxage, SOM@nyears+1))
        Historical@HistNjuv_HOS <- array(Historical@HistNjuv_HOS, c(dim(Historical@HistNjuv_HOS), 1))
      } else {
        Historical <- check_array(Historical, "HistNjuv_HOS", dims = c(nsim, maxage, SOM@nyears+1, Hatchery@n_r))
      }
    }

    SOM@Bio[[s]] <- Bio
    SOM@Habitat[[s]] <- Habitat
    SOM@Hatchery[[s]] <- Hatchery
    SOM@Harvest[[s]] <- Harvest
    SOM@Historical[[s]] <- Historical
  }

  # Check all maxage is identical
  maxage_s <- unique(vapply(SOM@Bio, slot, numeric(1), "maxage"))
  if (length(maxage_s) > 1) stop("Max age should be identical for all populations")

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

  if (length(slot(object, name)) != maxage) {
    stop(paste0("Need ", object_name, "@", name, " (maxage vector)"))
  }
  return(invisible(object))
}

check_maxage2matrix <- function(object, name, maxage, nsim) {
  object_name <- as.character(substitute(object))

  if (!is.array(slot(object, name))) {
    if (length(slot(object, name)) != maxage) {
      stop(paste0(object_name, "@", name, " needs to be length maxage (", maxage, ")"))
    }
    slot(object, name) <- matrix(slot(object, name), nsim, maxage, byrow = TRUE)
  }

  dim_i <- dim(slot(object, name))
  dim_check <- length(dim_i) == 2 && all(dim(slot(object, name)) == c(nsim, maxage))
  if (!dim_check) {
    stop(
      paste0(object_name, "@", name, " must be a matrix with dimension ",
             paste(c(nsim, maxage), collapse = ", "))
    )
  }
  return(invisible(object))
}

check_maxage2array <- function(object, name, maxage, nsim, years) {
  object_name <- as.character(substitute(object))

  if (!is.array(slot(object, name))) {
    if (length(slot(object, name)) != maxage) {
      stop(paste0(object_name, "@", name, " needs to be length maxage (", maxage, ")"))
    }
    slot(object, name) <- array(slot(object, name), c(maxage, nsim, years)) %>%
      aperm(c(2, 1, 3))
  }

  dim_i <- dim(slot(object, name))
  dim_check <- length(dim_i) == 3 && all(dim(slot(object, name)) == c(nsim, maxage, years))
  if (!dim_check) {
    stop(
      paste0(object_name, "@", name, " must be an array with dimension ",
             paste(c(nsim, maxage, years), collapse = ", "))
    )
  }
  return(invisible(object))
}

check_maxage2garray <- function(object, name, maxage, nsim, years, n_g) {
  object_name <- as.character(substitute(object))

  if (!is.array(slot(object, name))) {
    if (length(slot(object, name)) != maxage) {
      stop(paste0(object_name, "@", name, " needs to be length maxage (", maxage, ")"))
    }
    slot(object, name) <- array(slot(object, name), c(maxage, nsim, years, n_g)) %>%
      aperm(c(2, 1, 3, 4))
  }

  dim_i <- dim(slot(object, name))
  dim_check <- length(dim_i) == 4 && all(dim(slot(object, name)) == c(nsim, maxage, years, n_g))
  if (!dim_check) {
    stop(
      paste0(object_name, "@", name, " must be an array with dimension ",
             paste(c(nsim, maxage, years, n_g), collapse = ", "))
    )
  }
  return(invisible(object))
}

check_array <- function(object, name, dims) {
  object_name <- as.character(substitute(object))

  dim_i <- dim(slot(object, name))

  dim_check <- length(dim_i) == length(dims) && all(dim(slot(object, name)) == dims)
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
  if (length(slot(object, name)) != nsim) {
    stop(paste0(object_name, "@", name, " needs to length nsim (", nsim, ")"))
  }
  return(invisible(object))
}


solve_Herm <- function(h, p1 = 0) {
  p <- numeric(length(h))
  p[1] <- p1

  for (a in 2:length(p)) p[a] <- 1 - (1 - h[a])*(1 - p[a-1])
  #all(!MSEtool:::hrate(p) - h)
  return(p)
}
