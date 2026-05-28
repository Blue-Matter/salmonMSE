
NAor0 <- function(x) !length(x) || all(is.na(x)) || all(!x)

#' Check inputs to SOM object
#'
#' Ensures that the slots in the [salmonMSE::SOM-class] object have the correct dimensions. Function will
#' update some slots to their full dimensions.
#'
#' @param SOM [salmonMSE::SOM-class] object
#' @param silent Logical, whether to report progress in console
#' @returns Updated [salmonMSE::SOM-class] object with full dimensions in various slots
#' @examples
#' SOM_checked <- check_SOM(simple_SOM, silent = TRUE)
#'
#' @export
check_SOM <- function(SOM, silent = FALSE) {

  ### Base slots ----
  if (!length(SOM@nsim)) stop("Need SOM@nsim")
  #if (!length(SOM@nyears)) stop("Need SOM@nyears")
  if (!length(SOM@proyears)) stop("Need SOM@proyears")
  if (!length(SOM@seed)) stop("Need SOM@seed")
  nsim <- SOM@nsim
  if (nsim < 2) stop("Need SOM@nsim > 1")
  proyears <- SOM@proyears

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
    Bio <- update_S4(SOM@Bio[[s]])

    Bio <- check_numeric(Bio, "maxage")
    maxage <- Bio@maxage

    Bio <- check_numeric(Bio, "n_g", default = 1)
    Bio <- check_numeric(Bio, "p_LHG", size = Bio@n_g, default = rep(1/Bio@n_g, Bio@n_g))

    Bio <- check_maxage2array(Bio, "p_mature", maxage, nsim, proyears)
    if (length(Bio@p_female) == 1) Bio@p_female <- rep(Bio@p_female, maxage)
    Bio <- check_numeric(Bio, "p_female", size = maxage, default = rep(0.5, maxage))
    Bio <- check_numeric(Bio, "s_enroute", default = 1)
    Bio <- check_maxage2array(Bio, "fec", maxage, nsim, proyears)

    ### Check Habitat (initial) ----
    Habitat <- update_S4(SOM@Habitat[[s]])
    Habitat <- check_numeric(Habitat, "use_habitat", default = FALSE)

    if (length(dim(Bio@Mjuv_NOS)) == 3 && Bio@n_g == 1) {
      Bio <- check_maxage2array(Bio, "Mjuv_NOS", maxage-1, nsim, proyears)
      Bio@Mjuv_NOS <- array(Bio@Mjuv_NOS, c(dim(Bio@Mjuv_NOS), Bio@n_g))
    }
    Bio <- check_maxage2garray(Bio, "Mjuv_NOS", maxage-1, nsim, proyears, Bio@n_g)

    if (!Habitat@use_habitat) {
      # SRR pars
      slot(Bio, "SRrel") <- match.arg(slot(Bio, "SRrel"), choices = c("BH", "Ricker"))
      Bio <- check_numeric2nsim(Bio, "kappa", nsim)

      if (Bio@SRrel == "BH") {
        Bio <- check_numeric2nsim(Bio, "capacity", nsim)
      } else {
        Bio <- check_numeric2nsim(Bio, "Smax", nsim)
      }

      if (!length(Bio@phi) && length(Bio@tau)) {
        stop("Bio@tau is specified but not Bio@phi. Should specify both or none.")
      }
      if (length(Bio@phi) && !length(Bio@tau)) {
        stop("Bio@phi is specified but not Bio@tau. Should specify both or none.")
      }

      if (!length(Bio@phi)) {
        Bio@phi <- sapply(1:SOM@nsim, function(x) {
          calc_phi(
            Mjuv = matrix(Bio@Mjuv_NOS[x, , 1, ], Bio@maxage-1, Bio@n_g),
            p_mature = matrix(Bio@p_mature[x, , 1], Bio@maxage, Bio@n_g),
            p_female = Bio@p_female,
            fec = matrix(Bio@fec[x, , 1], Bio@maxage, Bio@n_g),
            s_enroute = Bio@s_enroute,
            n_g = Bio@n_g,
            p_LHG = Bio@p_LHG
          )
        })
      }
      if (!length(Bio@tau)) {
        Bio@tau <- sapply(1:SOM@nsim, function(x) {
          calc_phi(
            Mjuv = matrix(Bio@Mjuv_NOS[x, , 1, ], Bio@maxage-1, Bio@n_g),
            p_mature = matrix(Bio@p_mature[x, , 1], Bio@maxage, Bio@n_g),
            p_female = Bio@p_female,
            #fec = matrix(Bio@fec[x, , 1], Bio@maxage, Bio@n_g),
            s_enroute = Bio@s_enroute,
            n_g = Bio@n_g,
            p_LHG = Bio@p_LHG,
            output = "spawner"
          )
        })
      }
      Bio <- check_numeric2nsim(Bio, "phi", nsim)
      Bio <- check_numeric2nsim(Bio, "tau", nsim)
    } else {
      Habitat <- check_numeric(Habitat, "prespawn_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "prespawn_prod", default = 1)
      Habitat <- check_numeric(Habitat, "prespawn_capacity", default = Inf)

      Habitat <- check_numeric(Habitat, "egg_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "egg_prod", default = 1)
      Habitat <- check_numeric(Habitat, "egg_capacity", default = Inf)

      Habitat <- check_numeric(Habitat, "fry_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "fry_prod", default = 0.4)
      Habitat <- check_numeric(Habitat, "fry_capacity", default = Inf)
      if (!length(Habitat@fry_sdev)) {
        Habitat@fry_sdev <- matrix(1, nsim, proyears)
      }
      Habitat <- check_maxage2matrix(Habitat, "fry_sdev", proyears, nsim)

      Habitat <- check_numeric(Habitat, "smolt_rel", default = "BH")
      Habitat <- check_numeric(Habitat, "smolt_prod", default = 1)
      Habitat <- check_numeric(Habitat, "smolt_capacity", default = Inf)
      if (!length(Habitat@smolt_sdev)) {
        Habitat@smolt_sdev <- matrix(1, nsim, proyears)
      }
      Habitat <- check_maxage2matrix(Habitat, "smolt_sdev", proyears, nsim)
    }

    ### Check Hatchery ----
    Hatchery <- update_S4(SOM@Hatchery[[s]])

    Hatchery <- check_numeric(Hatchery, "n_r", default = 1)
    Hatchery <- check_numeric(Hatchery, "n_yearling", size = Hatchery@n_r, default = rep(0, Hatchery@n_r))
    Hatchery <- check_numeric(Hatchery, "n_subyearling", size = Hatchery@n_r, default = rep(0, Hatchery@n_r))
    if (!length(Hatchery@yearling_DD)) Hatchery@yearling_DD <- FALSE
    if (!length(Hatchery@subyearling_DD)) Hatchery@subyearling_DD <- FALSE

    do_hatchery <- sum(Hatchery@n_yearling, Hatchery@n_subyearling) > 0

    if (!length(Hatchery@stray_external)) Hatchery@stray_external <- matrix(0, maxage, Hatchery@n_r)
    Hatchery <- check_maxage2matrix(Hatchery, "stray_external", Hatchery@n_r, maxage)

    has_strays <- any(SOM@stray[-s, s] > 0) || sum(Hatchery@stray_external)

    if (do_hatchery || has_strays) {
      Hatchery <- check_numeric(Hatchery, "s_prespawn", default = 1)
      Hatchery <- check_numeric(Hatchery, "s_egg_smolt", default = 1)
      Hatchery <- check_numeric(Hatchery, "s_egg_subyearling", default = 1)

      if (length(dim(Hatchery@Mjuv_HOS)) == 3 && Hatchery@n_r == 1) {
        Hatchery <- check_maxage2array(Hatchery, "Mjuv_HOS", maxage-1, nsim, proyears)
        Hatchery@Mjuv_HOS <- array(Hatchery@Mjuv_HOS, c(dim(Hatchery@Mjuv_HOS), 1))
      }
      Hatchery <- check_maxage2garray(Hatchery, "Mjuv_HOS", maxage-1, nsim, proyears, Hatchery@n_r)

      if (!length(Hatchery@p_mature_HOS)) {
        Hatchery@p_mature_HOS <- array(Bio@p_mature, c(nsim, maxage, proyears, Hatchery@n_r))
      }
      if (length(dim(Hatchery@p_mature_HOS)) == 3 && Hatchery@n_r == 1) {
        Hatchery <- check_maxage2array(Hatchery, "p_mature_HOS", maxage, nsim, proyears)
        Hatchery@p_mature_HOS <- array(Hatchery@p_mature_HOS, c(nsim, maxage, proyears, 1))
      }
      Hatchery <- check_maxage2garray(Hatchery, "p_mature_HOS", maxage, nsim, proyears, Hatchery@n_r)
      Hatchery <- check_numeric(Hatchery, "gamma", default = 1)

      Hatchery <- check_numeric(Hatchery, "m", default = 0)

      Hatchery <- check_numeric(Hatchery, "brood_import", size = maxage, default = rep(0, maxage))
      Hatchery <- check_numeric(Hatchery, "pmax_esc", default = 0.75)
      Hatchery <- check_numeric(Hatchery, "pmax_NOB", default = 1)
      Hatchery <- check_numeric(Hatchery, "ptarget_NOB", default = 0.9)

      if (!length(Hatchery@fec_brood)) Hatchery@fec_brood <- Bio@fec
      Hatchery <- check_maxage2array(Hatchery, "fec_brood", maxage, nsim, proyears)

      if (!length(Hatchery@p_female_brood)) Hatchery@p_female_brood <- Bio@p_female
      Hatchery <- check_numeric(Hatchery, "p_female_brood", size = maxage, default = Bio@p_female)

      Hatchery <- check_numeric(Hatchery, "phatchery", default = NA)

      if (!is.function(Hatchery@premove_HOS) && is.numeric(Hatchery@premove_HOS)) {
        Hatchery <- check_numeric(Hatchery, "premove_HOS", default = 0)
      }
      if (!is.function(Hatchery@premove_NOS) && is.numeric(Hatchery@premove_NOS)) {
        Hatchery <- check_numeric(Hatchery, "premove_NOS", default = 0)
      }

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
      Hatchery <- check_numeric(Hatchery, "premove_HOS", default = 0)
      Hatchery <- check_numeric(Hatchery, "premove_NOS", default = 0)
      Hatchery <- check_numeric(Hatchery, "gamma", default = 1)
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
    Harvest <- update_S4(SOM@Harvest[[s]])

    Harvest <- check_numeric(Harvest, "type_PT", default = "u")
    Harvest <- check_numeric(Harvest, "type_T", default = "u")

    Harvest@type_PT <- match.arg(Harvest@type_PT, choices = c("u", "catch"))
    Harvest@type_T <- match.arg(Harvest@type_T, choices = c("u", "catch"))

    if (Harvest@type_PT == "u") {
      if (length(Harvest@u_preterminal) == 1) {
        Harvest <- check_numeric(Harvest, "u_preterminal", default = 0)
      } else if (is.matrix(Harvest@u_preterminal)) {
        Harvest <- check_maxage2matrix(Harvest, "u_preterminal", nsim, proyears)
      }
      Harvest <- check_numeric(Harvest, "K_PT", default = NA_real_)
    } else {
      Harvest <- check_numeric(Harvest, "u_preterminal", default = NA_real_)
      Harvest <- check_numeric(Harvest, "K_PT", default = 0)
    }

    if (Harvest@type_T == "u") {
      if (length(Harvest@u_terminal) == 1) {
        Harvest <- check_numeric(Harvest, "u_terminal", default = 0)
      } else if (is.matrix(Harvest@u_terminal)) {
        Harvest <- check_maxage2matrix(Harvest, "u_terminal", nsim, proyears)
      }
      Harvest <- check_numeric(Harvest, "K_T", default = NA_real_)
    } else {
      Harvest <- check_numeric(Harvest, "u_terminal", default = NA_real_)
      Harvest <- check_numeric(Harvest, "K_T", default = 0)
    }

    Harvest <- check_numeric(Harvest, "MSF_PT", default = FALSE)
    Harvest <- check_numeric(Harvest, "MSF_T", default = FALSE)

    Harvest <- check_numeric(Harvest, "release_mort", size = 2, default = c(0, 0))

    if (!length(Harvest@vulPT)) {
      if (NAor0(Harvest@u_preterminal > 0) || NAor0(Harvest@K_PT)) {
        Harvest@vulPT <- matrix(1, nsim, maxage)
      }
    }
    if (!length(Harvest@vulT)) {
      if (NAor0(Harvest@u_terminal > 0) || NAor0(Harvest@K_T)) {
        Harvest@vulT <- matrix(1, nsim, maxage)
      }
    }

    Harvest <- check_maxage2matrix(Harvest, "vulPT", maxage, nsim)
    Harvest <- check_maxage2matrix(Harvest, "vulT", maxage, nsim)

    ### Check Historical (initial) ----
    Historical <- update_S4(SOM@Historical[[s]])

    Njuv_NOS <- Njuv_HOS <- NA_real_

    if (!length(Historical@InitNjuv_NOS)) {
      Njuv_NOS <- 1000
    } else if (length(Historical@InitNjuv_NOS) == 1) {
      Njuv_NOS <- Historical@InitNjuv_NOS
    }

    if (!length(Historical@InitNjuv_NOS)) {
      Njuv_HOS <- 1000
    } else if (length(Historical@InitNjuv_HOS) == 1) {
      Njuv_HOS <- Historical@InitNjuv_HOS
    }

    if (!is.array(Historical@InitNjuv_NOS)) {
      Historical@InitNjuv_NOS <- array(0, c(nsim, maxage, Bio@n_g))
      Historical@InitNjuv_NOS[, maxage, ] <- Njuv_NOS/Bio@n_g
    }
    if (!is.array(Historical@InitNjuv_HOS)) {
      Historical@InitNjuv_HOS <- array(0, c(nsim, maxage, Hatchery@n_r))
      Historical@InitNjuv_HOS[, maxage, ] <- Njuv_HOS/Hatchery@n_r
    }
    Historical <- check_array(Historical, "InitNjuv_NOS", dims = c(nsim, maxage, Bio@n_g))
    Historical <- check_array(Historical, "InitNjuv_HOS", dims = c(nsim, maxage, Hatchery@n_r))

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


#' @importFrom methods .slotNames .hasSlot slot
update_S4 <- function(x) {
  snames <- .slotNames(class(x))
  slot_check <- sapply(snames, function(i) .hasSlot(x, i))

  if (any(!slot_check)) {
    xnew <- new(class(x))
    for (i in snames) {
      if (slot_check[i]) slot(xnew, i) <- slot(x, i)
    }
    return(xnew)
  } else {
    return(x)
  }
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
    if (output == "egg" && !is.matrix(fec)) fec <- matrix(fec, ncol = 1)
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

