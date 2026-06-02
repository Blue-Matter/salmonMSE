

#' Run salmonMSE
#'
#' @description
#' `salmonMSE()` runs a salmon management strategy evaluation from an operating model object (\linkS4class{SOM}), by checking
#' the operating model object with `check_SOM()`, running the projection in `ProjectSOM()` (parallel if called upon, then stitches
#' together the output in a single object), and calculates reference points with `calc_ref()`.
#' @param SOM An object of class \linkS4class{SOM}
#' @param ncores Integer, maximum number of processors to run projection with parallel processing
#' @param silent Logical, whether to report progress in console
#' @return
#' \linkS4class{SMSE} object
#' @examples
#' \dontrun{
#' SMSE <- salmonMSE(simple_SOM)
#' }
#'
#' @export
#' @importFrom parallel detectCores makeCluster stopCluster parLapplyLB
#' @importFrom abind abind
salmonMSE <- function(SOM, ncores = 1, silent = FALSE) {

  SOM <- check_SOM(SOM, silent = silent)

  if (ncores == 1) {
    SMSE <- ProjectSOM(SOM, check = FALSE)
  } else {
    nits <- split_sims(SOM@nsim, ncores = min(ncores, parallel::detectCores()))
    cores <- length(nits)

    if (cores == 1) {
      if (!silent) message("Running projection on one core (parallel processing not needed)")
      SMSE <- ProjectSOM(SOM, check = FALSE)
    } else {
      if (!silent) message("Running ", nsim, " simulations in parallel on ", cores, " cores")
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))

      SMSE_list <- parallel::parLapplyLB(cl, X = nits, ProjectSOM_parallel, SOM = SOM)

      # Stitch objects together
      SMSE <- new("SMSE")
      vars <- slotNames("SMSE")
      for (j in vars) {
        if (j %in% c("Name", "proyears", "nstocks", "Snames")) {
          slot(SMSE, j) <- slot(SMSE_list[[1]], j)
        } else if (j == "nsim") {
          slot(SMSE, j) <- sapply(SMSE_list, slot, name = j) %>% sum()
        } else if (j != "Misc") {
          slot(SMSE, j) <- lapply(SMSE_list, slot, name = j) %>% abind::abind(along = 1) %>% `dimnames<-`(NULL)
        }
      }
      if (length(SMSE_list[[1]]@Misc$RS)) {
        vars_RS <- names(SMSE_list[[1]]@Misc$RS)
        SMSE@Misc$RS <- lapply(vars_RS, function(i) {
          lapply(SMSE_list, function(j) j@Misc$RS[[i]]) %>% abind::abind(along = 1) %>% `dimnames<-`(NULL)
        }) %>%
          structure(names = vars_RS)
      }
      if (length(SMSE_list[[1]]@Misc$LHG)) {
        vars_LHG <- names(SMSE_list[[1]]@Misc$LHG)
        SMSE@Misc$LHG <- lapply(vars_LHG, function(i) {
          lapply(SMSE_list, function(j) j@Misc$LHG[[i]]) %>% abind::abind(along = 1) %>% `dimnames<-`(NULL)
        }) %>%
          structure(names = vars_LHG)
      }
    }
  }

  SMSE@Misc$SOM <- SOM
  SMSE@Misc$Ref <- calc_ref(SOM, check = FALSE)

  return(SMSE)
}

define_hatchery_args <- function(SOM) {
  ns <- length(SOM)

  output_s <- lapply(1:ns, function(s) {

    Hatchery <- SOM@Hatchery[[s]]

    egg_yearling <- ifelse(sum(Hatchery@n_yearling) > 0, sum(Hatchery@n_yearling)/Hatchery@s_egg_smolt, 0)
    egg_subyearling <- ifelse(sum(Hatchery@n_subyearling) > 0, sum(Hatchery@n_subyearling)/Hatchery@s_egg_subyearling, 0)
    egg_target <- egg_yearling + egg_subyearling

    if (egg_target > 0) {
      p_yearling <- Hatchery@n_yearling/sum(Hatchery@n_yearling, Hatchery@n_subyearling) # Vector by release strategy
      p_subyearling <- Hatchery@n_subyearling/sum(Hatchery@n_yearling, Hatchery@n_subyearling) # Vector by release strategy
    } else {
      p_yearling <- p_subyearling <- 0
    }

    output <- list(
      egg_target = egg_target,
      premove_NOS = Hatchery@premove_NOS
    )

    has_strays <- any(SOM@stray[-s, s] > 0) || sum(Hatchery@stray_external)

    if (egg_target > 0) {
      output_with_hatchery <- list(
        f_brood = Hatchery@f_brood,
        pmax_esc = Hatchery@pmax_esc,
        ptarget_NOB = Hatchery@ptarget_NOB,
        pmax_NOB = Hatchery@pmax_NOB,
        brood_import = Hatchery@brood_import,

        phatchery = Hatchery@phatchery,

        fec_brood = Hatchery@fec_brood,
        s_prespawn = Hatchery@s_prespawn,
        p_female = Hatchery@p_female_brood,

        p_yearling = p_yearling,
        p_subyearling = p_subyearling,

        s_yearling = Hatchery@s_egg_smolt,
        s_subyearling = Hatchery@s_egg_subyearling,

        yearling_DD = Hatchery@yearling_DD,
        subyearling_DD = Hatchery@subyearling_DD,

        premove_HOS = Hatchery@premove_HOS,
        gamma = Hatchery@gamma
      )

      output <- c(output, output_with_hatchery)

    } else if (has_strays) {

      output_with_strays <- list(
        premove_HOS = Hatchery@premove_HOS,
        gamma = Hatchery@gamma
      )
      output <- c(output, output_with_strays)

    }
    return(output)
  })

  return(output_s)
}

define_habitat_args <- function(SOM) slot(SOM, "Habitat")

define_fitness_args <- function(SOM) {
  ns <- length(SOM)

  output_s <- lapply(1:ns, function(s) {
    Hatchery <- SOM@Hatchery[[s]]

    do_hatchery <- sum(Hatchery@n_yearling, Hatchery@n_subyearling) > 0
    has_strays <- any(SOM@stray[-s, s] > 0) || sum(Hatchery@stray_external)
    has_HOS <- do_hatchery || has_strays
    do_fitness <- any(Hatchery@fitness_type == "Ford")

    output <- list(fitness_type = Hatchery@fitness_type)

    if (has_HOS && do_fitness) {
      output <- list(
        fitness_type = Hatchery@fitness_type,
        rel_loss = Hatchery@rel_loss,
        phenotype_variance = Hatchery@phenotype_variance,
        fitness_variance = Hatchery@fitness_variance,
        fitness_floor = Hatchery@fitness_floor,
        heritability = Hatchery@heritability,
        theta = Hatchery@theta
      )
    }
    return(output)
  })

  return(output_s)
}

define_SRRpars <- function(SOM) {
  ns <- length(SOM)

  output_s <- lapply(1:ns, function(s) {
    df <- data.frame()

    if (!SOM@Habitat[[s]]@use_habitat) {
      Bio <- SOM@Bio[[s]]
      SRrel <- Bio@SRrel
      df <- data.frame(
        kappa = Bio@kappa,
        phi = Bio@phi,
        tau = Bio@tau,
        SRrel = SRrel
      )
      if (SRrel == "Ricker") {
        df$Smax <- Bio@Smax
      } else {
        df$capacity <- Bio@capacity
      }
    }
    return(df)
  })

  return(output_s)
}


sapply2 <- base::sapply
formals(sapply2)$simplify <- "array"

ProjectSOM_parallel <- function(X, SOM, check = FALSE) ProjectSOM(SOM, sims = X, check = check)

split_sims <- function(nsim, ncores) {
  nits_min <- 2

  nits_prelim <- rep(ceiling(nsim/ncores), ncores)
  nits_prelim[nits_prelim < nits_min] <- nits_min

  nits <- lapply(1:ncores, function(i) {
    prev <- ifelse(i == 1, 0, sum(nits_prelim[seq(1, i - 1)]))
    sims <- prev + seq(1, nits_prelim[i])
    sims[sims <= nsim]
  })
  nits_use <- sapply(nits, length) > 0
  nits <- nits[nits_use]

  if (length(nits[[length(nits)]]) < nits_min) {
    nits[[length(nits) - 1]] <- c(nits[[length(nits) - 1]], nits[[length(nits)]])
    nits <- nits[-length(nits)]
  }

  return(nits)
}

#' @name salmonMSE
#' @description `ProjectSOM()` is the internal projection function.
#' @param check Logical, whether to check the structure of the input object with [check_SOM()]
#' @param sims Optional integer vector to run projection for a subset of simulations. Intended for parallel processing.
#' @export
ProjectSOM <- function(SOM, sims, check = FALSE) {
  if (check) SOM <- check_SOM(SOM)
  if (missing(sims)) sims <- seq(1, SOM@nsim)
  if (any(sims) > SOM@nsim) stop("There are `sims` > SOM@nsim")
  if (length(sims) < 2) stop("Need at least two simulations in `sims`")

  # Variables
  ns <- length(SOM@Bio) # Number of stocks
  nage <- sapply(SOM@Bio, slot, "maxage") %>% unique()

  n_r <- sapply(SOM@Hatchery, slot, "n_r") %>% unique()
  n_g <- sapply(SOM@Bio, slot, "n_g") %>% unique()

  if (length(n_r) > 1) stop("Number of release strategies vary by population.")
  if (length(n_g) > 1) stop("Number of life history groups vary by population.")

  nsim <- length(sims)
  proyears <- SOM@proyears

  # Hatchery arguments
  hatchery_args <- define_hatchery_args(SOM)
  m <- sapply(SOM@Hatchery, slot, "m") # Mark rate, need to make sure check_SOM default is length 1

  # Stock recruit parameters
  SRRpars <- define_SRRpars(SOM)

  # Fitness arguments
  fitness_args <- define_fitness_args(SOM)

  # Freshwater functions and arguments
  habitat_args <- define_habitat_args(SOM)

  #### Arrays of state variables ----
  # Marine life stages, brood, egg production by age
  Njuv_NOS <- Return_NOS <- Escapement_NOS <- NOB <- NOS <- Egg_NOS <- array(NA_real_, c(nsim, ns, nage, proyears, n_g))
  Njuv_HOS <- Return_HOS <- Escapement_HOS <- HOB <- HOB_stray <-
    HOS <- HOS_stray <- HOS_effective <- Egg_HOS <- array(NA_real_, c(nsim, ns, nage, proyears, n_r))
  HOB_import <- array(NA_real_, c(nsim, ns, nage, proyears))

  # Early freshwater life stages
  Fry_NOS <- Smolt_NOS <- Fry_HOS <- Smolt_HOS <- array(NA_real_, c(nsim, ns, proyears, n_g))
  Rel <- Smolt_Rel <- array(NA_real_, c(nsim, ns, proyears, n_r))

  # In-river removals
  IRR_NOS <- array(NA_real_, c(nsim, ns, nage, proyears, n_g))
  IRR_HOS <- array(NA_real_, c(nsim, ns, nage, proyears, n_r))

  # Marine catch and exploitation rate
  KPT_NOS <- KT_NOS <- DPT_NOS <- DDPT_NOS <- DT_NOS <- DDT_NOS <-
    UPT_NOS <- UT_NOS <- ExPT_NOS <- ExT_NOS <- array(NA_real_, c(nsim, ns, nage, proyears, n_g))
  KPT_HOS <- KT_HOS <- DPT_HOS <- DDPT_HOS <- DT_HOS <- DDT_HOS <-
    UPT_HOS <- UT_HOS <- ExPT_HOS <- ExT_HOS <- array(NA_real_, c(nsim, ns, nage, proyears, n_r))

  pNOB <- pHOS_census <- pHOS_effective <- PNI <- p_wild <- array(NA_real_, c(nsim, ns, proyears))

  zbar <- fitness <- array(NA_real_, c(nsim, ns, 2, proyears)) # proyears indexes brood year
  zbar_brood <- array(NA_real_, c(nsim, ns, nage, 2, proyears)) # proyears indexes return year
  fitness_loss <- array(NA_real_, c(nsim, ns, proyears, 2, 3)) # By brood year

  Mjuv_loss_NOS <- array(NA_real_, c(nsim, ns, nage-1, proyears, n_g))
  Mjuv_loss_HOS <- array(NA_real_, c(nsim, ns, nage-1, proyears, n_r))

  #### Data objects from SOM ----
  do_hatchery <- sapply(1:ns, function(s) hatchery_args[[s]]$egg_target > 0)
  has_strays <- sapply(1:ns, function(s) {
    any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)
  })

  # Maturity
  p_mature_NOS <- sapply2(1:ns, function(s) SOM@Bio[[s]]@p_mature[sims, , , drop = FALSE]) %>%
    array(c(nsim, nage, proyears, n_g, ns)) %>%
    aperm(c(1, 5, 2, 3, 4))
  p_mature_HOS <- sapply2(1:ns, function(s) {
    if (do_hatchery[s]) {
      SOM@Hatchery[[s]]@p_mature_HOS[sims, , , , drop = FALSE]
    } else {
      array(0, c(nsim, nage, proyears, n_r))
    }
  }) %>%
    aperm(c(1, 5, 2, 3, 4))

  # Juvenile natural mortality
  Mjuv_NOS <- array(NA_real_, c(nsim, ns, nage-1, proyears, n_g))
  Mjuv_HOS <- array(NA_real_, c(nsim, ns, nage-1, proyears, n_r))
  Mjuv_NOS[] <- sapply2(1:ns, function(s) SOM@Bio[[s]]@Mjuv_NOS[sims, , , , drop = FALSE]) %>%
    aperm(c(1, 5, 2, 3, 4))
  Mjuv_HOS[] <- sapply2(1:ns, function(s) {
    if (do_hatchery[s]) {
      SOM@Hatchery[[s]]@Mjuv_HOS[sims, , , , drop = FALSE]
    } else {
      array(0, c(nsim, nage-1, proyears, n_r))
    }
  }) %>%
    aperm(c(1, 5, 2, 3, 4))

  # Harvest/fishery settings - all from the first harvest object
  type_PT <- SOM@Harvest[[1]]@type_PT
  type_T <- SOM@Harvest[[1]]@type_T

  u_preterminal <- SOM@Harvest[[1]]@u_preterminal
  u_terminal <- SOM@Harvest[[1]]@u_terminal

  K_PT <- SOM@Harvest[[1]]@K_PT
  K_T <- SOM@Harvest[[1]]@K_T

  MSF_PT <- SOM@Harvest[[1]]@MSF_PT
  MSF_T <- SOM@Harvest[[1]]@MSF_T

  # Fishery vulnerability
  vulPT <- vulT <- array(NA_real_, c(nsim, ns, nage))
  vulPT[] <- sapply2(1:ns, function(s) SOM@Harvest[[s]]@vulPT[sims, ]) %>% aperm(c(1, 3, 2))
  vulT[] <- sapply2(1:ns, function(s) SOM@Harvest[[s]]@vulT[sims, ]) %>% aperm(c(1, 3, 2))

  # Release mortality
  release_mort <- array(NA_real_, c(2, ns))
  release_mort[] <- sapply(SOM@Harvest, slot, "release_mort") # Need to make sure check_SOM default is length 2

  # En-route survival
  s_enroute <- sapply(SOM@Bio, slot, "s_enroute")

  # External strays
  stray_external <- array(NA_real_, c(ns, nage, n_r))
  stray_external[] <- sapply2(SOM@Hatchery, slot, "stray_external") %>%
    aperm(c(3, 1, 2))

  #### Initialize population ----
  Njuv_NOS[, , , 1, ] <- sapply2(1:ns, function(s) SOM@Historical[[s]]@InitNjuv_NOS[sims, , , drop = FALSE]) %>%
    aperm(c(1, 4, 2, 3))
  Njuv_HOS[, , , 1, ] <- sapply2(1:ns, function(s) {
    if (do_hatchery[s]) {
      SOM@Historical[[s]]@InitNjuv_HOS[sims, , , drop = FALSE]
    } else {
      array(0, c(nsim, nage, n_r))
    }
  }) %>%
    aperm(c(1, 4, 2, 3))

  #### Projection ----
  for (y in 1:proyears) {

    # Preterminal catch - harvest management acts upon on all stocks simultaneously
    PT_Calcs <- lapply(1:nsim, function(x) {
      catch_func(
        NO = array(Njuv_NOS[x, , , y, ], c(ns, nage, n_g)),
        HO = array(Njuv_HOS[x, , , y, ], c(ns, nage, n_r)),
        type = type_PT,
        U = u_preterminal,
        K = K_PT,
        V = matrix(vulPT[x, , ], ns, nage),
        m = m,
        MSF = MSF_PT, # Need to make sure check_SOM default is length 1
        release_mort = release_mort[1, ]
      )
    })
    vars <- names(PT_Calcs[[1]])
    PT_Calcs_y <- lapply(vars, function(i) {
      sapply2(PT_Calcs, getElement, i)
    }) %>%
      structure(names = vars)

    KPT_NOS[, , , y, ] <- aperm(PT_Calcs_y$K_NO, c(4, 1:3))
    DPT_NOS[, , , y, ] <- aperm(PT_Calcs_y$D_NO, c(4, 1:3))
    DDPT_NOS[, , , y, ] <- aperm(PT_Calcs_y$DD_NO, c(4, 1:3))

    KPT_HOS[, , , y, ] <- aperm(PT_Calcs_y$K_HO, c(4, 1:3))
    DPT_HOS[, , , y, ] <- aperm(PT_Calcs_y$D_HO, c(4, 1:3))
    DDPT_HOS[, , , y, ] <- aperm(PT_Calcs_y$DD_HO, c(4, 1:3))

    UPT_NOS[, , , y, ] <- aperm(PT_Calcs_y$U_NO, c(4, 1:3))
    ExPT_NOS[, , , y, ] <- aperm(PT_Calcs_y$Ex_NO, c(4, 1:3))

    UPT_HOS[, , , y, ] <- aperm(PT_Calcs_y$U_HO, c(4, 1:3))
    ExPT_HOS[, , , y, ] <- aperm(PT_Calcs_y$Ex_HO, c(4, 1:3))

    # Maturity (begin second half)
    Return_NOS[, , , y, ] <- Njuv_NOS[, , , y, ] * (1 - ExPT_NOS[, , , y, ]) * p_mature_NOS[, , , y, ]
    Return_HOS[, , , y, ] <- Njuv_HOS[, , , y, ] * (1 - ExPT_HOS[, , , y, ]) * p_mature_HOS[, , , y, ]

    # Terminal marine catch - harvest management acts upon on all stocks simultaneously
    T_Calcs <- lapply(1:nsim, function(x) {
      catch_func(
        NO = array(Return_NOS[x, , , y, ], c(ns, nage, n_g)),
        HO = array(Return_HOS[x, , , y, ], c(ns, nage, n_r)),
        type = type_T,
        U = u_terminal,
        K = K_T,
        V = matrix(vulT[x, , ], ns, nage),
        m = m,
        MSF = MSF_T,
        release_mort = release_mort[2, ]
      )
    })
    vars <- names(T_Calcs[[1]])
    T_Calcs_y <- lapply(vars, function(i) {
      sapply2(T_Calcs, getElement, i)
    }) %>%
      structure(names = vars)

    KT_NOS[, , , y, ] <- aperm(T_Calcs_y$K_NO, c(4, 1:3))
    DT_NOS[, , , y, ] <- aperm(T_Calcs_y$D_NO, c(4, 1:3))
    DDT_NOS[, , , y, ] <- aperm(T_Calcs_y$DD_NO, c(4, 1:3))

    KT_HOS[, , , y, ] <- aperm(T_Calcs_y$K_HO, c(4, 1:3))
    DT_HOS[, , , y, ] <- aperm(T_Calcs_y$D_HO, c(4, 1:3))
    DDT_HOS[, , , y, ] <- aperm(T_Calcs_y$DD_HO, c(4, 1:3))

    UT_NOS[, , , y, ] <- aperm(T_Calcs_y$U_NO, c(4, 1:3))
    ExT_NOS[, , , y, ] <- aperm(T_Calcs_y$Ex_NO, c(4, 1:3))

    UT_HOS[, , , y, ] <- aperm(T_Calcs_y$U_HO, c(4, 1:3))
    ExT_HOS[, , , y, ] <- aperm(T_Calcs_y$Ex_HO, c(4, 1:3))

    # Escapement from marine fisheries
    Escapement_NOS[, , , y, ] <- Return_NOS[, , , y, ] * (1 - ExT_NOS[, , , y, ])
    Escapement_HOS[, , , y, ] <- Return_HOS[, , , y, ] * (1 - ExT_HOS[, , , y, ])

    # Move strays (internally)
    Stray_Calcs <- lapply(1:nsim, function(x) {
      stray_func(
        N = array(Escapement_HOS[x, , , y, ], c(ns, nage, n_r)),
        stray_matrix = SOM@stray,
        m = m
      )
    })

    # Straying, in-river return, brood, egg production, outmigrating in next year
    for (s in 1:ns) {

      # Mean phenotype by brood year of parents
      if (any(SOM@Hatchery[[s]]@fitness_type == "Ford") && (do_hatchery[s] || has_strays[s])) {
        for (a in 1:nage) {
          t <- y - a
          if (t <= 0) {
            zbar_brood[, s, a, , y] <- SOM@Hatchery[[s]]@zbar_start[sims, abs(t) + 1, ]
          } else {
            zbar_brood[, s, a, , y] <- zbar[, s, , t]
          }
        }
      }

      FW_Calcs <- lapply(1:nsim, function(x) {
        xx <- sims[x]

        # Calculate recipient strays (internal and external) and their mark rate
        Nage_stray <- matrix(stray_external[s, , ] + Stray_Calcs[[x]]$N_stray[s, , ], nage, n_r)
        m_stray <- sum(Stray_Calcs[[x]]$m_stray * Stray_Calcs[[x]]$N_stray, stray_external[s, , ])/sum(Nage_stray)
        m_stray[is.na(m_stray)] <- 0

        # Calculate broodtake, in-river removals, and spawners arriving at spawning grounds
        .hatchery_args <- hatchery_args[[s]]
        if (.hatchery_args$egg_target > 0) .hatchery_args$fec_brood <- hatchery_args[[s]]$fec_brood[xx, , y]

        .fitness_args <- fitness_args[[s]]
        if (!is.null(.fitness_args$heritability)) {
          .fitness_args$heritability <- fitness_args[[s]]$heritability[xx]
        }

        if (SOM@Habitat[[s]]@use_habitat) {
          .habitat_args <- habitat_args[[s]]
          .habitat_args@fry_sdev <- habitat_args[[s]]@fry_sdev[xx, y, drop = FALSE]
          .habitat_args@smolt_sdev <- habitat_args[[s]]@smolt_sdev[xx, y, drop = FALSE]
          SRRpars_x <- data.frame()
        } else {
          .habitat_args <- list()
          SRRpars_x <- SRRpars[[s]][xx, ]
        }

        FW_func(
          Nage_NOS = matrix(Escapement_NOS[x, s, , y, ], nage, n_g),
          Nage_HOS = matrix(Stray_Calcs[[x]]$N_remain[s, , ], nage, n_r),
          Nage_stray = Nage_stray,
          m = m[s],
          m_stray = m_stray,
          s_enroute = s_enroute[s],
          hatchery_args = .hatchery_args,
          zbar_brood = zbar_brood[x, s, , , y],
          fitness_args = .fitness_args,
          fec = SOM@Bio[[s]]@fec[xx, , y],
          p_female = SOM@Bio[[s]]@p_female,
          habitat_args = .habitat_args,
          SRRpars = SRRpars_x,
          p_LHG = SOM@Bio[[s]]@p_LHG
        )
      })

      # Assign FW_Calcs output to global variables
      vars <- names(FW_Calcs[[1]])
      FW_Calcs_y <- lapply(vars, function(i) sapply2(FW_Calcs, getElement, i)) %>%
        structure(names = vars)

      NOB[, s, , y, ] <- aperm(FW_Calcs_y$NOB, c(3, 1, 2))
      HOB[, s, , y, ] <- aperm(FW_Calcs_y$HOB_unmarked + FW_Calcs_y$HOB_marked, c(3, 1, 2))
      HOB_stray[, s, , y, ] <- aperm(FW_Calcs_y$HOB_stray, c(3, 1, 2))
      HOB_import[, s, , y] <- t(FW_Calcs_y$HOB_import)

      pNOB[, s, y] <- FW_Calcs_y$pNOB

      IRR_NOS[, s, , y, ] <- aperm(FW_Calcs_y$NO_remove, c(3, 1, 2))
      IRR_HOS[, s, , y, ] <- aperm(FW_Calcs_y$HO_remove, c(3, 1, 2))

      NOS[, s, , y, ] <- aperm(FW_Calcs_y$NOS, c(3, 1, 2))
      HOS[, s, , y, ] <- aperm(FW_Calcs_y$HOS, c(3, 1, 2))
      HOS_effective[, s, , y, ] <- aperm(FW_Calcs_y$HOS_effective, c(3, 1, 2))
      HOS_stray[, s, , y, ] <- aperm(FW_Calcs_y$HOS_stray, c(3, 1, 2))

      Egg_NOS[, s, , y, ] <- aperm(FW_Calcs_y$Egg_NOS, c(3, 1, 2))
      Egg_HOS[, s, , y, ] <- aperm(FW_Calcs_y$Egg_HOS, c(3, 1, 2))

      pHOS_census[, s, y] <- FW_Calcs_y$pHOScensus
      pHOS_effective[, s, y] <- FW_Calcs_y$pHOSeff

      fitness[, s, , y] <- t(FW_Calcs_y$fitness)
      zbar[, s, , y] <- t(FW_Calcs_y$zbar)
      fitness_loss[, s, y, , ] <- aperm(FW_Calcs_y$fitness_loss, c(3, 1, 2))

      if (y < proyears) {
        Fry_NOS[, s, y+1, ] <- if (n_g == 1) FW_Calcs_y$Fry_NOS else t(FW_Calcs_y$Fry_NOS)
        Fry_HOS[, s, y+1, ] <- if (n_g == 1) FW_Calcs_y$Fry_HOS else t(FW_Calcs_y$Fry_HOS)

        Smolt_NOS[, s, y+1, ] <- if (n_g == 1) FW_Calcs_y$Smolt_NOS else t(FW_Calcs_y$Smolt_NOS)
        Smolt_HOS[, s, y+1, ] <- if (n_g == 1) FW_Calcs_y$Smolt_HOS else t(FW_Calcs_y$Smolt_HOS)

        Rel[, s, y+1, ] <- if (n_r == 1) FW_Calcs_y$yearling + FW_Calcs_y$subyearling else t(FW_Calcs_y$yearling + FW_Calcs_y$subyearling)
        Smolt_Rel[, s, y+1, ] <- if (n_r == 1) FW_Calcs_y$Smolt_RelOut else t(FW_Calcs_y$Smolt_RelOut)

        Njuv_NOS[, s, 1, y+1, ] <- Smolt_NOS[, s, y+1, ] + Smolt_HOS[, s, y+1, ]
        Njuv_HOS[, s, 1, y+1, ] <- Smolt_Rel[, s, y+1, ]
      }
    }

    # Advance juvenile age classes to next year with maturity and natural mortality, penalized by fitness loss
    if (y < proyears) {
      Njuv_NOS_midpoint <- array(NA, c(nsim, ns, nage, n_g))
      Njuv_HOS_midpoint <- array(NA, c(nsim, ns, nage, n_r))
      Njuv_NOS_midpoint[] <- Njuv_NOS[, , , y, ] * (1 - ExPT_NOS[, , , y, ]) * (1 - p_mature_NOS[, , , y, ])
      Njuv_HOS_midpoint[] <- Njuv_HOS[, , , y, ] * (1 - ExPT_HOS[, , , y, ]) * (1 - p_mature_HOS[, , , y, ])

      for (a in seq(1, nage-1)) {
        t <- y - a
        if (t <= 0) {
          Mjuv_loss_NOS[, , a, y, ] <- Mjuv_NOS[, , a, y, ]
          Mjuv_loss_HOS[, , a, y, ] <- Mjuv_HOS[, , a, y, ]
        } else {
          Mjuv_loss_NOS[, , a, y, ] <- local({
            .M <- array(Mjuv_NOS[, , a, y, ], c(nsim, ns, n_g))
            surv_fitness <- exp(-.M) * array(fitness_loss[, , t, 1, 3], c(nsim, ns, n_g))
            .M[.M > .Machine$double.eps] <- -log(surv_fitness[.M > .Machine$double.eps])
            .M
          })
          Mjuv_loss_HOS[, , a, y, ] <- local({
            .M <- array(Mjuv_HOS[, , a, y, ], c(nsim, ns, n_r))
            surv_fitness <- exp(-.M) * array(fitness_loss[, , t, 1, 3], c(nsim, ns, n_r))
            .M[.M > .Machine$double.eps] <- -log(surv_fitness[.M > .Machine$double.eps])
            .M
          })
        }
      }

      Njuv_NOS[, , -1, y+1, ] <- Njuv_NOS_midpoint[, , seq(1, nage-1), ] * exp(-Mjuv_loss_NOS[, , , y, ])
      Njuv_HOS[, , -1, y+1, ] <- Njuv_HOS_midpoint[, , seq(1, nage-1), ] * exp(-Mjuv_loss_HOS[, , , y, ])
    }
  }

  # State variables after the projection
  for (s in 1:ns) {
    p_wild[, s, ] <- calc_pwild_age(
      NOS_a = apply(NOS[, s, , , , drop = FALSE], c(1, 3, 4), sum),
      HOS_a = apply(HOS[, s, , , , drop = FALSE], c(1, 3, 4), sum),
      fec = SOM@Bio[[s]]@fec[sims, , seq(1, proyears)],
      gamma = SOM@Hatchery[[s]]@gamma
    )

    # PNI
    PNI[, s, ] <- local({

      has_NOS <- apply(NOS[, s, , , , drop = FALSE], c(1, 4), sum) > 0

      any_HOS <- sum(HOS[, s, , , ]) > 0
      any_NOB <- sum(NOB[, s, , , ]) > 0
      any_HOB <- sum(HOB[, s, , , ], HOB_stray[, s, , , ], HOB_import[, s, , ]) > 0

      .PNI <- matrix(NA_real_, nsim, proyears)

      if (!any_HOS && !any_HOB) {  # No HOS, no HOB (natural system)

        .PNI[has_NOS] <- 1

      } else if (!any_NOB && any_HOS) { # One-way gene flow: no NOB but has HOS

        h2 <- fitness_args[[s]]$heritability
        fitness_variance <- fitness_args[[s]]$fitness_variance
        if (!is.null(h2) && !is.null(fitness_variance)) {
          .PNI[] <- h2/(h2 + (1 - h2 + fitness_variance) * pHOS_effective[, s, ])
        }

      } else {

        .PNI[] <- pNOB[, s, ]/(pNOB[, s, ] + pHOS_effective[, s, ])

      }
      return(.PNI)
    })

  }

  # Harvest rate and exploitation rate aggregated across life history groups and release strategies
  UPT_NOS_ <- apply(KPT_NOS, 1:4, sum)/apply(Njuv_NOS, 1:4, sum)
  UPT_NOS_[is.na(UPT_NOS_)] <- 0

  UT_NOS_ <- apply(KT_NOS, 1:4, sum)/apply(Return_NOS, 1:4, sum)
  UT_NOS_[is.na(UT_NOS_)] <- 0

  UPT_HOS_ <- apply(KPT_HOS, 1:4, sum)/apply(Njuv_HOS, 1:4, sum)
  UPT_HOS_[is.na(UPT_HOS_)] <- 0

  UT_HOS_ <- apply(KT_HOS, 1:4, sum)/apply(Return_HOS, 1:4, sum)
  UT_HOS_[is.na(UT_HOS_)] <- 0

  ExPT_NOS_ <- apply(KPT_NOS + DDPT_NOS, 1:4, sum)/apply(Njuv_NOS, 1:4, sum)
  ExPT_NOS_[is.na(ExPT_NOS_)] <- 0

  ExT_NOS_ <- apply(KT_NOS + DDT_NOS, 1:4, sum)/apply(Return_NOS, 1:4, sum)
  ExT_NOS_[is.na(ExT_NOS_)] <- 0

  ExPT_HOS_ <- apply(KPT_HOS + DDPT_HOS, 1:4, sum)/apply(Njuv_HOS, 1:4, sum)
  ExPT_HOS_[is.na(ExPT_HOS_)] <- 0

  ExT_HOS_ <- apply(KT_HOS + DDT_HOS, 1:4, sum)/apply(Return_HOS, 1:4, sum)
  ExT_HOS_[is.na(ExT_HOS_)] <- 0

  # Output
  SMSE <- new(
    "SMSE",
    Name = SOM@Name,
    proyears = proyears,
    nsim = nsim,
    nstocks = ns,
    Snames = sapply(1:ns, function(s) if (length(SOM@Bio[[s]]@Name)) SOM@Bio[[s]]@Name else paste("Population", s)),
    Egg_NOS = apply(Egg_NOS, c(1, 2, 4), sum),
    Egg_HOS = apply(Egg_NOS, c(1, 2, 4), sum),
    Fry_NOS = apply(Fry_NOS, 1:3, sum),
    Fry_HOS = apply(Fry_HOS, 1:3, sum),
    Smolt_NOS = apply(Smolt_NOS, 1:3, sum),
    Smolt_HOS = apply(Smolt_HOS, 1:3, sum),
    Smolt_Rel = apply(Smolt_Rel, 1:3, sum),
    Njuv_NOS = apply(Njuv_NOS, 1:4, sum),
    Njuv_HOS = apply(Njuv_HOS, 1:4, sum),
    Return_NOS = apply(Return_NOS, 1:4, sum),
    Return_HOS = apply(Return_HOS, 1:4, sum),
    Escapement_NOS = apply(Escapement_NOS, 1:4, sum),
    Escapement_HOS = apply(Escapement_HOS, 1:4, sum),
    NOB = apply(NOB, c(1:2, 4), sum),
    HOB = apply(HOB, c(1:2, 4), sum),
    HOB_stray = apply(HOB_stray, c(1:2, 4), sum),
    HOB_import = apply(HOB_import, c(1:2, 4), sum),
    NOS = apply(NOS, 1:4, sum),
    HOS = apply(HOS, 1:4, sum),
    HOS_stray = apply(HOS_stray, 1:4, sum),
    HOS_effective = apply(HOS_effective, 1:4, sum),
    KPT_NOS = apply(KPT_NOS, c(1:2, 4), sum),
    KT_NOS = apply(KT_NOS, c(1:2, 4), sum),
    KPT_HOS = apply(KPT_HOS, c(1:2, 4), sum),
    KT_HOS = apply(KT_HOS, c(1:2, 4), sum),
    DPT_NOS = apply(DPT_NOS, c(1:2, 4), sum),
    DT_NOS = apply(DT_NOS, c(1:2, 4), sum),
    DPT_HOS = apply(DPT_HOS, c(1:2, 4), sum),
    DT_HOS = apply(DT_HOS, c(1:2, 4), sum),
    UPT_NOS = UPT_NOS_,
    UT_NOS = UT_NOS_,
    UPT_HOS = UPT_HOS_,
    UT_HOS = UT_HOS_,
    ExPT_NOS = ExPT_NOS_,
    ExT_NOS = ExT_NOS_,
    ExPT_HOS = ExPT_HOS_,
    ExT_HOS = ExT_HOS_,
    fitness = fitness,
    pNOB = pNOB,
    pHOS_census = pHOS_census,
    pHOS_effective = pHOS_effective,
    PNI = PNI,
    p_wild = p_wild,
    Mjuv_loss = apply(Mjuv_loss_NOS[, , , , 1, drop = FALSE], 1:4, identity)
  )

  if (n_r > 1) {
    SMSE@Misc$RS <- list(
      Smolt = Smolt_Rel, Esc = Escapement_HOS, HOS = HOS, Egg = apply(Egg_HOS, c(1, 2, 4, 5), sum)
    )
  }

  if (n_g > 1) {
    SMSE@Misc$LHG <- list(
      Fry = Fry_NOS + Fry_HOS,
      Smolt = Smolt_NOS + Smolt_HOS,
      Esc = Escapement_NOS,
      NOS = NOS,
      Egg = apply(Egg_NOS, c(1, 2, 4, 5), sum) # sum over ages
    )
  }

  return(SMSE)
}
