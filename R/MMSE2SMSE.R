
#' @rdname salmonMSE-int
#' @param MMSE Object of class [MSEtool::MMSE-class] returned from MSEtool
#' @param Harvest_MMP Optional harvest function created by [make_Harvest_MMP()]
#' @param N Data frame of natural origin abundance at age saved in the [salmonMSE_env] environment during the simulation
#' @param stateN Data frame of natural origin state variables saved in the [salmonMSE_env] environment during the simulation
#' @param Ford Data frame of phenotypic trait values saved in the [salmonMSE_env] environment during the simulation
#' @param H Data frame of hatchery origin abundance at age saved in the [salmonMSE_env] environment during the simulation
#' @param stateH Data frame of hatchery origin state variables saved in the [salmonMSE_env] environment during the simulation
#' @return
#' `MMSE2SMSE`: \linkS4class{SMSE} object
#' @export
MMSE2SMSE <- function(MMSE, SOM, Harvest_MMP, N, stateN, Ford, H, stateH) {
  ns <- length(SOM@Bio) # Number of stocks
  nage <- SOM@Bio[[1]]@maxage
  nyears_real <- 2

  # Declare arrays
  Njuv_NOS <- Njuv_HOS <- Escapement_NOS <- Escapement_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  Return_NOS <- Return_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  NOB <- HOB <- HOB_stray <- HOB_import <- array(0, c(SOM@nsim, ns, SOM@proyears))
  KPT_NOS <- KT_NOS <- KPT_HOS <- KT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  DPT_NOS <- DT_NOS <- DPT_HOS <- DT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  UPT_NOS <- UT_NOS <- UPT_HOS <- UT_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  ExPT_NOS <- ExT_NOS <- ExPT_HOS <- ExT_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))

  Egg_NOS <- Egg_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  Fry_NOS <- Fry_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  Smolt_NOS <- Smolt_HOS <- Smolt_Rel <- array(0, c(SOM@nsim, ns, SOM@proyears))

  NOS <- HOS <- HOS_stray <- HOS_effective <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  fitness <- array(NA_real_, c(SOM@nsim, ns, 2, SOM@proyears))

  pNOB <- pHOS_census <- pHOS_effective <- PNI <- p_wild <- array(NA_real_, c(SOM@nsim, ns, SOM@proyears))

  Mjuv_loss <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))

  LHG <- RS <- vector("list", length = ns)

  # openMSE year and age indices
  t1 <- seq(1, 2 * SOM@proyears, 2)
  t2 <- seq(2, 2 * SOM@proyears, 2)

  a1 <- seq(1, 2 * nage + 1, 2)
  a2 <- seq(2, 2 * nage + 1, 2)

  a_imm <- a1[-length(a1)]
  a_return <- a2
  a_esc <- a1[-1]

  f <- 1 # Fleet
  mp <- 1 # Only MP run per OM

  pindex <- make_stock_index(SOM)

  # Loop over s
  for (s in 1:ns) {

    # NOS state variables from MMSE object (length n_g)
    p_NOS_imm <- pindex$p[pindex$s == s & pindex$origin == "natural" & pindex$stage == "juvenile"]
    p_NOS_return <- pindex$p[pindex$s == s & pindex$origin == "natural" & pindex$stage == "recruitment"] # MSEtool population index for returning NOS
    p_NOS_escapement <- pindex$p[pindex$s == s & pindex$origin == "natural" & pindex$stage == "escapement"] # NOS escapement

    #y_spawnOM <- which(colSums(MMSE@SSB[, p_NOS_escapement[1], mp, ]) > 0) # Odd time steps, indicated by presence of escapement
    #if (y_spawnOM[1] == 1) {
    #  y_spawnOM <- y_spawnOM[-1] # This is actually the spawning from the last historical year
    #}
    #y_spawn <- 0.5 * (y_spawnOM - 1)
    y_spawnOM <- t1
    if (y_spawnOM[1] == 1) {
      y_spawnOM <- y_spawnOM[-1] # This is actually the spawning from the last historical year
    }
    y_spawn <- seq(1, length(unique(N$t)))

    Njuv_NOS[, s, , ] <- apply(MMSE@N[, p_NOS_imm, a_imm, mp, t1, , drop = FALSE], c(1, 3, 5), sum)
    Return_NOS[, s, , ] <- apply(MMSE@N[, p_NOS_return, a_return, mp, t2, , drop = FALSE], c(1, 3, 5), sum)
    Escapement_NOS[, s, , y_spawn] <- apply(MMSE@N[, p_NOS_escapement, a_esc, mp, y_spawnOM, , drop = FALSE], c(1, 3, 5), sum)

    # Kept catch
    KPT_NOS[, s, ] <- apply(MMSE@Catch[, p_NOS_imm, f, mp, t1, drop = FALSE], c(1, 5), sum)
    KT_NOS[, s, ] <- apply(MMSE@Catch[, p_NOS_return, f, mp, t2, drop = FALSE], c(1, 5), sum)

    # Total discards (live + dead)
    DPT_NOS[, s, ] <- apply(MMSE@Removals[, p_NOS_imm, f, mp, t1, drop = FALSE], c(1, 5), sum) - KPT_NOS[, s, ]
    DT_NOS[, s, ] <- apply(MMSE@Removals[, p_NOS_return, f, mp, t2, drop = FALSE], c(1, 5), sum) - KT_NOS[, s, ]

    # Exploitation rate from kept + dead discards (DD)
    vulPT <- array(SOM@Harvest[[s]]@vulPT, c(SOM@nsim, nage, length(t1)))
    DDPT_NOS <- SOM@Harvest[[s]]@release_mort[1] * DPT_NOS[, s, ]
    ExPT_NOS[, s, , ] <- local({
      Ka_NOS <- sapply(1:length(p_NOS_imm), function(i) { # sim x age x year
        FM <- MMSE@FM[, p_NOS_imm[i], f, mp, t1]
        Fa <- sapply(1:length(t1), function(t) vulPT[, , t] * FM[, t], simplify = "array")
        N <- apply(MMSE@N[, p_NOS_imm[i], a_imm, mp, t1, ], 1:3, sum)
        (1 - exp(-Fa)) * N
      }, simplify = "array") %>%
        apply(1:3, sum)
      Ka_NOS/Njuv_NOS[, s, , ]
    })
    ExPT_NOS[is.na(ExPT_NOS)] <- 0

    vulT <- array(SOM@Harvest[[s]]@vulT, c(SOM@nsim, nage, length(t2)))
    DDT_NOS <- SOM@Harvest[[s]]@release_mort[2] * DT_NOS[, s, ]
    ExT_NOS[, s, , ] <- local({
      Ka_NOS <- sapply(1:length(p_NOS_return), function(i) { # sim x age x year
        FM <- MMSE@FM[, p_NOS_return[i], f, mp, t2]
        Fa <- sapply(1:length(t2), function(t) vulT[, , t] * FM[, t], simplify = "array")
        N <- apply(MMSE@N[, p_NOS_return[i], a_return, mp, t2, ], 1:3, sum)
        (1 - exp(-Fa)) * N
      }, simplify = "array") %>%
        apply(1:3, sum)
      Ka_NOS/Return_NOS[, s, , ]
    })
    ExT_NOS[is.na(ExT_NOS)] <- 0

    # Harvest rate from kept catch
    if (all(!DDPT_NOS)) UPT_NOS[, s, , ] <- ExPT_NOS[, s, , ]
    if (all(!DDT_NOS)) UT_NOS[, s, , ] <- ExT_NOS[, s, , ]

    do_hatchery <- sum(SOM@Hatchery[[s]]@n_subyearling > 0, SOM@Hatchery[[s]]@n_yearling) > 0
    has_strays <- any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)
    n_g <- length(unique(pindex$g[pindex$s == s & pindex$origin == "natural"]))
    n_r <- length(unique(pindex$r[pindex$s == s & pindex$origin == "hatchery"]))

    s_enroute <- SOM@Bio[[s]]@s_enroute
    use_smolt_func <- TRUE

    if (use_smolt_func) {
      # Sum across LHG
      NOS[, s, , y_spawn] <- get_salmonMSE_agevar(N, var = "NOS", s)
      Egg_NOS[, s, y_spawn] <- get_salmonMSE_var(stateN, var = "Egg_NOS", s)
      Fry_NOS[, s, y_spawn + 1] <- get_salmonMSE_var(stateN, var = "Fry_NOS", s)
      Smolt_NOS[, s, y_spawn + 1] <- get_salmonMSE_var(stateN, var = "smolt_NOS", s)

    } else {
      # If no hatchery, the NOS escapement is also the NOS, Egg_NOS is the spawning output
      NOS[, s, , ] <- Escapement_NOS[, s, , ]
      Egg_NOS[, s, y_spawn] <- apply(MMSE@SSB[, p_NOS_escapement, mp, y_spawnOM, drop = FALSE], c(1, 4), sum) # -1 from 1-year lag
      Fry_NOS[, s, seq(2, SOM@proyears)] <- Egg_NOS[, s, seq(2, SOM@proyears) - 1]

      a_smolt <- 1
      Smolt_NOS[, s, ] <- Njuv_NOS[, s, a_smolt, ]
    }

    if (do_hatchery || has_strays) {

      # HOS state variables from MMSE object
      p_HOS_imm <- pindex$p[pindex$s == s & pindex$origin == "hatchery" & pindex$stage == "juvenile"]
      p_HOS_return <- pindex$p[pindex$s == s & pindex$origin == "hatchery" & pindex$stage == "recruitment"] # MSEtool population index for returning HOS
      p_HOS_escapement <- pindex$p[pindex$s == s & pindex$origin == "hatchery" & pindex$stage == "escapement"] # HOS escapement

      Njuv_HOS[, s, , ] <- apply(MMSE@N[, p_HOS_imm, a_imm, mp, t1, , drop = FALSE], c(1, 3, 5), sum)
      Return_HOS[, s, , ] <- apply(MMSE@N[, p_HOS_return, a_return, mp, t2, , drop = FALSE], c(1, 3, 5), sum)
      Escapement_HOS[, s, , y_spawn] <- apply(MMSE@N[, p_HOS_escapement, a_esc, mp, y_spawnOM, , drop = FALSE], c(1, 3, 5), sum)

      # Kept catch
      KPT_HOS[, s, ] <- apply(MMSE@Catch[, p_HOS_imm, f, mp, t1, drop = FALSE], c(1, 5), sum)
      KT_HOS[, s, ] <- apply(MMSE@Catch[, p_HOS_return, f, mp, t2, drop = FALSE], c(1, 5), sum)

      # Total discards (live + dead)
      DPT_HOS[, s, ] <- apply(MMSE@Removals[, p_HOS_imm, f, mp, t1, drop = FALSE], c(1, 5), sum) - KPT_HOS[, s, ]
      DT_HOS[, s, ] <- apply(MMSE@Removals[, p_HOS_return, f, mp, t2, drop = FALSE], c(1, 5), sum) - KT_HOS[, s, ]

      # Exploitation rate from kept + dead discards (DD)
      DDPT_HOS <- SOM@Harvest[[s]]@release_mort[1] * DPT_HOS[, s, ]
      ExPT_HOS[, s, , ] <- local({
        Ka_HOS <- sapply(1:length(p_HOS_imm), function(i) { # sim x age x year
          FM <- MMSE@FM[, p_HOS_imm[i], f, mp, t1]
          Fa <- sapply(1:length(t1), function(t) vulPT[, , t] * FM[, t], simplify = "array")
          N <- apply(MMSE@N[, p_HOS_imm[i], a_imm, mp, t1, ], 1:3, sum)
          (1 - exp(-Fa)) * N
        }, simplify = "array") %>%
          apply(1:3, sum)
        Ka_HOS/Njuv_HOS[, s, , ]
      })
      ExPT_HOS[is.na(ExPT_HOS)] <- 0

      DDT_HOS <- SOM@Harvest[[s]]@release_mort[2] * DT_HOS[, s, ]
      ExT_HOS[, s, , ] <- local({
        Ka_HOS <- sapply(1:length(p_HOS_return), function(i) { # sim x age x year
          FM <- MMSE@FM[, p_HOS_return[i], f, mp, t2]
          Fa <- sapply(1:length(t2), function(t) vulT[, , t] * FM[, t], simplify = "array")
          N <- apply(MMSE@N[, p_HOS_return[i], a_return, mp, t2, ], 1:3, sum)
          (1 - exp(-Fa)) * N
        }, simplify = "array") %>%
          apply(1:3, sum)
        Ka_HOS/Return_HOS[, s, , ]
      })
      ExT_HOS[is.na(ExT_HOS)] <- 0

      # Harvest rate from kept catch
      if (all(!DDPT_HOS)) {
        UPT_HOS[, s, , ] <- ExPT_HOS[, s, , ]
      } else {
        stop("Dead discards of hatchery fish found")
      }

      if (all(!DDT_HOS)) {
        UT_HOS[, s, , ] <- ExT_HOS[, s, , ]
      } else {
        stop("Dead discards of hatchery fish found")
      }

      # NOS + HOS state variables from salmonMSE
      ngen <- length(unique(salmonMSE_env$N$t))
      if (length(y_spawn) != ngen) warning("Number of generations in salmonMSE state variables does not match generations in openMSE")

      # Sum across RS
      HOS[, s, , y_spawn] <- get_salmonMSE_agevar(H, var = "HOS", s)
      HOS_stray[, s, , y_spawn] <- get_salmonMSE_agevar(H, var = "HOS_stray", s)
      HOS_effective[, s, , y_spawn] <- get_salmonMSE_agevar(H, var = "HOS_effective", s)

      Egg_HOS[, s, y_spawn] <- get_salmonMSE_var(stateH, var = "Egg_HOS", s)
      Fry_HOS[, s, y_spawn + 1] <- get_salmonMSE_var(stateN, var = "Fry_HOS", s)
      Smolt_HOS[, ns, y_spawn + 1] <- get_salmonMSE_var(stateN, var = "smolt_HOS", s)

      # Broodtake & fitness
      NOB[, s, y_spawn] <- get_salmonMSE_var(N, var = "NOB", s)
      HOB[, s, y_spawn] <- get_salmonMSE_var(H, var = "HOB", s)
      HOB_stray[, s, y_spawn] <- get_salmonMSE_var(H, var = "HOB_stray", s)
      HOB_import[, s, y_spawn] <- get_salmonMSE_var(H, var = "HOB_import", s)

      fitness[, s, 1, y_spawn + 1] <- filter(Ford, .data$type == "natural", .data$t > 2 * nyears_real) %>%
        get_salmonMSE_var(var = "fitness", s)
      fitness[, s, 2, y_spawn + 1] <- filter(Ford, .data$type == "hatchery", .data$t > 2 * nyears_real) %>%
        get_salmonMSE_var(var = "fitness", s)

      # Convert fitness to NA for empty cohorts
      y_empty <- colSums(Smolt_NOS[, s, ] + Smolt_Rel[, s, ]) == 0
      fitness[, s, , y_empty] <- NA

      # Smolt releases and SAR loss from openMSE
      a_smolt <- 1

      smolt_rel_openmse <- apply(MMSE@N[, p_HOS_imm, a_smolt, mp, y_spawnOM, , drop = FALSE], c(1, 5), sum)
      #smolt_rel_salmonmse <- get_salmonMSE_var(stateH, var = "smolt_rel", s) # debugging purposes
      Smolt_Rel[, s, y_spawn + 1] <- smolt_rel_openmse

      if (!is.null(MMSE@Misc$MICE$M_ageArray)) {
        Mjuv_loss[, s, , ] <- MMSE@Misc$MICE$M_ageArray[, p_NOS_imm[1], a2, mp, t2] # Report first LHG only
      }

      pNOB[, s, y_spawn] <- get_salmonMSE_var(stateN, var = "pNOB", s, FUN = unique)
      pHOS_effective[, s, y_spawn] <- get_salmonMSE_var(stateN, var = "pHOSeff", s, FUN = unique)
      pHOS_census[, s, y_spawn] <- get_salmonMSE_var(stateN, var = "pHOScensus", s, FUN = unique)

      # Withler et al. 2018, page 17, 21
      h2 <- SOM@Hatchery[[s]]@heritability
      fitness_variance <- SOM@Hatchery[[s]]@fitness_variance

      # Exact PNI, see HSRG 2009, Appendix C, Eq. 33, Columbia River Hatchery Reform System-Wide Report
      PNI[, s, y_spawn] <- ifelse(
        pNOB[, s, y_spawn] > 0,
        pNOB[, s, y_spawn]/(pNOB[, s, y_spawn] + pHOS_effective[, s, y_spawn]),
        h2/(h2 + (1 - h2 + fitness_variance) * pHOS_effective[, s, y_spawn])
      )

      #PNI[, s, ] <- (h2 + (1 - h2 + fitness_variance) * pNOB[, s, ])/(h2 + (1 - h2 + fitness_variance) * (pHOS_effective[, s, ] + pNOB[, s, ]))

      p_wild[, s, ] <- calc_pwild_age(NOS[, s, , ], HOS[, s, , ], SOM@Bio[[s]]@fec[, , seq(1, SOM@proyears)], SOM@Hatchery[[s]]@gamma)

      if (n_r > 1) {
        Smolt_r <- array(NA_real_, c(SOM@nsim, n_r, SOM@proyears))
        Esc_r <- HOS_r <- array(NA_real_, c(SOM@nsim, n_r, nage, SOM@proyears))

        x <- a <- r <- t <- NULL

        a_smolt <- 1
        Smolt_r[] <- apply(MMSE@N[, p_HOS_imm, a_smolt, mp, t1, , drop = FALSE], c(1, 2, 5), sum)
        Esc_r[, , , y_spawn] <- apply(MMSE@N[, p_HOS_escapement, a_esc, mp, y_spawnOM, , drop = FALSE], c(1, 2, 3, 5), sum)
        HOS_r[, , , y_spawn] <- dplyr::filter(H, .data$s == .env$s) %>%
          summarise(value = sum(.data$HOS), .by = c(x, r, a, t)) %>%
          reshape2::acast(list("x", "r", "a", "t"), value.var = "value")

        RS[[s]] <- list(
          Smolt = Smolt_r,
          Esc = Esc_r,
          HOS = HOS_r
        )
      }

    } else {
      PNI[, s, y_spawn] <- p_wild[, s, y_spawn] <- 1
      pHOS_census[, s, y_spawn] <- pHOS_effective[, s, y_spawn] <- 0
    }

    if (n_g > 1) {
      Egg_g <- Fry_g <- Smolt_g <- array(NA_real_, c(SOM@nsim, n_g, SOM@proyears))
      Esc_g <- NOS_g <- array(NA_real_, c(SOM@nsim, n_g, nage, SOM@proyears))

      x <- a <- g <- t <- NULL
      Egg_g[, , y_spawn] <- dplyr::filter(stateN, .data$s == .env$s) %>%
        summarise(value = sum(.data$Egg_NOS), .by = c(x, g, t)) %>%
        reshape2::acast(list("x", "g", "t"), value.var = "value")

      Fry_g[, , y_spawn + 1] <- dplyr::filter(stateN, .data$s == .env$s) %>%
        summarise(value = sum(.data$Fry_NOS), .by = c(x, g, t)) %>%
        reshape2::acast(list("x", "g", "t"), value.var = "value")

      a_smolt <- 1
      Smolt_g[] <- apply(MMSE@N[, p_NOS_imm, a_smolt, mp, t1, , drop = FALSE], c(1, 2, 5), sum)
      Esc_g[, , , y_spawn] <- apply(MMSE@N[, p_NOS_escapement, a_esc, mp, y_spawnOM, , drop = FALSE], c(1, 2, 3, 5), sum)
      NOS_g[, , , y_spawn] <- dplyr::filter(N, .data$s == .env$s) %>%
        summarise(value = sum(.data$NOS), .by = c(x, g, a, t)) %>%
        reshape2::acast(list("x", "g", "a", "t"), value.var = "value")

      LHG[[s]] <- list(
        Egg = Egg_g,
        Fry = Fry_g,
        Smolt = Smolt_g,
        Esc = Esc_g,
        NOS = NOS_g
      )
    }
  }

  SMSE <- new(
    "SMSE",
    Name = "salmon MSE results",
    #nyears = SOM@nyears,
    proyears = SOM@proyears,
    nsim = SOM@nsim,
    nstocks = ns,
    Snames = sapply(1:ns, function(s) if (length(SOM@Bio[[s]]@Name)) SOM@Bio[[s]]@Name else paste("Population", s)),
    Egg_NOS = Egg_NOS,
    Egg_HOS = Egg_HOS,
    Fry_NOS = Fry_NOS,
    Fry_HOS = Fry_HOS,
    Smolt_NOS = Smolt_NOS,
    Smolt_HOS = Smolt_HOS,
    Smolt_Rel = Smolt_Rel,
    Njuv_NOS = Njuv_NOS,
    Njuv_HOS = Njuv_HOS,
    Return_NOS = Return_NOS,
    Return_HOS = Return_HOS,
    Escapement_NOS = Escapement_NOS,
    Escapement_HOS = Escapement_HOS,
    NOB = NOB,
    HOB = HOB,
    HOB = HOB_stray,
    HOB_import = HOB_import,
    NOS = NOS,
    HOS = HOS,
    HOS_stray = HOS_stray,
    HOS_effective = HOS_effective,
    KPT_NOS = KPT_NOS,
    KT_NOS = KT_NOS,
    KPT_HOS = KPT_HOS,
    KT_HOS = KT_HOS,
    DPT_NOS = DPT_NOS,
    DT_NOS = DT_NOS,
    DPT_HOS = DPT_HOS,
    DT_HOS = DT_HOS,
    UPT_NOS = UPT_NOS,
    UT_NOS = UT_NOS,
    UPT_HOS = UPT_HOS,
    UT_HOS = UT_HOS,
    ExPT_NOS = ExPT_NOS,
    ExT_NOS = ExT_NOS,
    ExPT_HOS = ExPT_HOS,
    ExT_HOS = ExT_HOS,
    fitness = fitness,
    pNOB = pNOB,
    pHOS_census = pHOS_census,
    pHOS_effective = pHOS_effective,
    PNI = PNI,
    p_wild = p_wild,
    Mjuv_loss = Mjuv_loss
  )
  if (!missing(Harvest_MMP)) SMSE@Misc$Harvest_MMP <- Harvest_MMP
  if (!missing(Ford) && nrow(Ford)) SMSE@Misc$Ford <- salmonMSE_env$Ford

  SMSE@Misc$SOM <- SOM
  SMSE@Misc$Ref <- calc_ref(SOM, check = FALSE)
  SMSE@Misc$LHG <- LHG
  SMSE@Misc$RS <- RS

  return(SMSE)
}


#' @importFrom dplyr filter summarise
#' @importFrom reshape2 acast
get_salmonMSE_var <- function(d, var = "Egg_NOS", s = 1, FUN = function(x) sum(x, na.rm = TRUE)) {
  x <- t <- NULL
  dplyr::filter(d, .data$s %in% .env$s) %>%
    summarise(value = FUN(.data[[var]]), .by = c(x, t)) %>%
    reshape2::acast(list("x", "t"), value.var = "value")
}

get_salmonMSE_agevar <- function(d, var = "Egg_NOS", s = 1, FUN = function(x) sum(x, na.rm = TRUE)) {
  x <- a <- t <- NULL
  dplyr::filter(d, .data$s %in% .env$s) %>%
    summarise(value = FUN(.data[[var]]), .by = c(x, a, t)) %>%
    reshape2::acast(list("x", "a", "t"), value.var = "value")
}


#' Proportion wild spawners
#'
#' @description Calculate the proportion of wild spawners from a time series of spawners
#' - `calc_pwild()` is the simple calculation based on the proportion of hatchery spawners
#' - `calc_pwild_age()` performs the calculation weighted by age class fecundity
#'
#' @param NOS_a Array `[nsim, maxage, years]` for natural spawners
#' @param HOS_a Array `[nsim, maxage, years]` for hatchery spawners
#' @param fec Array `[nsim, maxage, years]` for age class fecundity
#' @param gamma Numeric, reduced reproductive success of hatchery spawners
#' @param pHOS_cur Numeric, proportion of hatchery spawners in current generation
#' @param pHOS_prev Numeric, proportion of hatchery spawners in previous generation
#' @return `calc_pwild_age()` a matrix of pWILD by simulation and year. `calc_pwild()` returns a numeric
#' @keywords internal
calc_pwild_age <- function(NOS_a, HOS_a, fec, gamma) {

  nsim <- dim(NOS_a)[1]
  maxage <- dim(NOS_a)[2]
  proyears <- dim(NOS_a)[3]

  prob <- array(NA_real_, c(nsim, proyears))

  for (y in 1:proyears) {
    if (sum(NOS_a[, , y], HOS_a[, , y])) {
      pNOS_y <- matrix(1, nsim, maxage)
      pnatural_b <- phetero_b <- phatch_b <- matrix(0, nsim, maxage)

      pNOS_y[] <- NOS_a[, , y]/rowSums(NOS_a[, , y] + HOS_a[, , y])

      for (a in 1:maxage) {
        pNOS_b <- matrix(1, nsim, maxage)
        pHOS_b <- matrix(0, nsim, maxage)

        b <- y - a # brood_year
        if (b > 0) {
          pNOS_b[] <- NOS_a[, , b]/rowSums(fec[, , b] * NOS_a[, , b] + fec[, , b] * HOS_a[, , b])
          pHOS_b[] <- HOS_a[, , b]/rowSums(fec[, , b] * NOS_a[, , b] + fec[, , b] * HOS_a[, , b])
        }
        pnatural_b[, a] <- rowSums(pNOS_b, na.rm = TRUE)^2
        phetero_b[, a] <- 2 * gamma * rowSums(pNOS_b, na.rm = TRUE) * rowSums(pHOS_b, na.rm = TRUE)
        phatch_b[, a] <- gamma^2 * rowSums(pHOS_b, na.rm = TRUE)^2
      }
      denom <- pnatural_b + phetero_b + phatch_b

      prob[, y] <- rowSums(pNOS_y * pnatural_b/denom, na.rm = TRUE)
    }
  }

  return(prob)
}

