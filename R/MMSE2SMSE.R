
#' @rdname salmonMSE-int
#' @param MMSE Object of class \linkS4class{MMSE} returned from MSEtool
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

  # Declare arrays
  Njuv_NOS <- Njuv_HOS <- Escapement_NOS <- Escapement_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  Return_NOS <- Return_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  NOB <- HOB <- array(0, c(SOM@nsim, ns, SOM@proyears))
  KPT_NOS <- KT_NOS <- KPT_HOS <- KT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  DPT_NOS <- DT_NOS <- DPT_HOS <- DT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  UPT_NOS <- UT_NOS <- UPT_HOS <- UT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  ExPT_NOS <- ExT_NOS <- ExPT_HOS <- ExT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))

  Egg_NOS <- Egg_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  Smolt_NOS <- Smolt_HOS <- Smolt_Rel <- array(0, c(SOM@nsim, ns, SOM@proyears))

  NOS <- HOS <- HOS_effective <- array(0, c(SOM@nsim, ns, SOM@proyears))
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

    y_spawnOM <- which(MMSE@SSB[1, p_NOS_escapement[1], mp, ] > 0) # Odd time steps, indicated by presence of escapement
    if (y_spawnOM[1] == 1) {
      y_spawnOM <- y_spawnOM[-1] # This is actually the spawning from the last historical year
    }
    y_spawn <- 0.5 * (y_spawnOM - 1)

    Njuv_NOS[, s, , ] <- apply(MMSE@N[, p_NOS_imm, a_imm, mp, t1, , drop = FALSE], c(1, 3, 5), sum)
    Return_NOS[, s, , ] <- apply(MMSE@N[, p_NOS_return, a_return, mp, t2, , drop = FALSE], c(1, 3, 5), sum)
    Escapement_NOS[, s, , y_spawn] <- apply(MMSE@N[, p_NOS_escapement, a_esc, mp, y_spawnOM, , drop = FALSE], c(1, 3, 5), sum)

    # Kept catch
    KPT_NOS[, s, ] <- apply(MMSE@Catch[, p_NOS_imm, f, mp, t1, drop = FALSE], c(1, 5), sum)
    KT_NOS[, s, ] <- apply(MMSE@Catch[, p_NOS_return, f, mp, t2, drop = FALSE], c(1, 5), sum)

    # Total discards (live + dead)
    DPT_NOS[, s, ] <- apply(MMSE@Removals[, p_NOS_imm, f, mp, t1, drop = FALSE], c(1, 5), sum) - KPT_NOS[, s, ]
    DT_NOS[, s, ] <- apply(MMSE@Removals[, p_NOS_return, f, mp, t2, drop = FALSE], c(1, 5), sum) - KT_NOS[, s, ]

    # Harvest rate from kept catch
    vulPT <- array(SOM@Harvest[[s]]@vulPT, c(SOM@nsim, nage, length(t1)))
    vulNOS_imm <- apply(vulPT * Njuv_NOS[, s, , ], c(1, 3), sum)
    UPT_NOS[, s, ] <- KPT_NOS[, s, ]/vulNOS_imm
    UPT_NOS[is.na(UPT_NOS)] <- 0

    vulT <- array(SOM@Harvest[[s]]@vulT, c(SOM@nsim, nage, length(t2)))
    vulNOS_ret <- apply(vulT * Return_NOS[, s, , ], c(1, 3), sum)
    UT_NOS[, s, ] <- KT_NOS[, s, ]/vulNOS_ret
    UT_NOS[is.na(UT_NOS)] <- 0

    # Exploitation rate from kept + dead discards (DD)
    DDPT_NOS <- SOM@Harvest[[s]]@release_mort[1] * DPT_NOS[, s, ]
    ExPT_NOS[, s, ] <- (KPT_NOS[, s, ] + DDPT_NOS)/vulNOS_imm
    ExPT_NOS[is.na(ExPT_NOS)] <- 0

    DDT_NOS <- SOM@Harvest[[s]]@release_mort[2] * DT_NOS[, s, ]
    ExT_NOS[, s, ] <- (KT_NOS[, s, ] + DDT_NOS)/vulNOS_ret
    ExT_NOS[is.na(ExT_NOS)] <- 0

    do_hatchery <- sum(SOM@Hatchery[[s]]@n_subyearling > 0, SOM@Hatchery[[s]]@n_yearling) > 0
    has_strays <- any(SOM@stray[-s, s] > 0)
    n_g <- length(unique(pindex$g[pindex$s == s & pindex$origin == "natural"]))
    n_r <- length(unique(pindex$r[pindex$s == s & pindex$origin == "hatchery"]))

    s_enroute <- SOM@Bio[[s]]@s_enroute

    if (do_hatchery || has_strays || s_enroute < 1) {

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

      # Harvest rate from kept catch
      vulHOS_imm <- apply(vulPT * Njuv_HOS[, s, , ], c(1, 3), sum)
      UPT_HOS[, s, ] <- KPT_HOS[, s, ]/vulHOS_imm
      UPT_HOS[is.na(UPT_HOS)] <- 0

      vulHOS_ret <- apply(vulT * Return_HOS[, s, , ], c(1, 3), sum)
      UT_HOS[, s, ] <- KT_HOS[, s, ]/vulHOS_ret
      UT_HOS[is.na(UT_HOS)] <- 0

      # Exploitation rate from kept + dead discards (DD)
      DDPT_HOS <- SOM@Harvest[[s]]@release_mort[1] * DPT_HOS[, s, ]
      ExPT_HOS[, s, ] <- (KPT_HOS[, s, ] + DDPT_HOS)/vulHOS_imm
      ExPT_HOS[is.na(ExPT_HOS)] <- 0

      DDT_HOS <- SOM@Harvest[[s]]@release_mort[2] * DT_HOS[, s, ]
      ExT_HOS[, s, ] <- (KT_HOS[, s, ] + DDT_HOS)/vulHOS_ret
      ExT_HOS[is.na(ExT_HOS)] <- 0

      # NOS + HOS state variables from salmonMSE
      ngen <- length(unique(salmonMSE_env$N$t))
      if (length(y_spawn) != ngen) warning("Number of generations in salmonMSE state variables does not match generations in openMSE")

      # Sum across LHG
      NOS[, s, y_spawn] <- get_salmonMSE_var(N, var = "NOS", s)
      Egg_NOS[, s, y_spawn] <- get_salmonMSE_var(stateN, var = "Egg_NOS", s)
      Smolt_NOS[, s, y_spawn + 1] <- get_salmonMSE_var(stateN, var = "smolt_NOS", s)

      # Sum across RS
      HOS[, s, y_spawn] <- get_salmonMSE_var(H, var = "HOS", s)
      HOS_effective[, s, y_spawn] <- get_salmonMSE_var(H, var = "HOS_effective", s)

      Egg_HOS[, s, y_spawn] <- get_salmonMSE_var(stateH, var = "Egg_HOS", s)
      Smolt_HOS[, ns, y_spawn + 1] <- get_salmonMSE_var(stateN, var = "smolt_HOS", s)

      # Broodtake & fitness
      NOB[, s, y_spawn] <- get_salmonMSE_var(N, var = "NOB", s)
      HOB[, s, y_spawn] <- get_salmonMSE_var(H, var = "HOB", s)

      fitness[, s, 1, y_spawn + 1] <- filter(Ford, type == "natural") %>%
        get_salmonMSE_var(var = "fitness", s)
      fitness[, s, 2, y_spawn + 1] <- filter(Ford, type == "hatchery") %>%
        get_salmonMSE_var(var = "fitness", s)

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

      PNI[, s, y_spawn] <- pNOB[, s, y_spawn]/(pNOB[, s, y_spawn] + pHOS_effective[, s, y_spawn]) # Withler et al. 2018, page 17

      NOS_a <- HOScensus_a <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
      NOS_a[, s, , y_spawn] <- get_salmonMSE_agevar(N, "NOS", s)
      HOScensus_a[, s, , y_spawn] <- get_salmonMSE_agevar(H, "HOS", s)

      p_wild[, s, ] <- calc_pwild_age(NOS_a[, s, , ], HOScensus_a[, s, , ], SOM@Bio[[s]]@fec, SOM@Hatchery[[s]]@gamma)

    } else {
      # If no hatchery, the NOS escapement is also the NOS, Egg_NOS is the spawning output
      NOS[, s, ] <- apply(Escapement_NOS[, s, , ], c(1, 3), sum)
      Egg_NOS[, s, y_spawn] <- apply(MMSE@SSB[, p_NOS_escapement, mp, y_spawnOM, drop = FALSE], c(1, 4), sum) # -1 from 1-year lag

      a_smolt <- 1
      Smolt_NOS[, s, ] <- Njuv_NOS[, s, a_smolt, ]

      PNI[, s, y_spawn] <- p_wild[, s, y_spawn] <- 1
      pHOS_census[, s, y_spawn] <- pHOS_effective[, s, y_spawn] <- 0
    }

    if (n_g > 1) {
      Egg_g <- Smolt_g <- array(NA_real_, c(SOM@nsim, n_g, SOM@proyears))
      Esc_g <- NOS_g <- array(NA_real_, c(SOM@nsim, n_g, nage, SOM@proyears))

      x <- a <- g <- t <- NULL
      Egg_g[, , y_spawn] <- dplyr::filter(stateN, .data$s == .env$s) %>%
        summarise(value = sum(.data$Egg_NOS_g), .by = c(x, g, t)) %>%
        reshape2::acast(list("x", "g", "t"), value.var = "value")

      a_smolt <- 1
      Smolt_g[] <- apply(MMSE@N[, p_NOS_imm, a_smolt, mp, t1, , drop = FALSE], c(1, 2, 5), sum)
      Esc_g[, , , y_spawn] <- apply(MMSE@N[, p_NOS_escapement, a_esc, mp, y_spawnOM, , drop = FALSE], c(1, 2, 3, 5), sum)
      NOS_g[, , , y_spawn] <- dplyr::filter(N, .data$s == .env$s) %>%
        summarise(value = sum(.data$NOS), .by = c(x, g, a, t)) %>%
        reshape2::acast(list("x", "g", "a", "t"), value.var = "value")

      LHG[[s]] <- list(
        Egg = Egg_g,
        Smolt = Smolt_g,
        Esc = Esc_g,
        NOS = NOS_g
      )
    }

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

  }

  SMSE <- new(
    "SMSE",
    Name = "salmon MSE results",
    nyears = SOM@nyears,
    proyears = SOM@proyears,
    nsim = SOM@nsim,
    nstocks = ns,
    Snames = sapply(1:ns, function(s) if (length(SOM@Bio[[s]]@Name)) SOM@Bio[[s]]@Name else paste("Population", s)),
    Egg_NOS = Egg_NOS,
    Egg_HOS = Egg_HOS,
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
    NOS = NOS,
    HOS = HOS,
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
get_salmonMSE_var <- function(d, var = "Egg_NOS", s = 1, FUN = sum) {
  x <- t <- NULL
  dplyr::filter(d, .data$s %in% .env$s) %>%
    summarise(value = FUN(.data[[var]]), .by = c(x, t)) %>%
    reshape2::acast(list("x", "t"), value.var = "value")
}

get_salmonMSE_agevar <- function(d, var = "Egg_NOS", s = 1, FUN = sum) {
  x <- a <- t <- NULL
  dplyr::filter(d, .data$s %in% .env$s) %>%
    summarise(value = FUN(.data[[var]]), .by = c(x, a, t)) %>%
    reshape2::acast(list("x", "a", "t"), value.var = "value")
}

calc_pwild_age <- function(NOS_a, HOS_a, fec, gamma) {

  nsim <- dim(NOS_a)[1]
  maxage <- dim(NOS_a)[2]
  proyears <- dim(NOS_a)[3]

  fec <- matrix(fec, nsim, maxage, byrow = TRUE)

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
          pNOS_b[] <- NOS_a[, , b]/rowSums(fec * NOS_a[, , b] + fec * HOS_a[, , b])
          pHOS_b[] <- HOS_a[, , b]/rowSums(fec * NOS_a[, , b] + fec * HOS_a[, , b])
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

