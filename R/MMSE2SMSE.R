
#' @rdname salmonMSE-int
#' @param MMSE Object of class \linkS4class{MMSE} returned from MSEtool
#' @param Harvest_MMP Optional harvest function created by [make_Harvest_MMP()]
#' @param N Data frame of abundance at age saved in the [salmonMSE_env] environment during the simulation
#' @param Ford Data frame of phenotypic trait values saved in the [salmonMSE_env] environment during the simulation
#' @param state Data frame of state variables saved in the [salmonMSE_env] environment during the simulation
#' @return
#' `MMSE2SMSE`: \linkS4class{SMSE} object
#' @export
MMSE2SMSE <- function(MMSE, SOM, Harvest_MMP, N, Ford, state) {
  ns <- 1 # Number of stocks
  nage <- SOM@maxage

  # Declare arrays
  Escapement_NOS <- Escapement_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  Return_NOS <- Return_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  NOB <- HOB <- array(0, c(SOM@nsim, ns, SOM@proyears))
  KPT_NOS <- KT_NOS <- KPT_HOS <- KT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  DPT_NOS <- DT_NOS <- DPT_HOS <- DT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  UPT_NOS <- UT_NOS <- UPT_HOS <- UT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  ExPT_NOS <- ExT_NOS <- ExPT_HOS <- ExT_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))

  Fry_NOS <- Fry_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  Smolt_NOS <- Smolt_HOS <- Smolt_Rel <- array(0, c(SOM@nsim, ns, SOM@proyears))

  NOS <- HOS <- HOS_effective <- array(0, c(SOM@nsim, ns, SOM@proyears))
  fitness <- array(NA_real_, c(SOM@nsim, ns, 2, SOM@proyears))

  Mocean_loss <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  PNI <- array(NA_real_, c(SOM@nsim, ns, SOM@proyears))
  p_wild <- array(NA_real_, c(SOM@nsim, ns, SOM@proyears))

  # NOS state variables from MMSE object
  p_NOS_imm <- 1
  p_NOS_return <- 2 # MSEtool population index for returning NOS
  p_NOS_escapement <- 3 # NOS escapement

  # openMSE year and age indices
  t1 <- seq(1, 2 * SOM@proyears, 2)
  t2 <- seq(2, 2 * SOM@proyears, 2)

  a1 <- seq(1, 2 * SOM@maxage + 1, 2)
  a2 <- seq(2, 2 * SOM@maxage + 1, 2)

  a_imm <- a1[-length(a1)]
  a_return <- a2
  a_esc <- a1[-1]
  #age_escapement <- age_bio + 1 # MSEtool age class

  f <- 1 # Fleet
  mp <- 1 # Only MP run per OM

  y_spawnOM <- which(MMSE@SSB[1, p_NOS_escapement, mp, ] > 0)
  y_spawn <- 0.5 * (y_spawnOM - 1)

  Return_NOS[, ns, , ] <- apply(MMSE@N[, p_NOS_return, a_return, mp, t2, ], 1:3, sum)
  Escapement_NOS[, ns, , -SOM@proyears] <- apply(MMSE@N[, p_NOS_escapement, a_esc, mp, t1, ], 1:3, sum)[, , - 1]

  # Kept catch
  KPT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_imm, f, mp, t1]
  KT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_return, f, mp, t2]

  # Total discards (live + dead)
  DPT_NOS[, ns, ] <- MMSE@Removals[, p_NOS_imm, f, mp, t1] - MMSE@Catch[, p_NOS_imm, f, mp, t1]
  DT_NOS[, ns, ] <- MMSE@Removals[, p_NOS_return, f, mp, t2] - MMSE@Catch[, p_NOS_return, f, mp, t2]

  # Harvest rate from kept catch
  UPT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_imm, f, mp, t1]/apply(MMSE@N[, p_NOS_imm, a_imm, mp, t1, ], c(1, 3), sum)
  UPT_NOS[is.na(UPT_NOS)] <- 0

  UT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_return, f, mp, t2]/apply(MMSE@N[, p_NOS_return, a_return, mp, t2, ], c(1, 3), sum)
  UT_NOS[is.na(UT_NOS)] <- 0

  # Exploitation rate from kept + dead discards (DD)
  DDPT_NOS <- SOM@release_mort[1] * DPT_NOS[, ns, ]
  ExPT_NOS[, ns, ] <- (KPT_NOS[, ns, ] + DDPT_NOS)/apply(MMSE@N[, p_NOS_imm, a_imm, mp, t1, ], c(1, 3), sum)
  ExPT_NOS[is.na(ExPT_NOS)] <- 0

  DDT_NOS <- SOM@release_mort[2] * DT_NOS[, ns, ]
  ExT_NOS[, ns, ] <- (KT_NOS[, ns, ] + DDT_NOS)/apply(MMSE@N[, p_NOS_return, a_return, mp, t2, ], c(1, 3), sum)
  ExT_NOS[is.na(ExT_NOS)] <- 0

  do_hatchery <- SOM@n_subyearling > 0 || SOM@n_yearling > 0
  if (do_hatchery) {
    # HOS state variables from MMSE object
    p_HOS_imm <- 4
    p_HOS_return <- 5 # Population index for returning HOS
    p_HOS_escapement <- 6 # NOS escapement

    Return_HOS[, ns, , ] <- apply(MMSE@N[, p_HOS_return, a_return, mp, t2, ], 1:3, sum)
    Escapement_HOS[, ns, , -SOM@proyears] <- apply(MMSE@N[, p_HOS_escapement, a_esc, mp, t1, ], 1:3, sum)[, , -1]

    # Kept catch
    KPT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_imm, f, mp, t1]
    KT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_return, f, mp, t2]

    # Total discards (live + dead)
    DPT_HOS[, ns, ] <- MMSE@Removals[, p_HOS_imm, f, mp, t1] - MMSE@Catch[, p_HOS_imm, f, mp, t1]
    DT_HOS[, ns, ] <- MMSE@Removals[, p_HOS_return, f, mp, t2] - MMSE@Catch[, p_HOS_return, f, mp, t2]

    # Harvest rate from kept catch
    UPT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_imm, f, mp, t1]/apply(MMSE@N[, p_HOS_imm, a_imm, mp, t1, ], c(1, 3), sum)
    UPT_HOS[is.na(UPT_HOS)] <- 0

    UT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_return, f, mp, t2]/apply(MMSE@N[, p_HOS_return, a_return, mp, t2, ], c(1, 3), sum)
    UT_HOS[is.na(UT_HOS)] <- 0

    # Exploitation rate from kept + dead discards (DD)
    DDPT_HOS <- SOM@release_mort[1] * DPT_HOS[, ns, ]
    ExPT_HOS[, ns, ] <- (KPT_HOS[, ns, ] + DDPT_HOS)/apply(MMSE@N[, p_HOS_imm, a_imm, mp, t1, ], c(1, 3), sum)
    ExPT_HOS[is.na(ExPT_HOS)] <- 0

    DDT_HOS <- SOM@release_mort[2] * DT_HOS[, ns, ]
    ExT_HOS[, ns, ] <- (KT_HOS[, ns, ] + DDT_HOS)/apply(MMSE@N[, p_HOS_return, a_return, mp, t2, ], c(1, 3), sum)
    ExT_HOS[is.na(ExT_HOS)] <- 0

    # NOS + HOS state variables from salmonMSE
    ngen <- length(unique(salmonMSE_env$N$t))
    if (length(y_spawn) != ngen) warning("Number of generations in salmonMSE state variables does not match generations in openMSE")

    NOS[, ns, y_spawn] <- get_salmonMSE_var(N, var = "NOS")
    Fry_NOS[, ns, y_spawn] <- get_salmonMSE_var(state, var = "fry_NOS")
    Smolt_NOS[, ns, y_spawn + 1] <- get_salmonMSE_var(state, var = "smolt_NOS")

    HOS[, ns, y_spawn] <- get_salmonMSE_var(N, var = "HOS")
    HOS_effective[, ns, y_spawn] <- get_salmonMSE_var(N, var = "HOS_effective")

    Fry_HOS[, ns, y_spawn] <- get_salmonMSE_var(state, var = "fry_HOS")
    Smolt_HOS[, ns, y_spawn + 1] <- get_salmonMSE_var(state, var = "smolt_HOS")

    # Broodtake & fitness
    NOB[, ns, y_spawn] <- get_salmonMSE_var(N, var = "NOB")
    HOB[, ns, y_spawn] <- get_salmonMSE_var(N, var = "HOB")

    fitness[, ns, 1, y_spawn + 1] <- get_salmonMSE_var(state, var = "fitness_natural")
    fitness[, ns, 2, y_spawn + 1] <- get_salmonMSE_var(state, var = "fitness_hatchery")

    # Smolt releases and SAR loss from openMSE
    p_smolt_rel <- 4
    a_smolt <- 1
    Smolt_Rel[, ns, y_spawn + 1] <- apply(MMSE@N[, p_smolt_rel, a_smolt, mp, y_spawnOM, ], 1:2, sum)

    p_smolt <- 1
    Mocean_loss[, ns, , ] <- MMSE@Misc$MICE$M_ageArray[, p_smolt, a2, mp, t2]

    pNOB <- get_salmonMSE_var(state, var = "pNOB")
    pHOSeff <- get_salmonMSE_var(state, var = "pHOSeff")
    pHOScensus <- get_salmonMSE_var(state, var = "pHOScensus")

    PNI[, ns, y_spawn] <- pNOB/(pNOB + pHOSeff) # Withler et al. 2018, page 17

    NOS_a <- HOScensus_a <- array(0, c(SOM@nsim, ns, SOM@maxage, SOM@proyears))
    NOS_a[, ns, , y_spawn] <- get_salmonMSE_agevar(N, "NOS")
    HOScensus_a[, ns, , y_spawn] <- get_salmonMSE_agevar(N, "HOS")

    p_wild[, ns, ] <- calc_pwild_age(NOS_a[, ns, , ], HOScensus_a[, ns, , ], SOM@fec, SOM@gamma)

  } else {
    # If no hatchery, the NOS escapement is also the NOS, fry_NOS is the spawning output
    NOS[] <- apply(Escapement_NOS, c(1, 2, 4), sum)

    p_smolt <- 1 # MSEtool population index for immature NOS
    p_spawn <- 3 # NOS escapement
    a_smolt <- 1

    Fry_NOS[, ns, y_spawn] <- MMSE@SSB[, p_spawn, mp, y_spawnOM] # -1 from 1-year lag
    Smolt_NOS[, ns, y_spawn + 1] <- apply(MMSE@N[, p_smolt, a_smolt, mp, y_spawnOM, ], 1:2, sum)

    PNI[, ns, y_spawn] <- p_wild[, ns, y_spawn] <- 1
  }

  SMSE <- new(
    "SMSE",
    Name = "salmon MSE results",
    nyears = SOM@nyears,
    proyears = SOM@proyears,
    nsim = SOM@nsim,
    nstocks = ns,
    Snames = "Single CU with hatchery",
    Fry_NOS = Fry_NOS,
    Fry_HOS = Fry_HOS,
    Smolt_NOS = Smolt_NOS,
    Smolt_HOS = Smolt_HOS,
    Smolt_Rel = Smolt_Rel,
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
    PNI = PNI,
    p_wild = p_wild,
    Mocean_loss = Mocean_loss
  )
  if (!missing(Harvest_MMP)) SMSE@Misc$Harvest_MMP <- Harvest_MMP
  if (!missing(Ford) && nrow(Ford)) SMSE@Misc$Ford <- salmonMSE_env$Ford

  return(SMSE)
}


#' @importFrom dplyr filter summarise
#' @importFrom reshape2 acast
get_salmonMSE_var <- function(d, var = "fry_NOS", p_smolt = 1) {
  x <- t <- NULL
  dplyr::filter(d, .data$p_smolt == .env$p_smolt) %>%
    summarise(value = sum(.data[[var]]), .by = c(x, t)) %>%
    reshape2::acast(list("x", "t"), value.var = "value")
}

get_salmonMSE_agevar <- function(d, var = "fry_NOS", p_smolt = 1) {
  dplyr::filter(d, .data$p_smolt == .env$p_smolt) %>%
    reshape2::acast(list("x", "a", "t"), value.var = var)
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

