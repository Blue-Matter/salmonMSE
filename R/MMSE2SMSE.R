
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
  fitness <- array(0, c(SOM@nsim, ns, SOM@proyears))

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

  # Total catch
  KPT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_imm, f, mp, t1]
  KT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_return, f, mp, t2]

  DPT_NOS[, ns, ] <- MMSE@Removals[, p_NOS_imm, f, mp, t1] - MMSE@Catch[, p_NOS_imm, f, mp, t1]
  DT_NOS[, ns, ] <- MMSE@Removals[, p_NOS_return, f, mp, t2] - MMSE@Catch[, p_NOS_return, f, mp, t2]

  ExPT_NOS[, ns, ] <- 1 - exp(-MMSE@FM[, p_NOS_imm, f, mp, t1])
  UPT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_imm, f, mp, t1]/apply(MMSE@N[, p_NOS_imm, a_imm, mp, t1, ], c(1, 3), sum)
  UPT_NOS[is.na(UPT_NOS)] <- 0

  ExT_NOS[, ns, ] <- 1 - exp(-MMSE@FM[, p_NOS_return, f, mp, t2])
  UT_NOS[, ns, ] <- MMSE@Catch[, p_NOS_return, f, mp, t2]/apply(MMSE@N[, p_NOS_return, a_return, mp, t2, ], c(1, 3), sum)
  UT_NOS[is.na(UT_NOS)] <- 0

  do_hatchery <- SOM@n_subyearling > 0 || SOM@n_yearling > 0
  if (do_hatchery) {
    # HOS state variables from MMSE object
    p_HOS_imm <- 4
    p_HOS_return <- 5 # Population index for returning HOS
    p_HOS_escapement <- 6 # NOS escapement

    Return_HOS[, ns, , ] <- apply(MMSE@N[, p_HOS_return, a_return, mp, t2, ], 1:3, sum)
    Escapement_HOS[, ns, , -SOM@proyears] <- apply(MMSE@N[, p_HOS_escapement, a_esc, mp, t1, ], 1:3, sum)[, , -1]

    KPT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_imm, f, mp, t1]
    KT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_return, f, mp, t2]

    DPT_HOS[, ns, ] <- MMSE@Removals[, p_HOS_imm, f, mp, t1] - MMSE@Catch[, p_HOS_imm, f, mp, t1]
    DT_HOS[, ns, ] <- MMSE@Removals[, p_HOS_return, f, mp, t2] - MMSE@Catch[, p_HOS_return, f, mp, t2]

    ExPT_HOS[, ns, ] <- 1 - exp(-MMSE@FM[, p_HOS_imm, f, mp, t1])
    UPT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_imm, f, mp, t1]/apply(MMSE@N[, p_HOS_imm, a_imm, mp, t1, ], c(1, 3), sum)
    UPT_HOS[is.na(UPT_HOS)] <- 0

    ExT_HOS[, ns, ] <- 1 - exp(-MMSE@FM[, p_HOS_return, f, mp, t1])
    UT_HOS[, ns, ] <- MMSE@Catch[, p_HOS_return, f, mp, t2]/apply(MMSE@N[, p_HOS_return, a_return, mp, t2, ], c(1, 3), sum)
    UT_HOS[is.na(UT_HOS)] <- 0

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

    fitness[, ns, y_spawn] <- get_salmonMSE_var(state, var = "fitness") # Lag by generation?

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
    p_wild[, ns, y_spawn] <- sapply(1:ncol(pHOScensus), function(t) {
      if (t == 1) {
        rep(NA_real_, nrow(pHOScensus))
      } else {
        calc_pwild(pHOScensus[, t], pHOScensus[, t-1], SOM@gamma)
      }
    })

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
get_salmonMSE_var <- function(N, var = "fry_NOS", p_smolt = 1) {
  x <- t <- NULL
  filter(N, .data$p_smolt == .env$p_smolt) %>%
    summarise(value = sum(.data[[var]]), .by = c(x, t)) %>%
    reshape2::acast(list("x", "t"), value.var = "value")
}


# Withler et al. 2018, page 27
calc_pwild <- function(pHOS_cur, pHOS_prev, gamma) {
  num <- (1 - pHOS_prev)^2
  denom <- num + 2 * gamma * pHOS_prev * (1 - pHOS_prev) + gamma * gamma * pHOS_prev^2
  (1 - pHOS_cur) * num/denom
}
