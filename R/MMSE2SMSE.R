
#' @rdname salmonMSE-int
#' @param MMSE Object of class \linkS4class{MMSE} returned from MSEtool
#' @param Harvest_MMP Optional harvest function created by [make_Harvest_MMP()]
#' @param N Data frame of state variables saved in the [salmonMSE_env] environment during the simulation
#' @param Ford Data frame of phenotypic trait values saved in the [salmonMSE_env] environment during the simulation
#' @return
#' `MMSE2SMSE`: \linkS4class{SMSE} object
#' @export
MMSE2SMSE <- function(MMSE, SOM, Harvest_MMP, N, Ford) {

  ns <- 1 # Number of stocks
  nage <- SOM@maxage

  # Declare arrays
  Escapement_NOS <- Escapement_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  Return_NOS <- Return_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@proyears))
  NOB <- HOB <- array(0, c(SOM@nsim, ns, SOM@proyears))
  Catch_NOS <- Catch_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  U_NOS <- U_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))

  Fry_NOS <- Fry_HOS <- array(0, c(SOM@nsim, ns, SOM@proyears))
  Smolt_NOS <- Smolt_HOS <- Smolt_Rel <- array(0, c(SOM@nsim, ns, SOM@proyears))

  NOS <- HOS <- HOS_effective <- array(0, c(SOM@nsim, ns, SOM@proyears))
  fitness <- array(0, c(SOM@nsim, ns, SOM@proyears))

  SAR_loss <- array(0, c(SOM@nsim, ns, SOM@proyears))
  PNI <- array(NA_real_, c(SOM@nsim, ns, SOM@proyears))
  p_wild <- array(NA_real_, c(SOM@nsim, ns, SOM@proyears))

  # NOS state variables from MMSE object
  p_NOS_return <- 2 # MSEtool population index for returning NOS
  p_NOS_escapement <- 3 # NOS escapement
  age_escapement <- 4 # MSEtool age class

  Return_NOS[, ns, nage, ] <- apply(MMSE@N[, p_NOS_return, nage, 1, , ], 1:2, sum)

  y_spawn <- Return_NOS[1, 1, nage, ] > 0 # Years for which there is spawning

  Escapement_NOS[, ns, nage, -SOM@proyears] <- apply(MMSE@N[, p_NOS_escapement, age_escapement, 1, -1, ], 1:2, sum)
  Catch_NOS[, ns, ] <- MMSE@Catch[, p_NOS_return, 1, 1, ]
  U_NOS[, ns, ] <- 1 - exp(-MMSE@FM[, p_NOS_return, 1, 1, ])

  do_hatchery <- SOM@n_subyearling > 0 || SOM@n_yearling > 0
  if (do_hatchery) {
    # HOS state variables from MMSE object
    p_HOS_return <- 5 # Population index for returning HOS
    p_HOS_escapement <- 6 # NOS escapement

    Return_HOS[, ns, nage, ] <- apply(MMSE@N[, p_HOS_return, nage, 1, , ], 1:2, sum)
    Escapement_HOS[, ns, nage, -SOM@proyears] <- apply(MMSE@N[, p_HOS_escapement, age_escapement, 1, -1, ], 1:2, sum)
    Catch_HOS[, ns, ] <- MMSE@Catch[, p_HOS_return, 1, 1, ]
    U_HOS[, ns, ] <- 1 - exp(-MMSE@FM[, p_HOS_return, 1, 1, ])

    # NOS + HOS state variables from salmonMSE
    ngen <- length(unique(salmonMSE_env$N$t))
    if (sum(y_spawn) != ngen) warning("Number of generations in salmonMSE state variables does not match generations in openMSE")

    NOS[, ns, y_spawn] <- get_salmonMSE_var(N, var = "NOS")
    Fry_NOS[, ns, y_spawn] <-  get_salmonMSE_var(N, var = "fry_NOS")
    Smolt_NOS[, ns, which(y_spawn) + 1] <- get_salmonMSE_var(N, var = "smolt_NOS")

    HOS[, ns, y_spawn] <- get_salmonMSE_var(N, var = "HOS")
    HOS_effective[, ns, y_spawn] <- get_salmonMSE_var(N, var = "HOS_effective")

    Fry_HOS[, ns, y_spawn] <- get_salmonMSE_var(N, var = "fry_HOS")
    Smolt_HOS[, ns, which(y_spawn) + 1] <- get_salmonMSE_var(N, var = "smolt_HOS")

    # Broodtake & fitness
    NOB[, ns, y_spawn] <- get_salmonMSE_var(N, var = "NOB")
    HOB[, ns, y_spawn] <- get_salmonMSE_var(N, var = "HOB")

    fitness[, ns, y_spawn] <- get_salmonMSE_var(N, var = "fitness") # Lag by generation?

    # Smolt releases and SAR loss from openMSE
    p_smolt_rel <- 4
    Smolt_Rel[, ns, ]  <- apply(MMSE@N[, p_smolt_rel, 1, 1, , ], 1:2, sum)

    p_smolt <- 1
    SAR_loss[, ns, ] <- exp(-MMSE@Misc$MICE$M_ageArray[, p_smolt, age_escapement - 2, 1, ])

    pNOB <- get_salmonMSE_var(N, var = "pNOB")
    pHOS <- get_salmonMSE_var(N, var = "pHOS")

    PNI[, ns, y_spawn] <- pNOB/(pNOB + pHOS) # Withler et al. 2018, page 17
    p_wild[, ns, y_spawn] <- sapply(1:ncol(pHOS), function(t) {
      if (t == 1) {
        rep(NA_real_, nrow(pHOS))
      } else {
        calc_pwild(pHOS[, t], pHOS[, t-1], SOM@gamma)
      }
    })

  } else {
    # If no hatchery, the NOS escapement is the NOS, fry_NOS is the spawning output
    NOS[] <- apply(Escapement_NOS, c(1, 2, 4), sum)

    p_smolt <- 1 # MSEtool population index for immature NOS
    p_spawn <- 3 # NOS escapement

    Fry_NOS[, ns, -SOM@proyears] <- MMSE@SSB[, p_spawn, 1, -1]
    Smolt_NOS[, ns, ] <- apply(MMSE@N[, p_smolt, 1, 1, , ], 1:2, sum)

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
    Catch_NOS = Catch_NOS,
    Catch_HOS = Catch_HOS,
    U_NOS = U_NOS,
    U_HOS = U_HOS,
    fitness = fitness,
    PNI = PNI,
    p_wild = p_wild,
    SAR_loss = SAR_loss
  )
  if (!missing(Harvest_MMP)) SMSE@Misc$Harvest_MMP <- Harvest_MMP
  if (!missing(Ford) && nrow(Ford)) SMSE@Misc$Ford <- salmonMSE_env$Ford # t is generation

  return(SMSE)
}


#' @importFrom dplyr filter
#' @importFrom reshape2 acast
get_salmonMSE_var <- function(N, var = "fry_NOS", p_smolt = 1) {
  filter(N, .data$p_smolt == .env$p_smolt) %>% reshape2::acast(list("x", "t"), value.var = var)
}


# Withler et al. 2018, page 27
calc_pwild <- function(pHOS_cur, pHOS_prev, gamma) {
  num <- (1 - pHOS_prev)^2
  denom <- num + 2 * gamma * pHOS_prev * (1 - pHOS_prev) + gamma * gamma * pHOS_prev^2
  (1 - pHOS_cur) * num/denom
}
