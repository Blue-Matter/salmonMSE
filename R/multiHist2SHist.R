
#' @rdname salmonMSE-int
#' @param multiHist Class multiHist object returned from MSEtool
#' @return
#' `multiHist2SHist`: \linkS4class{SHist} object
#' @export
multiHist2SHist <- function(multiHist, SOM, check = TRUE) {

  if (check) SOM <- check_SOM(SOM)

  ns <- 1 # Number of stocks
  nage <- SOM@maxage

  # Declare arrays
  Njuv_NOS <- Njuv_HOS <- Escapement_NOS <- Escapement_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@nyears))
  Return_NOS <- Return_HOS <- array(0, c(SOM@nsim, ns, nage, SOM@nyears))
  KPT_NOS <- KT_NOS <- KPT_HOS <- KT_HOS <- array(0, c(SOM@nsim, ns, SOM@nyears))
  DPT_NOS <- DT_NOS <- DPT_HOS <- DT_HOS <- array(0, c(SOM@nsim, ns, SOM@nyears))
  UPT_NOS <- UT_NOS <- UPT_HOS <- UT_HOS <- array(0, c(SOM@nsim, ns, SOM@nyears))
  ExPT_NOS <- ExT_NOS <- ExPT_HOS <- ExT_HOS <- array(0, c(SOM@nsim, ns, SOM@nyears))

  Egg_NOS <- Egg_HOS <- array(0, c(SOM@nsim, ns, SOM@nyears))
  Smolt <- Smolt_Rel <- array(0, c(SOM@nsim, ns, SOM@nyears))

  NOS <- HOS <- HOS_effective <- array(0, c(SOM@nsim, ns, SOM@nyears))

  # NOS state variables from MMSE object
  p_NOS_imm <- 1
  p_NOS_return <- 2 # MSEtool population index for returning NOS
  p_NOS_escapement <- 3 # NOS escapement

  # openMSE year and age indices
  t1 <- seq(1, 2 * SOM@nyears, 2)
  t2 <- seq(2, 2 * SOM@nyears, 2)

  a1 <- seq(1, 2 * SOM@maxage + 1, 2)
  a2 <- seq(2, 2 * SOM@maxage + 1, 2)

  a_imm <- a1[-length(a1)]
  a_return <- a2
  a_esc <- a1[-1]
  #age_escapement <- age_bio + 1 # MSEtool age class

  f <- 1 # Fleet

  y_spawnOM <- which(rowSums(multiHist[[p_NOS_escapement]][[f]]@TSdata$SBiomass[1, , ]) > 0)

  # Reported spawning is delayed by one year
  t1_sp <- t1[-1]
  t2_sp <- t2[-1]
  y_spawn <- 0.5 * (y_spawnOM - 1)

  Njuv_NOS[, ns, , ] <- apply(multiHist[[p_NOS_imm]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
  Return_NOS[, ns, , ] <- apply(multiHist[[p_NOS_return]][[f]]@AtAge$Number[, a_return, t2, ], 1:3, sum)
  Escapement_NOS[, ns, , 1:length(t1_sp)] <- apply(multiHist[[p_NOS_escapement]][[f]]@AtAge$Number[, a_esc, t1_sp, ], 1:3, sum)

  # Kept catch
  KPT_NOS[, ns, ] <- apply(multiHist[[p_NOS_imm]][[f]]@TSdata$Landings[, t1, ], 1:2, sum)
  KT_NOS[, ns, ] <- apply(multiHist[[p_NOS_return]][[f]]@TSdata$Landings[, t2, ], 1:2, sum)

  # Total discards (live + dead)
  DPT_NOS[, ns, ] <- apply(multiHist[[p_NOS_imm]][[f]]@TSdata$Discards[, t1, ], 1:2, sum)
  DT_NOS[, ns, ] <- apply(multiHist[[p_NOS_return]][[f]]@TSdata$Discards[, t2, ], 1:2, sum)

  # Harvest rate from kept catch
  vulPT <- array(SOM@vulPT, c(SOM@maxage, SOM@nsim, length(t1))) %>% aperm(c(2, 1, 3))
  NOS_imm_a <- apply(multiHist[[p_NOS_imm]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
  vulNOS_imm <- apply(vulPT * NOS_imm_a, c(1, 3), sum)
  UPT_NOS[, ns, ] <- KPT_NOS[, ns, ]/vulNOS_imm
  UPT_NOS[is.na(UPT_NOS)] <- 0

  vulT <- array(SOM@vulT, c(SOM@maxage, SOM@nsim, length(t2))) %>% aperm(c(2, 1, 3))
  NOS_ret_a <- apply(multiHist[[p_NOS_return]][[f]]@AtAge$Number[, a_return, t2, ], 1:3, sum)
  vulNOS_ret <- apply(vulT * NOS_ret_a, c(1, 3), sum)
  UT_NOS[, ns, ] <- KT_NOS[, ns, ]/vulNOS_ret
  UT_NOS[is.na(UT_NOS)] <- 0

  # Exploitation rate from kept + dead discards (DD)
  DDPT_NOS <- SOM@release_mort[1] * DPT_NOS[, ns, ]
  ExPT_NOS[, ns, ] <- (KPT_NOS[, ns, ] + DDPT_NOS)/vulNOS_imm
  ExPT_NOS[is.na(ExPT_NOS)] <- 0

  DDT_NOS <- SOM@release_mort[2] * DT_NOS[, ns, ]
  ExT_NOS[, ns, ] <- (KT_NOS[, ns, ] + DDT_NOS)/vulNOS_ret
  ExT_NOS[is.na(ExT_NOS)] <- 0

  Egg_NOS[, ns, 1:length(t1_sp)] <- apply(multiHist[[p_NOS_escapement]][[f]]@TSdata$SBiomass[, t1_sp, ], 1:2, sum)
  Smolt[, ns, ] <- apply(multiHist[[p_NOS_imm]][[f]]@AtAge$Number[, 1, t1, ], c(1, 2), sum)

  do_hatchery <- SOM@n_subyearling > 0 || SOM@n_yearling > 0
  if (do_hatchery) {
    # HOS state variables from MMSE object
    p_HOS_imm <- 4
    p_HOS_return <- 5 # Population index for returning HOS
    p_HOS_escapement <- 6 # NOS escapement

    Njuv_HOS[, ns, , ] <- apply(multiHist[[p_HOS_imm]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
    Return_HOS[, ns, , ] <- apply(multiHist[[p_HOS_return]][[f]]@AtAge$Number[, a_return, t2, ], 1:3, sum)
    Escapement_HOS[, ns, , 1:length(t1_sp)] <- apply(multiHist[[p_HOS_escapement]][[f]]@AtAge$Number[, a_esc, t1_sp, ], 1:3, sum)

    # Kept catch
    KPT_HOS[, ns, ] <- apply(multiHist[[p_HOS_imm]][[f]]@TSdata$Landings[, t1, ], 1:2, sum)
    KT_HOS[, ns, ] <- apply(multiHist[[p_HOS_return]][[f]]@TSdata$Landings[, t2, ], 1:2, sum)

    # Total discards (live + dead)
    DPT_HOS[, ns, ] <- apply(multiHist[[p_HOS_imm]][[f]]@TSdata$Discards[, t1, ], 1:2, sum)
    DT_HOS[, ns, ] <- apply(multiHist[[p_HOS_return]][[f]]@TSdata$Discards[, t2, ], 1:2, sum)

    # Harvest rate from kept catch
    HOS_imm_a <- apply(multiHist[[p_HOS_imm]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
    vulHOS_imm <- apply(vulPT * HOS_imm_a, c(1, 3), sum)
    UPT_HOS[, ns, ] <- KPT_HOS[, ns, ]/vulHOS_imm
    UPT_HOS[is.na(UPT_HOS)] <- 0

    HOS_ret_a <- apply(multiHist[[p_HOS_return]][[f]]@AtAge$Number[, a_imm, t2, ], 1:3, sum)
    vulHOS_ret <- apply(vulT * HOS_ret_a, c(1, 3), sum)
    UT_HOS[, ns, ] <- KT_HOS[, ns, ]/vulHOS_ret
    UT_HOS[is.na(UT_HOS)] <- 0

    # Exploitation rate from kept + dead discards (DD)
    DDPT_HOS <- SOM@release_mort[1] * DPT_HOS[, ns, ]
    ExPT_HOS[, ns, ] <- (KPT_HOS[, ns, ] + DDPT_HOS)/vulHOS_imm
    ExPT_HOS[is.na(ExPT_HOS)] <- 0

    DDT_HOS <- SOM@release_mort[2] * DT_HOS[, ns, ]
    ExT_HOS[, ns, ] <- (KT_HOS[, ns, ] + DDT_HOS)/vulHOS_ret
    ExT_HOS[is.na(ExT_HOS)] <- 0

    Egg_HOS[, ns, 1:length(t1_sp)] <- apply(multiHist[[p_HOS_escapement]][[f]]@TSdata$SBiomass[, t1_sp, ], 1:2, sum)
    Smolt_Rel[, ns, ] <- apply(multiHist[[p_HOS_imm]][[f]]@AtAge$Number[, 1, t1, ], c(1, 2), sum)
  }

  SHist <- new(
    "SHist",
    Name = SOM@Name,
    nyears = SOM@nyears,
    nsim = SOM@nsim,
    nstocks = ns,
    Snames = "character",
    Egg_NOS = Egg_NOS,
    Egg_HOS = Egg_HOS,
    Smolt = Smolt,
    Smolt_Rel = Smolt_Rel,
    Njuv_NOS = Njuv_NOS,
    Njuv_HOS = Njuv_HOS,
    Return_NOS = Return_NOS,
    Return_HOS = Return_HOS,
    Escapement_NOS = Escapement_NOS,
    Escapement_HOS = Escapement_HOS,
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
    ExT_HOS = ExT_HOS
  )

  return(SHist)
}
