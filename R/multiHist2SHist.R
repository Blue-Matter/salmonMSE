
#' @rdname salmonMSE-int
#' @param multiHist Class multiHist object returned from MSEtool
#' @return
#' `multiHist2SHist`: \linkS4class{SHist} object
#' @export
multiHist2SHist <- function(multiHist, SOM, check = TRUE) {

  if (check) SOM <- check_SOM(SOM)

  ns <- length(SOM@Bio) # Number of stocks
  nage <- SOM@Bio[[1]]@maxage

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

  # openMSE year and age indices
  t1 <- seq(1, 2 * SOM@nyears, 2)
  t2 <- seq(2, 2 * SOM@nyears, 2)

  a1 <- seq(1, 2 * nage + 1, 2)
  a2 <- seq(2, 2 * nage + 1, 2)

  a_imm <- a1[-length(a1)]
  a_return <- a2
  a_esc <- a1[-1]
  #age_escapement <- age_bio + 1 # MSEtool age class

  f <- 1 # Fleet

  pindex <- make_stock_index(SOM)

  for (s in 1:ns) {

    # NOS state variables from MMSE object
    p_NOS_imm <- pindex$p[pindex$s == s & pindex$origin == "natural" & pindex$stage == "juvenile"]
    p_NOS_return <- pindex$p[pindex$s == s & pindex$origin == "natural" & pindex$stage == "recruitment"] # MSEtool population index for returning NOS
    p_NOS_escapement <- pindex$p[pindex$s == s & pindex$origin == "natural" & pindex$stage == "escapement"] # NOS escapement

    y_spawnOM <- which(rowSums(multiHist[[p_NOS_escapement[1]]][[f]]@TSdata$SBiomass[1, , ]) > 0)

    # Reported spawning is delayed by one year
    t1_sp <- t1[-1]
    t2_sp <- t2[-1]
    y_spawn <- 0.5 * (y_spawnOM - 1)

    Njuv_NOS[, s, , ] <- sapply(p_NOS_imm, function(p) {
      apply(multiHist[[p]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
    }, simplify = "array") %>%
      apply(1:3, sum)

    Return_NOS[, s, , ] <- sapply(p_NOS_return, function(p) {
      apply(multiHist[[p]][[f]]@AtAge$Number[, a_return, t2, ], 1:3, sum)
    }, simplify = "array") %>%
      apply(1:3, sum)

    Escapement_NOS[, s, , 1:length(t1_sp)] <- sapply(p_NOS_escapement, function(p) {
      apply(multiHist[[p]][[f]]@AtAge$Number[, a_esc, t1_sp, , drop = FALSE], 1:3, sum)
    }, simplify = "array") %>%
      apply(1:3, sum)

    # Kept catch
    KPT_NOS[, s, ] <- sapply(p_NOS_imm, function(p) {
      apply(multiHist[[p]][[f]]@TSdata$Landings[, t1, ], 1:2, sum)
    }, simplify = "array") %>%
      apply(1:2, sum)

    KT_NOS[, s, ] <- sapply(p_NOS_return, function(p) {
      apply(multiHist[[p]][[f]]@TSdata$Landings[, t2, ], 1:2, sum)
    }, simplify = "array") %>%
      apply(1:2, sum)

    # Total discards (live + dead)
    DPT_NOS[, s, ] <- sapply(p_NOS_imm, function(p) {
      apply(multiHist[[p]][[f]]@TSdata$Discards[, t1, ], 1:2, sum)
    }, simplify = "array") %>%
      apply(1:2, sum)

    DT_NOS[, s, ] <- sapply(p_NOS_return, function(p) {
      apply(multiHist[[p]][[f]]@TSdata$Discards[, t2, ], 1:2, sum)
    }, simplify = "array") %>%
      apply(1:2, sum)

    # Harvest rate from kept catch
    vulPT <- array(SOM@Harvest[[s]]@vulPT, c(nage, SOM@nsim, length(t1))) %>% aperm(c(2, 1, 3))
    vulNOS_imm <- apply(vulPT * Njuv_NOS[, s, , ], c(1, 3), sum)
    UPT_NOS[, s, ] <- KPT_NOS[, s, ]/vulNOS_imm
    UPT_NOS[is.na(UPT_NOS)] <- 0

    vulT <- array(SOM@Harvest[[s]]@vulT, c(nage, SOM@nsim, length(t2))) %>% aperm(c(2, 1, 3))
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

    Egg_NOS[, s, 1:length(t1_sp)] <- sapply(p_NOS_escapement, function(p) {
      apply(multiHist[[p]][[f]]@TSdata$SBiomass[, t1_sp, , drop = FALSE], 1:2, sum)
    }, simplify = "array") %>%
      apply(1:2, sum)

    Smolt[, s, ] <- sapply(p_NOS_imm, function(p) {
      apply(multiHist[[p]][[f]]@AtAge$Number[, 1, t1, ], 1:2, sum)
    }, simplify = "array") %>%
      apply(1:2, sum)

    do_hatchery <- SOM@Hatchery[[s]]@n_subyearling > 0 || SOM@Hatchery[[s]]@n_yearling > 0
    if (do_hatchery) {

      # HOS state variables from MMSE object
      p_HOS_imm <- pindex$p[pindex$s == s & pindex$origin == "hatchery" & pindex$stage == "juvenile"]
      p_HOS_return <- pindex$p[pindex$s == s & pindex$origin == "hatchery" & pindex$stage == "recruitment"] # MSEtool population index for returning HOS
      p_HOS_escapement <- pindex$p[pindex$s == s & pindex$origin == "hatchery" & pindex$stage == "escapement"] # HOS escapement

      Njuv_HOS[, s, , ] <- apply(multiHist[[p_HOS_imm]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
      Return_HOS[, s, , ] <- apply(multiHist[[p_HOS_return]][[f]]@AtAge$Number[, a_return, t2, ], 1:3, sum)
      Escapement_HOS[, s, , 1:length(t1_sp)] <- apply(multiHist[[p_HOS_escapement]][[f]]@AtAge$Number[, a_esc, t1_sp, , drop = FALSE], 1:3, sum)

      # Kept catch
      KPT_HOS[, s, ] <- apply(multiHist[[p_HOS_imm]][[f]]@TSdata$Landings[, t1, ], 1:2, sum)
      KT_HOS[, s, ] <- apply(multiHist[[p_HOS_return]][[f]]@TSdata$Landings[, t2, ], 1:2, sum)

      # Total discards (live + dead)
      DPT_HOS[, s, ] <- apply(multiHist[[p_HOS_imm]][[f]]@TSdata$Discards[, t1, ], 1:2, sum)
      DT_HOS[, s, ] <- apply(multiHist[[p_HOS_return]][[f]]@TSdata$Discards[, t2, ], 1:2, sum)

      # Harvest rate from kept catch
      HOS_imm_a <- apply(multiHist[[p_HOS_imm]][[f]]@AtAge$Number[, a_imm, t1, ], 1:3, sum)
      vulHOS_imm <- apply(vulPT * HOS_imm_a, c(1, 3), sum)
      UPT_HOS[, s, ] <- KPT_HOS[, s, ]/vulHOS_imm
      UPT_HOS[is.na(UPT_HOS)] <- 0

      HOS_ret_a <- apply(multiHist[[p_HOS_return]][[f]]@AtAge$Number[, a_imm, t2, ], 1:3, sum)
      vulHOS_ret <- apply(vulT * HOS_ret_a, c(1, 3), sum)
      UT_HOS[, s, ] <- KT_HOS[, s, ]/vulHOS_ret
      UT_HOS[is.na(UT_HOS)] <- 0

      # Exploitation rate from kept + dead discards (DD)
      DDPT_HOS <- SOM@Harvest[[s]]@release_mort[1] * DPT_HOS[, s, ]
      ExPT_HOS[, s, ] <- (KPT_HOS[, s, ] + DDPT_HOS)/vulHOS_imm
      ExPT_HOS[is.na(ExPT_HOS)] <- 0

      DDT_HOS <- SOM@Harvest[[s]]@release_mort[2] * DT_HOS[, s, ]
      ExT_HOS[, s, ] <- (KT_HOS[, s, ] + DDT_HOS)/vulHOS_ret
      ExT_HOS[is.na(ExT_HOS)] <- 0

      Egg_HOS[, s, 1:length(t1_sp)] <- apply(multiHist[[p_HOS_escapement]][[f]]@TSdata$SBiomass[, t1_sp, , drop = FALSE], 1:2, sum)
      Smolt_Rel[, s, ] <- apply(multiHist[[p_HOS_imm]][[f]]@AtAge$Number[, 1, t1, ], c(1, 2), sum)
    }

  }

  SHist <- new(
    "SHist",
    Name = SOM@Name,
    nyears = SOM@nyears,
    nsim = SOM@nsim,
    nstocks = ns,
    Snames = sapply(1:ns, function(s) if (length(SOM@Bio[[s]]@Name)) SOM@Bio[[s]]@Name else paste("Population", s)),
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
