
#' Reference points
#'
#' Calculate MSY and Sgen reference points for the operating model. Uses the biological parameters (maturity, natural mortality)
#' in the last year of the projection.
#'
#' @param SOM An object of class [salmonMSE::SOM-class]
#' @param rel_F Numeric length 2, indicates the relative effort in the preterminal and terminal fisheries, with a maximum value of 1.
#' For example, `c(1, 0)` indicates a yield calculation for only the preterminal fishery. By default, uses the ratio of the harvest rates in the operating model.
#' @param check Logical, whether to check the SOM object using [check_SOM()]
#' @returns Matrix of various state variables (catch, exploitation rate, egg production, spawners) at MSY and Sgen by simulation
#' @importFrom stats optimize
#' @export
calc_ref <- function(SOM, rel_F, check = TRUE) {
  if (check) SOM <- check_SOM(SOM)

  ns <- length(SOM@Bio)

  y <- SOM@nyears + SOM@proyears

  if (missing(rel_F)) {
    rel_F_s <- lapply(1:ns, function(s) {
      x <- c(SOM@Harvest[[s]]@u_preterminal, SOM@Harvest[[s]]@u_terminal)
      x/max(x)
    })
  } else {
    rel_F_s <- lapply(1:ns, function(...) rel_F)
  }

  ref_s <- lapply(1:ns, function(s) {
    SRR <- make_SRR(SOM@Bio[[s]], SOM@Habitat[[s]])

    Mjuv_NOS <- SOM@Bio[[s]]@Mjuv_NOS
    fec <- SOM@Bio[[s]]@fec
    p_female <- SOM@Bio[[s]]@p_female
    vulPT <- SOM@Harvest[[s]]@vulPT
    vulT <- SOM@Harvest[[s]]@vulT
    p_mature <- SOM@Bio[[s]]@p_mature

    val <- sapply(1:SOM@nsim, function(x) {
      opt <- optimize(
        .calc_eq, interval = c(1e-8, 3),
        Mjuv = Mjuv_NOS[x, , y], fec = fec, p_female = p_female, rel_F = rel_F_s[[s]],
        vulPT = vulPT, vulT = vulT, p_mature = p_mature[x, , y], s_enroute = 1,
        SRRpars = SRR$SRRpars[x, ],
        maximum = TRUE
      )
      ref <- .calc_eq(
        .F = opt$maximum,
        Mjuv = Mjuv_NOS[x, , y], fec = fec, p_female = p_female, rel_F = rel_F_s[[s]],
        vulPT = vulPT, vulT = vulT, p_mature = p_mature[x, , y], s_enroute = 1,
        SRRpars = SRR$SRRpars[x, ], opt = FALSE
      ) %>% unlist()

      names(ref) <- paste0(names(ref), "_MSY")
      if (ref["KPT_MSY"] == 0) ref["UPT_MSY"] <- 0
      if (ref["KT_MSY"] == 0) ref["UT_MSY"] <- 0

      if (any(ref < 0, na.rm = TRUE)) {

        ref[!is.na(ref)] <- 0
        ref["Sgen"] <- 0

      } else {

        opt_Sgen <- try(
          uniroot(
            .calc_Sgen, interval = c(1, 100) * opt$maximum,
            Mjuv = Mjuv_NOS[x, , y], fec = fec, p_female = p_female, rel_F = rel_F_s[[s]],
            vulPT = vulPT, vulT = vulT, p_mature = p_mature[x, , y], s_enroute = 1,
            SRRpars = SRR$SRRpars[x, ], SMSY = ref["Spawners_MSY"]
          ),
          silent = TRUE
        )

        if (is.character(opt_Sgen)) {
          ref["Sgen"] <- NA
        } else {
          Sgen_out <- .calc_Sgen(
            .F = opt_Sgen$root,
            Mjuv = Mjuv_NOS[x, , y], fec = fec, p_female = p_female, rel_F = rel_F_s[[s]],
            vulPT = vulPT, vulT = vulT, p_mature = p_mature[x, , y], s_enroute = 1,
            SRRpars = SRR$SRRpars[x, ], SMSY = ref["Spawners"], opt = FALSE
          )

          ref["Sgen"] <- sum(Sgen_out$Spawners)
        }

      }
      return(ref)
    })

    return(val)

  })

  return(ref_s)
}


.calc_Sgen <- function(.F, Mjuv, fec, p_female, gamma = 1, rel_F, vulPT, vulT, p_mature, s_enroute = 1,
                       SRRpars, popt = 1, SMSY, opt = TRUE) {

  Sgen_ref <- .calc_eq(
    .F, Mjuv, fec, p_female, gamma, rel_F, vulPT, vulT, p_mature, s_enroute,
    SRRpars, opt = FALSE, aggregate_age = FALSE
  )

  if (Sgen_ref$Egg <= 0) {
    Sgen_plusone <- 0
  } else {
    nyears <- maxage <- length(Mjuv)

    Njuv <- Spawner <- matrix(0, maxage, nyears)
    Egg <- rep(0, nyears)
    Njuv[, 1] <- Sgen_ref$Njuv

    for (y in 1:nyears) {
      Spawner[, y] <- Njuv[, y] * p_mature * s_enroute
      Egg[y] <- sum(Spawner[, y] * gamma * p_female * fec)

      if (y < nyears) {
        Njuv[1, y+1] <- calc_smolt(
          Egg[y],
          kappa = SRRpars["kappa"], capacity = SRRpars["capacity_smolt"], Smax = SRRpars["Smax"], phi = SRRpars["phi"],
          kappa_improve = SRRpars["kappa_improve"], capacity_improve = SRRpars["capacity_smolt_improve"], fitness_loss = 1,
          SRrel = as.character(SRRpars["SRrel"])
        ) %>%
          as.numeric()
        Njuv[2:maxage, y+1] <- Njuv[2:maxage - 1, y] * exp(-Mjuv[2:maxage - 1]) * (1 - p_mature[2:maxage - 1])
      }
    }

    Sgen_plusone <- sum(Spawner[, nyears])
  }

  if (opt) {
    return(Sgen_plusone - SMSY)
  } else {
    return(Sgen_ref)
  }
}


.calc_eq <- function(.F, Mjuv, fec, p_female, gamma = 1, rel_F, vulPT, vulT, p_mature, s_enroute = 1,
                     SRRpars, opt = TRUE, aggregate_age = TRUE) {

  FPT <- vulPT * .F * rel_F[1]
  FT <- vulT * .F * rel_F[2]

  surv_juv <- calc_survival(Mjuv + FPT, p_mature)
  surv_return <- surv_juv * p_mature
  surv_esc <- surv_return * exp(-FT)
  surv_spawn <- surv_esc * s_enroute

  EPR <- sum(surv_spawn * gamma * p_female * fec)

  Smolt <- calc_smolt(
    EPR,
    kappa = SRRpars["kappa"], capacity = SRRpars["capacity_smolt"], Smax = SRRpars["Smax"], phi = SRRpars["phi"],
    kappa_improve = SRRpars["kappa_improve"], capacity_improve = SRRpars["capacity_smolt_improve"], fitness_loss = 1,
    SRrel = as.character(SRRpars["SRrel"]), per_recruit = TRUE
  ) %>%
    as.numeric()

  Njuv <- Smolt * surv_juv
  Return <- Smolt * surv_return

  KPT <- Njuv * (1 - exp(-FPT))
  KT <- Return * (1 - exp(-FT))

  if (opt) {
    return(sum(KPT, KT))
  } else {

    output <- list(
      KPT = KPT,
      KT = KT,
      FPT = max(FPT),
      FT = max(FT),
      UPT = sum(KPT)/sum(vulPT * Njuv),
      UT = sum(KT)/sum(vulT * Return),
      Smolt = Smolt,
      Njuv = Njuv,
      Return = Return,
      Escapement = Smolt * surv_esc,
      Spawners = Smolt * surv_spawn,
      Egg = Smolt * EPR
    )

    if (aggregate_age) output <- lapply(output, sum)

    return(output)
  }
}
