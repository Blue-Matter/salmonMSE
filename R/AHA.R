

#' All-H Analyzer
#'
#' Wrapper function for an implementation of All-H Analyzer ([AHA](https://www.streamnet.org/home/data-maps/hatchery-reform/hsrg-tools/)) in R. Can be used to compare outputs between AHA
#' and salmonMSE.
#'
#' @param SOM An object of class \linkS4class{SOM}
#' @param ngen Integer, the number of generations for which to run the simulation
#' @param silent Logical, indicates whether to silence messages to the R console
#' @return A named list containing vectors of state variables (by simulation, population, and generation). See \linkS4class{SMSE} object description.
#'
#' @references Hatchery Scientific Review Group. 2020. All-H Analyzer Tool Guide and Documentation. May 2020.
#' @export
AHA <- function(SOM, ngen = 100, silent = FALSE) {
  SOM <- check_SOM(SOM)
  ns <- length(SOM@Bio)

  output_s <- list()

  for (s in 1:ns) {

    if (!silent) message("Running AHA for population ", s)

    if (sum(SOM@Harvest[[s]]@u_preterminal)) {
      warning("Pre-terminal fishing is not modeled in AHA. Setting SOM@Harvest[[", s, "]]@u_preterminal = 0")
      SOM@Harvest[[s]]@u_preterminal <- 0
    }

    if (SOM@Bio[[s]]@SRrel != "BH") {
      warning("Only Beverton-Holt smolt production is used in AHA. Setting SOM@Bio[[", s, "]]@SRrel = \"BH\"")
      SOM@Bio[[s]]@SRrel <- "BH"
    }

    if (!length(SOM@Bio[[s]]@capacity)) {
      stop("Need to specify SOM@Bio[[", s, "]]@capacity")
    }

    if (SOM@Hatchery[[s]]@fitness_type[2] != "none") {
      warning("Only fitness dynamics in the natural environment will be modeled.")
    }

    if (SOM@Hatchery[[s]]@pmax_esc < 1) {
      warning("SOM@Hatchery[[", s, "]]@pmax_esc = ", SOM@Hatchery@pmax_esc,
              ". AHA compatibility is maintained only when SOM@Hatchery@pmax_esc = 1")
    }

    if (SOM@Hatchery[[s]]@m < 1) {
      warning("SOM@Hatchery[[", s, "]]@m = ", SOM@Hatchery@m, ". AHA compatibility is maintained only when mark rate m = 1 for non-selective broodstock.")
    }

    age_mature <- SOM@Bio[[s]]@p_mature[1, , 1] > 0
    age_mature <- which(age_mature)[1]
    message("Age of maturity assumed to be: ", age_mature)

    fec <- SOM@Bio[[s]]@fec[age_mature]
    SOM@Bio[[s]]@fec <- fec
    message("Fecundity of spawners assumed to be: ", fec)

    do_hatchery <- sum(SOM@Hatchery[[s]]@n_subyearling, SOM@Hatchery[[s]]@n_yearling) > 0
    if (do_hatchery) {
      fec_brood <- SOM@Hatchery[[s]]@fec_brood[age_mature]
      message("Fecundity of broodtake assumed to be: ", fec_brood)

      if (is.na(SOM@Hatchery[[s]]@phatchery)) stop("Need value for SOM@Hatchery[[", s, "]]@phatchery")
    } else {
      fec_brood <- 0
    }
    SOM@Hatchery[[s]]@fec_brood <- fec_brood

    message("SAR calculated from survival from Mjuv to age ", age_mature)

    if (SOM@Bio[[s]]@n_g > 1) {
      warning("Multiple life history groups detected: SOM@Bio[[", s, "]]@n_g > 1. Using parameters for group 1")
    }
    if (SOM@Hatchery[[s]]@n_r > 1) {
      warning("Multiple hatchery release groups detected: SOM@Hatchery[[", s, "]]@n_r > 1. Using parameters for group 1")
    }

    surv_NOS <- surv_HOS <- matrix(1, SOM@nsim, age_mature)
    for (a in 2:age_mature) {
      surv_NOS[, a] <- surv_NOS[, a-1] * exp(-SOM@Bio[[s]]@Mjuv_NOS[, a-1, 1, 1])
      if (do_hatchery) {
        surv_HOS[, a] <- surv_HOS[, a-1] * exp(-SOM@Hatchery[[s]]@Mjuv_HOS[, a-1, 1, 1])
      }
    }

    if (SOM@Harvest[[s]]@MSF_T) {
      .F <- get_F(u = SOM@Harvest[[s]]@u_terminal, M = 1e-8, ret = SOM@Hatchery[[s]]@m, release_mort = SOM@Harvest[[s]]@release_mort[2])
      Frel <- (1 - SOM@Hatchery[[s]]@m) * SOM@Harvest[[s]]@release_mort[2] * .F
      Fret <- SOM@Hatchery[[s]]@m * .F
      u_NOR <- 1 - exp(-Frel)
      u_HOR <- 1 - exp(-Frel - Fret)
      message(
        "Mark-selective fishing for terminal fishery detected. For AHA, setting harvest rate of natural and hatchery fish = ",
        round(u_NOR, 3), " and ", round(u_HOR, 3), ", respectively. Catch represents kept catch + dead releases."
      )
    }

    output_sim <- lapply(1:SOM@nsim, .AHA_wrapper, s = s, SOM = SOM, ngen = ngen, SAR_NOS = surv_NOS[, age_mature], SAR_HOS = surv_HOS[, age_mature])
    var_report <- names(output_sim[[1]])

    output_s[[s]] <- lapply(var_report, function(i) {
      sapply(output_sim, getElement, i) %>% t()
    }) %>%
      structure(names = var_report)

  }

  output <- lapply(var_report, function(i) {
    mat <- sapply(output_s, getElement, i, simplify = "array") %>%
      aperm(c(1, 3, 2))
  }) %>%
    structure(names = var_report)

  return(output)
}


.AHA_wrapper <- function(x, s, SOM, ngen, SAR_NOS, SAR_HOS, silent = FALSE) {

  Bio <- SOM@Bio[[s]]
  Habitat <- SOM@Habitat[[s]]
  Hatchery <- SOM@Hatchery[[s]]
  Harvest <- SOM@Harvest[[s]]

  if (SOM@Harvest[[s]]@MSF_T) {
    .F <- get_F(u = Harvest@u_terminal, M = 1e-8, ret = Hatchery@m, release_mort = Harvest@release_mort[2])
    Frel <- (1 - Hatchery@m) * Harvest@release_mort[2] * .F
    Fret <- Hatchery@m * .F
    u_NOR <- 1 - exp(-Frel)
    u_HOR <- 1 - exp(-Frel - Fret)
  } else {
    u_NOR <- u_HOR <- Harvest@u_terminal
  }

  output <- .AHA(
    prod_adult = Bio@kappa[x],
    capacity_adult = Bio@capacity[x] * SAR_NOS[x],
    fec_spawn = sum(Bio@fec),
    p_female = Bio@p_female,
    surv_ocean = SAR_NOS[x],
    surv_passage_juv = 1,
    surv_passage_adult = Bio@s_enroute,
    u_HOR = c(u_HOR, 0, 0, 0),
    u_NOR = c(u_NOR, 0, 0, 0),
    surv_pre_spawn = Hatchery@s_prespawn,
    fec_brood = sum(Hatchery@fec_brood),
    surv_egg_smolt = Hatchery@s_egg_smolt,
    surv_egg_subyearling = Hatchery@s_egg_subyearling,
    capacity_spawn_em = 1e12,
    capacity_smolt_adult = 1e12,
    surv_adult_return_of_yearling = SAR_HOS[x],
    surv_adult_return_of_subyearling = SAR_HOS[x],
    RRS_HOS = Hatchery@gamma,
    p_weir_efficiency = Hatchery@premove_HOS,
    p_return_hatchery = Hatchery@phatchery,
    n_release_yearling = sum(Hatchery@n_yearling),
    n_release_subyearling = sum(Hatchery@n_subyearling),
    brood_import = 0,
    brood_export = 0,
    p_NOR_brood_max = Hatchery@pmax_NOB,
    p_HOS_goal = 0, # Not used
    p_NOB_goal = Hatchery@ptarget_NOB,
    SAR_vary = "none",
    fitness_type = Hatchery@fitness_type[1], # Fitness
    theta = Hatchery@theta,
    rel_loss_egg = Hatchery@rel_loss[1],
    rel_loss_fry = Hatchery@rel_loss[2],
    rel_loss_smolt = Hatchery@rel_loss[3],
    Zpop_start = c(Hatchery@zbar_start[x, 1, ], 100),
    phenotype_variance = Hatchery@phenotype_variance,
    selection_strength = sqrt(Hatchery@fitness_variance/Hatchery@phenotype_variance),
    heritability = Hatchery@heritability[x],
    fitness_floor = Hatchery@fitness_floor,
    strays_total = 0,
    ngen = ngen
  )

  var_out <- c(
    "NOR" = "Return_NOS",
    "HOR" = "Return_HOS",
    "harvest_NOR" = "KT_NOS",
    "harvest_HOR" = "KT_HOS",
    "escapement_NOR" = "Escapement_NOS",
    "escapement_HOR" = "Escapement_HOS",
    "NOB" = "NOB",
    "HOB" = "HOB",
    "NOS" = "NOS",
    "HOS" = "HOS",
    "HOS_effective" = "HOS_effective",
    "fry_NOS" = "Egg_NOS",
    "fry_HOS" = "Egg_HOS",
    "fry_HOR" = "Egg_Rel",
    "smolt_NOS" = "Smolt_NOS",
    "smolt_HOS" = "Smolt_HOS",
    "smolt_HOR" = "Smolt_Rel",
    "fitness" = "fitness",
    "SAR_loss" = "SAR_loss",
    "pNOB" = "pNOB",
    "pHOSeff" = "pHOSeff",
    "alpha" = "alpha",
    "beta" = "beta"
  )

  out <- lapply(names(var_out), function(x) {
    sapply(output, getElement, x)
  }) %>%
    structure(names = var_out)
  out$PNI <- out$pNOB/(out$pNOB + out$pHOSeff)
  out$p_wild <- sapply(2:length(out$pHOSeff), function(g) {
    calc_pwild(out$pHOSeff[g], out$pHOSeff[g-1], Hatchery@gamma)
  })
  out$Ref <- output[[1]][["Ref"]]

  return(out)
}

.AHA <- function(prod_adult = 2.7,                     # AHA habitat fields
                 capacity_adult = 17250,
                 #capacity_smolt = 1.917e6,
                 fec_spawn = 5040,
                 p_female = 0.49,
                 surv_ocean = 0.01,                    # Ocean dynamics
                 surv_passage_juv = 1,
                 surv_passage_adult = 1,
                 u_HOR = 1e-2 * c(3.8, 2.5, 0, 14),   # Harvest dynamics
                 u_NOR = 1e-2 * c(3.8, 2.5, 0, 14),
                 surv_pre_spawn = 1,                   ## Survival prior to spawning
                 fec_brood = 5040,                       # In-hatchery
                 surv_egg_smolt = 1e-6,
                 surv_egg_subyearling = 0.92,
                 capacity_spawn_em = 1e12,                 ## Hard wired hatchery settings
                 capacity_smolt_adult = 1e12,
                 surv_adult_return_of_yearling = 1e-4,       # Post-hatchery
                 surv_adult_return_of_subyearling = 8.5e-3,
                 RRS_HOS = 0.8,
                 p_weir_efficiency = 0,
                 p_return_hatchery = 0.04,
                 n_release_yearling = 0,
                 n_release_subyearling = 2e6,
                 brood_import = 0,
                 brood_export = 0,
                 p_NOR_brood_max = 0.7,
                 p_HOS_goal = 0,
                 p_NOB_goal = 0.51,
                 SAR_vary = c("none", "PDO", "random"),
                 fitness_type = c("Ford", "Busack", "none"), # Fitness
                 theta = c(100, 80, 70),
                 rel_loss_egg = 0.5,
                 rel_loss_fry = 0.4,
                 rel_loss_smolt = 0.1,
                 Zpop_start = c(93.1, 92, 90),
                 phenotype_variance = 10,
                 selection_strength = 3,
                 heritability = 0.5,
                 fitness_floor = 0.5,
                 strays_total = 0,
                 ngen = 100) {

  SAR_vary <- match.arg(SAR_vary)
  stopifnot(SAR_vary == "none")

  fitness_type <- match.arg(fitness_type)
  if (fitness_type == "Busack") stop("Busack fitness calcs not completed yet")
  # Habitat section
  prod_spawn_em <- fec_spawn * p_female

  #habitat_h8 <- 0.009 # See Excel spreadsheet
  habitat_h8 <- surv_ocean
  surv_smolt_adult_obs <- surv_passage_juv * surv_ocean * surv_passage_adult # SAR, survival from release to adult return (22)
  surv_smolt_adult_base <- ifelse(habitat_h8 > 0, habitat_h8, surv_smolt_adult_obs) # (23)

  prod_smolt_adult <- surv_smolt_adult_base # (13)
  prod_em_smolt <- prod_adult/prod_spawn_em/prod_smolt_adult # (11)

  capacity_em_smolt <- local({ # (12)
    denom <- 1/capacity_adult-1/capacity_smolt_adult
    1/prod_smolt_adult/denom
  })

  # Reporting variables (5,6)
  adjust_prod_adult <- prod_adult * surv_smolt_adult_obs / surv_smolt_adult_base
  adjust_capacity_adult <- prod_em_smolt * surv_smolt_adult_obs/
    (1/capacity_spawn_em + prod_em_smolt/capacity_em_smolt + prod_em_smolt * surv_smolt_adult_obs/capacity_smolt_adult)

  # Ref pt (15-18)
  Neq <- adjust_capacity_adult * (1 - 1/adjust_prod_adult)
  MSY <- ((sqrt(adjust_prod_adult)-1)^2) * adjust_capacity_adult / adjust_prod_adult
  SMSY <- (sqrt(adjust_prod_adult) - 1) * adjust_capacity_adult /adjust_prod_adult
  umsy <- MSY/(MSY + SMSY)

  # (41-42)
  surv_adult_return_of_yearling_applied <- surv_adult_return_of_yearling *
    surv_smolt_adult_base/surv_smolt_adult_obs
  surv_adult_return_of_subyearling_applied <- surv_adult_return_of_subyearling *
    surv_smolt_adult_base/surv_smolt_adult_obs

  # (49-50)
  p_yearling <- n_release_yearling/(n_release_yearling + n_release_subyearling)
  p_subyearling <- 1 - p_yearling

  # (52)
  brood_local <- n_release_yearling/(surv_pre_spawn * p_female * fec_brood * surv_egg_smolt) +
    n_release_subyearling/(surv_pre_spawn * p_female * fec_brood * surv_egg_subyearling) - brood_import

  #SAR_rel (#56) Relative survival of hatchery origin fish after release from hatchery
  surv_adult_smolt_rel <- (surv_adult_return_of_yearling_applied * p_yearling +
                             surv_adult_return_of_subyearling_applied * p_subyearling)/surv_smolt_adult_base

  # Calc fitness variables (80,81)
  omega <- sqrt(phenotype_variance) * selection_strength
  omega2 <- omega * omega
  A <- 1 - heritability * phenotype_variance/(omega2 + phenotype_variance)

  egg_per_spawner <- fec_brood * surv_pre_spawn * p_female # E
  surv_egg_release <- surv_egg_subyearling * p_subyearling + surv_egg_smolt * p_yearling # F

  # BC
  surv_adult_return_apply <- numeric(ngen)
  if (SAR_vary == "none") {
    surv_adult_return_apply[] <- surv_ocean
  } else if (SAR_vary == "random") {
    #surv_random <- MIN(LOGINV(RAND(),LN(BB59),MAX(LN(2*$I$17+1)/2,0.0000000001)),0.99999)
    #surv_adult_return_apply[] <- surv_random
  } else { #PDO
    #surv_adult_return_apply[] <- surv_ocean * SAR_V8
  }

  # Geometric mean BC55
  geomean <- function(x) exp(mean(log(x)))
  surv_adult_return_geo <- geomean(surv_adult_return_apply[1:ngen])
  #surv_adult_return_geo <- geomean(surv_adult_return_apply[21:100])

  # G (hatchery survival from release to adulthood)
  surv_release_adult <- surv_adult_return_apply * surv_adult_smolt_rel * surv_passage_juv

  ###### AHA loop
  AHA_loop <- list()

  # Fitness calculation
  pbar <- matrix(NA, 2, ngen)

  for(g in 1:ngen) {

    if (g == 1) {
      HOR_total <- NOR_total <- 1e3
      #HOR_total <- 0
      #phi0 <- surv_smolt_adult_base * p_female * fec_spawn
      #SRRpars <- calc_SRRpars(prod_spawn_em * prod_em_smolt, capacity_em_smolt, fec_spawn, p_female)
      #R0 <- MSEtool::R0conv(SRRpars[1], SRRpars[2], phi0)
      #NOR_total <- R0 * phi0
    } else {
      pbar[, g-1] <- AHA_loop[[g-1]]$pbar
      HOR_total <- AHA_loop[[g-1]]$adult_HOR
      NOR_total <- AHA_loop[[g-1]]$adult_HOS + AHA_loop[[g-1]]$adult_NOS
    }

    AHA_loop[[g]] <- .AHA_loop(
      g,
      HOR_total, NOR_total,
      surv_passage_adult, surv_passage_juv,
      u_HOR, u_NOR,
      p_NOB_goal,
      brood_local, brood_import,
      p_NOR_brood_max,
      p_return_hatchery,
      surv_adult_return_apply[g], SAR_ratio = surv_adult_return_apply[g]/surv_adult_return_geo,
      p_weir_efficiency,
      strays_total,
      RRS_HOS,
      p_HOS_goal,
      fitness_type, Zpop_start, theta, omega2, A, phenotype_variance, fitness_floor, heritability,
      pbar_prev = if (g == 1) NA else pbar[, g-1],
      rel_loss_egg, rel_loss_fry, rel_loss_smolt,
      prod_spawn_em, prod_em_smolt, prod_smolt_adult,
      capacity_spawn_em, capacity_em_smolt, capacity_smolt_adult,
      egg_per_spawner, surv_egg_release, surv_release_adult[g]
    )
  }

  # Reference points in the first entry in the loop
  AHA_loop[[1]]$Ref <- c(Neq = Neq, MSY = MSY, SMSY = SMSY, umsy = umsy)

  return(AHA_loop)
}


.AHA_loop <- function(g = 1, HOR_total = 1e3, NOR_total = 1e3,
                      surv_passage_adult, surv_passage_juv,
                      u_HOR, u_NOR, p_NOB_goal, brood_local, brood_import,
                      p_NOR_brood_max, p_return_hatchery,
                      surv_adult_return_apply,
                      SAR_ratio,
                      p_weir_efficiency,
                      strays_total,
                      RRS_HOS,
                      p_HOS_goal,
                      fitness_type,
                      Zpop_start,
                      theta, omega2, A, phenotype_variance, fitness_floor, heritability, pbar_prev,
                      rel_loss_egg, rel_loss_fry, rel_loss_smolt,
                      prod_spawn_em, prod_em_smolt, prod_smolt_adult,
                      capacity_spawn_em, capacity_em_smolt, capacity_smolt_adult,
                      egg_per_spawner, surv_egg_release, surv_release_adult, ...) {

  # Fishery harvest
  nfishery <- 4
  do_harvest_NOR <- catch_fn(NOR_total, u_NOR, surv_passage = c(1, 1, 1, surv_passage_adult), nfishery = nfishery)
  do_harvest_HOR <- catch_fn(HOR_total, u_HOR, surv_passage = c(1, 1, 1, surv_passage_adult), nfishery = nfishery)

  harvest_NOR <- sum(do_harvest_NOR$catch)
  harvest_HOR <- sum(do_harvest_HOR$catch)

  #AHA is always programmed so that NOR escapement is a minimum of 1.
  #NOB is limited by the max % of the NOR run specified in AHA
  escapement_NOR <- do_harvest_NOR$N[nfishery + 1]
  escapement_HOR <- do_harvest_HOR$N[nfishery + 1]

  brood_NOB <- ifelse(escapement_NOR > 1,
                      min(p_NOB_goal * brood_local, p_NOR_brood_max * (escapement_NOR - 1)),
                      0)

  brood_HOB <- local({
    aa <- min(brood_NOB/p_NOB_goal, brood_local, escapement_HOR)
    bb <- max(brood_local, aa)
    cc <- min(brood_local, escapement_HOR)
    dd <- ifelse(p_NOB_goal > 0, bb - brood_NOB, cc)
    ee <- max(0, dd)

    min(ee, p_return_hatchery * escapement_HOR)
  })

  total_brood <- brood_NOB + brood_HOB + brood_import

  NOS <- escapement_NOR - brood_NOB

  HOS_total <- max(0, escapement_HOR * (1 - p_return_hatchery)) * (1 - p_weir_efficiency) + SAR_ratio * strays_total
  HOS_effective <- max(0, escapement_HOR * (1 - p_return_hatchery)) * (1 - p_weir_efficiency) * RRS_HOS + SAR_ratio * strays_total

  #HOS_goal <- min(HOS_effective, NOS * p_HOS_goal/(1 - p_HOS_goal))
  #HOS_surplus <- HOS_total - HOS_goal

  total_escapement <- HOS_total + NOS
  total_spawners <- HOS_effective + NOS

  pNOB <- ifelse(total_brood > 0, brood_NOB/total_brood,
                 ifelse(brood_import > 0, 0, 1))
  pNOS <- NOS/total_spawners

  # Fitness this generation
  if (fitness_type == "none") {

    fitness <- 1

  } else if (g == 1) {

    if (fitness_type == "Ford") {

      pbar <- Zpop_start[1:2]
      fitness <- local({
        x <- exp(-0.5 * (pbar[1] - theta[1])^2/(omega2 + phenotype_variance))
        max(fitness_floor, x)
      })

    } else { # Busack

      Zpop <- Zpop_start
      fitness <- abs(Zpop[1, g] - mean(theta[-1]))/abs(theta[1] - mean(theta[-1]))

    }

  } else if (fitness_type == "Ford") {

    do_fitness <- calc_Ford_fitness(pNOB, pHOS = 1 - pNOS, pbar_prev,
                                    theta, omega2, phenotype_variance, A, fitness_floor, heritability)
    pbar <- do_fitness$pbar
    fitness <- do_fitness$fitness

  } else { # Busack - incomplete
    #do_fitness <- calc_Busack_fitness()
    #Zpop <- do_fitness$Zpop
    #fitness <- do_fitness$fitness
  }

  # AQ,AR,AS
  fitness_egg <- fitness^rel_loss_egg
  fitness_fry <- fitness^rel_loss_fry
  fitness_smolt <- fitness^rel_loss_smolt

  # AT, AU, AV, AW
  fitprod_spawn_em <- prod_spawn_em * fitness_egg
  fitprod_em_smolt <- prod_em_smolt * fitness_fry
  fitprod_smolt_adult <- surv_passage_juv * surv_adult_return_apply * fitness_smolt
  fitprod_cumulative <- fitprod_spawn_em * fitprod_em_smolt * fitprod_smolt_adult

  # AX, AY, AZ, BA
  fitcap_spawn_em <- capacity_spawn_em * fitness_egg
  fitcap_em_smolt <- capacity_em_smolt * fitness_fry
  fitcap_smolt_adult <- capacity_smolt_adult * fitness_smolt
  fitcap_cumulative <-
    fitprod_cumulative/(fitprod_spawn_em/fitcap_spawn_em+fitprod_spawn_em*
                         fitprod_em_smolt/fitcap_em_smolt+fitprod_spawn_em*
                         fitprod_em_smolt*fitprod_smolt_adult/fitcap_smolt_adult)

  ## Next generation - natural production
  # AI, AJ - offspring from hatchery origin spawners + natural origin spawners
  fry_HOS <- calc_SRR(HOS_effective, total_spawners, fitprod_spawn_em, fitcap_spawn_em)
  fry_NOS <- calc_SRR(NOS, total_spawners, fitprod_spawn_em, fitcap_spawn_em)
  total_fry <- fry_HOS + fry_NOS

  smolt_HOS <- calc_SRR(fry_HOS, total_fry, fitprod_em_smolt, fitcap_em_smolt)
  smolt_NOS <- calc_SRR(fry_NOS, total_fry, fitprod_em_smolt, fitcap_em_smolt)
  total_smolt <- smolt_HOS + smolt_NOS

  SRRpar <- calc_SRRpars(fitprod_em_smolt, fitcap_em_smolt)

  ## Hatchery origin release
  fry_HOR <- total_brood * egg_per_spawner
  smolt_HOR <- total_brood * egg_per_spawner * surv_egg_release

  # Next generation return
  adult_HOS <- calc_SRR(smolt_HOS, total_smolt + smolt_HOR, fitprod_smolt_adult, fitcap_smolt_adult)
  adult_NOS <- calc_SRR(smolt_NOS, total_smolt + smolt_HOR, fitprod_smolt_adult, fitcap_smolt_adult)

  adult_HOR <- local({
    smolt_temp <- smolt_HOR * surv_release_adult
    smolt_temp2 <- smolt_temp + total_smolt * surv_adult_return_apply * surv_passage_juv
    smolt_temp/(1 + smolt_temp2/capacity_smolt_adult)
  })

  out <- list(
    NOR = NOR_total, # Natural origin return (function input)
    HOR = HOR_total, # Hatchery origin return (function input)
    harvest_NOR = harvest_NOR,
    harvest_HOR = harvest_HOR,
    escapement_NOR = escapement_NOR,
    escapement_HOR = escapement_HOR,
    NOB = brood_NOB,
    HOB = brood_HOB,
    NOS = NOS,
    HOS = HOS_total,
    HOS_effective = HOS_effective,
    fitness = fitness,
    fry_NOS = fry_NOS, # Next generation
    fry_HOS = fry_HOS,
    fry_HOR = fry_HOR, # Hatchery origin release
    smolt_NOS = smolt_NOS,
    smolt_HOS = smolt_HOS,
    smolt_HOR = smolt_HOR,
    adult_HOS = adult_HOS, # Actually the return
    adult_NOS = adult_NOS,
    adult_HOR = adult_HOR, # Hatchery origin return
    SAR_loss = fitprod_smolt_adult,
    pNOB = pNOB,
    pHOSeff = 1 - pNOS,
    alpha = SRRpar[1],
    beta = SRRpar[2]
  )

  if (fitness_type == "Ford") out$pbar <- pbar

  return(out)
}
