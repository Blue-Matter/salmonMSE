

#' Internal stray function
#'
#' Calculates the number of fish that move from donor to recipient populations from a movement matrix
#'
#' @param N Array `[ns, nage, n_r]` of hatchery-origin fish
#' @param stray_matrix Matrix `[ns, ns]` that specifies proportion that move from population in the i-th row to the j-th column. The diagonal informs non-straying proportion
#' @param m Vector `[ns]` mark rater of hatchery-origin fish by donor population
#' @returns Named list:
#' - `N_remain` Abundance of hatchery-origin fish that do not stray, same dimension as `N`
#' - `N_stray` Abundance of strays by recipient population, same dimension as `N`
#' - `m_stray` Mark rate of strays by recipient population, vector `[ns]`
#' @keywords internal
stray_func <- function(N, stray_matrix, m) {

  ns <- dim(N)[1]
  N_stray_marked <- N_stray_unmarked <- N_remain <- array(0, dim(N))

  for (s_recipient in 1:ns) {
    for (s_donor in 1:ns) {
      if (s_recipient == s_donor) {
        N_remain[s_donor, , ] <- stray_matrix[s_donor, s_recipient] * N[s_donor, , ]
      } else {
        N_stray_marked[s_recipient, , ] <- N_stray_marked[s_recipient, , ] +
          stray_matrix[s_donor, s_recipient] * N[s_donor, , ] * m[s_donor]

        N_stray_unmarked[s_recipient, , ] <- N_stray_unmarked[s_recipient, , ] +
          stray_matrix[s_donor, s_recipient] * N[s_donor, , ] * (1 - m[s_donor])
      }
    }
  }
  N_stray <- N_stray_marked + N_stray_unmarked
  m_stray <- apply(N_stray_marked, 1, sum)/apply(N_stray, 1, sum)
  m_stray[is.na(m_stray)] <- 0

  list(N_remain = N_remain, N_stray = N_stray, m_stray = m_stray)
}



FW_func <- function(Nage_NOS, Nage_HOS, Nage_stray, m, m_stray, s_enroute,
                    hatchery_args = list(),
                    zbar_brood, fitness_args = list(),
                    fec, p_female, habitat_args,
                    SRRpars = data.frame(), p_LHG) {

  # Calculate broodtake, in-river removals, and spawners arriving at spawning grounds
  Brood_Calcs <- brood_func(
    Nage_NOS,
    Nage_HOS,
    Nage_stray,
    m,
    m_stray,
    s_enroute,
    hatchery_args = hatchery_args
  )

  if (inherits(habitat_args, "Habitat")) {

    NextGen_Calcs <- nextgen_habitat_func(
      Brood_Calcs = Brood_Calcs,
      Habitat = habitat_args,
      fec = fec,
      p_female = p_female,
      hatchery_args = hatchery_args,
      fitness_args = fitness_args,
      zbar_brood = zbar_brood,
      p_LHG = p_LHG
    )

  } else {

    NextGen_Calcs <- nextgen_SRR_func(
      Brood_Calcs = Brood_Calcs,
      fec = fec,
      p_female = p_female,
      hatchery_args = hatchery_args,
      fitness_args = fitness_args,
      zbar_brood = zbar_brood,
      SRRpars = SRRpars,
      p_LHG = p_LHG
    )
  }

  output <- c(
    Brood_Calcs$broodtake[c("pNOB", "NOB", "HOB_unmarked", "HOB_marked", "HOB_import", "HOB_stray")],
    Brood_Calcs$hatchery_production,
    Brood_Calcs$spawners[c("NO_remove", "HO_remove")],
    NextGen_Calcs
  )

  return(output)
}


#' Brood function
#'
#' Calculates broodtake and in-river removals from escapement of marine fisheries. This function also applies en-route mortality.
#'
#' @param Nage_NOS Array `[nage, n_g]` of natural-origin fish
#' @param Nage_HOS Array `[nage, n_r]` of hatchery-origin fish
#' @param stray_external Array `[nage, n_r]` of hatchery-origin strays
#' @param m Numeric, mark rate of `Nage_HOS`
#' @param m_stray Numeric, mark rate of `stray_external`
#' @param s_enroute Numeric, en-route survival
#' @param hatchery_args List of various hatchery arguments. See details.
#' @details
#' `hatchery_args` is a named list containing:
#' -
#' -
#' @returns Named list:
#' - `broodtake` list returned by [calc_broodtake()]
#' - `hatchery_production` list returned by [calc_yearling()]
#' - `spawners` list returned by [calc_spawners()]
#' @keywords internal
#' @seealso [calc_broodtake()]
brood_func <- function(Nage_NOS, Nage_HOS, stray_external, m, m_stray, s_enroute, hatchery_args) {

  # En-route mortality
  Nage_NOS_enroute <- Nage_NOS * s_enroute
  Nage_HOS_enroute <- Nage_HOS * s_enroute
  stray_external_enroute <- stray_external * s_enroute

  # Broodtake
  custom_brood_rule <- is.function(hatchery_args$f_brood) && !is.null(formals(hatchery_args$f_brood))
  if (!custom_brood_rule || !hatchery_args$egg_target) { # AHA approach
    if (hatchery_args$egg_target > 0) {
      Nage_NOS_avail_brood <- Nage_NOS_enroute * hatchery_args$pmax_esc
      Nage_HOS_avail_brood <- Nage_HOS_enroute * hatchery_args$pmax_esc
    } else {
      Nage_NOS_avail_brood <- Nage_NOS_enroute
      Nage_HOS_avail_brood <- Nage_HOS_enroute
    }

    broodtake <- calc_broodtake(
      NO = Nage_NOS_avail_brood,
      HO = Nage_HOS_avail_brood,
      stray_external_enroute,
      hatchery_args$brood_import,
      hatchery_args$ptarget_NOB,
      hatchery_args$pmax_NOB,
      hatchery_args$phatchery,
      hatchery_args$egg_target,
      hatchery_args$p_female,
      hatchery_args$fec_brood,
      hatchery_args$s_prespawn,
      m
    )
  } else { # Custom brood approach
    broodtake <- calc_broodtake_custom(
      f_brood = hatchery_args$f_brood,
      NO = Nage_NOS_enroute,
      HO = Nage_HOS_enroute,
      stray_external_enroute,
      hatchery_args$p_female,
      hatchery_args$fec_brood,
      hatchery_args$s_prespawn,
      m,
      hatchery_args$egg_target
    )
  }

  hatchery_production <- calc_yearling(
    sum(broodtake$egg_NOB, broodtake$egg_HOB_unmarked, broodtake$egg_HOB_marked, broodtake$egg_HOB_import),
    hatchery_args$s_yearling,
    hatchery_args$s_subyearling,
    hatchery_args$p_yearling,
    hatchery_args$p_subyearling
  )

  # Spawners arriving at spawning grounds and in-river removals after broodtake
  spawners <- calc_spawners(
    broodtake, Nage_NOS_enroute, Nage_HOS_enroute, stray_external_enroute,
    ifelse(custom_brood_rule, NA, hatchery_args$phatchery),
    hatchery_args$premove_HOS, hatchery_args$premove_NOS, hatchery_args$m
  )

  output <- list(
    broodtake = broodtake,
    hatchery_production = hatchery_production,
    spawners = spawners
  )

  return(output)
}

#' Spawning and early life stage function
#'
#' Calculates egg production from spawners arriving at spawning grounds, broodtake and in-river removals from escapement of marine fisheries. This function also applies en-route mortality.
#'
#' @param Brood_Calcs List, returned by [brood_func()]
#' @param fec Vector `[nage]`, fecundity at age of female fish
#' @param p_female Vector `[nage]`, proportion female by age class
#' @param fitness_args List
#' @returns Named list:
#' - `NOS` Matrix `[nage, n_g]`
#' - `HOS` Matrix `[nage, n_r]`
#' - `HOS_effective` Matrix `[nage, n_r]`
#' - `HOS_stray` Matrix `[nage, n_r]`
#' - `pHOSeff` Numeric
#' - `pHOScensus` Numeric
#' - `Egg_NOS` `[nage, n_g]`
#' - `Egg_HOS` `[nage, n_r]`
#' - `Smolt_RelOut` Vector `[n_r]`
#' - `Fry_NOS` Vector `[n_g]`
#' - `Fry_HOS` Vector `[n_g]`
#' - `Smolt_NOS` Vector `[n_g]`
#' - `Smolt_HOS` Vector `[n_g]`
#' @keywords internal
nextgen_SRR_func <- function(Brood_Calcs, fec, p_female,
                             hatchery_args = list(), fitness_args = list(),
                             zbar_brood, SRRpars, p_LHG) {

  NOS <- Brood_Calcs$spawners$NOS
  HOS <- Brood_Calcs$spawners$HOS + Brood_Calcs$spawners$HOS_stray
  if (!is.null(hatchery_args$gamma)) {
    HOS_effective <- HOS * hatchery_args$gamma
  } else {
    HOS_effective <- HOS
  }

  #### Spawners weighted by fecundity ----
  pHOSeff <- sum(fec * HOS_effective)/sum(fec * NOS, fec * HOS_effective)
  pHOScensus <- sum(fec * HOS)/sum(fec * NOS, fec * HOS)

  #### Natural egg production
  # maxage x life cycle group (NOS) and release strategy (HOS)
  Egg_NOS <- NOS * p_female * fec
  Egg_HOS <- HOS_effective * p_female * fec

  total_egg <- sum(Egg_NOS, Egg_HOS)

  #### Calculate fitness ----
  Egg_NOB <- Brood_Calcs$broodtake$egg_NOB
  Egg_HOB <- Brood_Calcs$broodtake$egg_HOB_unmarked + Brood_Calcs$broodtake$egg_HOB_marked
  Egg_HOB[, 1] <- Egg_HOB[, 1] + Brood_Calcs$broodtake$egg_HOB_import
  fitness_calcs <- fitness_func(
    Egg_NOS = rowSums(Egg_NOS),
    Egg_HOS = rowSums(Egg_HOS),
    Egg_NOB = rowSums(Egg_NOB),
    Egg_HOB = rowSums(Egg_HOB),
    zbar_brood = zbar_brood,
    fitness_args
  )
  fitness_loss <- fitness_calcs$fitness_loss

  #### Egg-fry survival ----
  Fry_NOS <- Egg_NOS * fitness_loss[1, 1]
  Fry_HOS <- Egg_HOS * fitness_loss[1, 1]

  #### Fry-smolt survival ----
  # This can be density-dependent in competition with hatchery releases
  # Hatchery releases (after fitness loss) by release strategy (next generation)
  yearling <- Brood_Calcs$hatchery_production$yearling * fitness_loss[2, 1] * fitness_loss[2, 2]
  subyearling <- Brood_Calcs$hatchery_production$subyearling * fitness_loss[2, 1]

  # Juveniles (ostensibly fry) in freshwater environment that experience density-dependent survival together
  total_fry_DD <- sum(Fry_NOS, Fry_HOS)
  if (sum(yearling) && hatchery_args$yearling_DD) total_fry_DD <- sum(total_fry_DD, yearling)
  if (sum(subyearling) && hatchery_args$subyearling_DD) total_fry_DD <- sum(total_fry_DD, subyearling)

  # Calculate hatchery-origin density-dependent survival
  if (sum(yearling) && hatchery_args$yearling_DD) {
    smolt_yearling <- calc_smolt(
      yearling, total_fry_DD,
      SRRpars[["kappa"]], SRRpars[["capacity"]], SRRpars[["Smax"]], SRRpars[["phi"]], SRRpars[["tau"]],
      fitness_loss[2, 1] * fitness_loss[2, 2],
      SRRpars[["SRrel"]]
    )
  } else {
    smolt_yearling <- yearling
  }

  if (sum(subyearling) && hatchery_args$subyearling_DD) {
    smolt_subyearling <- calc_smolt(
      subyearling, total_fry_DD,
      SRRpars[["kappa"]], SRRpars[["capacity"]], SRRpars[["Smax"]], SRRpars[["phi"]], SRRpars[["tau"]],
      fitness_loss[2, 1] * fitness_loss[2, 2],
      SRRpars[["SRrel"]]
    )
  } else {
    smolt_subyearling <- subyearling
  }

  # Total outmigrating hatchery releases
  Smolt_RelOut <- smolt_yearling + smolt_subyearling

  # Natural-origin fry production: Calculate total number and re-distribute among life history group in next generation
  Fry_NOS_g <- sum(Fry_NOS) * p_LHG
  Fry_HOS_g <- sum(Fry_HOS) * p_LHG

  # Natural smolt production for life history group g (next generation)
  Smolt_NOS <- calc_smolt(
    Fry_NOS_g, total_fry_DD,
    SRRpars[["kappa"]], SRRpars[["capacity"]], SRRpars[["Smax"]], SRRpars[["phi"]], SRRpars[["tau"]],
    fitness_loss[1, 2],
    SRRpars[["SRrel"]]
  )
  Smolt_HOS <- calc_smolt(
    Fry_HOS_g, total_fry_DD,
    SRRpars[["kappa"]], SRRpars[["capacity"]], SRRpars[["Smax"]], SRRpars[["phi"]], SRRpars[["tau"]],
    fitness_loss[1, 2],
    SRRpars[["SRrel"]]
  )

  output <- list(
    NOS = NOS,
    HOS = HOS,
    HOS_effective = HOS_effective,
    HOS_stray = Brood_Calcs$spawners$HOS_stray,
    pHOSeff = pHOSeff,
    pHOScensus = pHOScensus,
    Egg_NOS = Egg_NOS,
    Egg_HOS = Egg_HOS,
    Smolt_RelOut = Smolt_RelOut,
    Fry_NOS = Fry_NOS_g,
    Fry_HOS = Fry_HOS_g,
    Smolt_NOS = Smolt_NOS,
    Smolt_HOS = Smolt_HOS
  )

  return(c(output, fitness_calcs))
}

#' @name nextgen_SRR_func
#' @param Habitat \linkS4class{Habitat} object, modified by [ProjectSOM()]
#' @keywords internal
nextgen_habitat_func <- function(Brood_Calcs, Habitat, fec, p_female,
                                 hatchery_args = list(), fitness_args = list(),
                                 zbar_brood, p_LHG) {

  NOS_init <- Brood_Calcs$spawners$NOS
  HOS_init <- Brood_Calcs$spawners$HOS
  HOS_stray_init <- Brood_Calcs$spawners$HOS_stray

  #### Pre-spawn mortality ----
  # By parental life cycle group x age
  NOS <- calc_SRR(
    NOS_init, sum(NOS_init, HOS_init, HOS_stray_init),
    p = Habitat@prespawn_prod, capacity = Habitat@prespawn_capacity,
    type = Habitat@prespawn_rel
  )

  # By parental release strategy x age
  HOS_ <- calc_SRR(
    HOS_init, sum(NOS_init, HOS_init, HOS_stray_init),
    p = Habitat@prespawn_prod, capacity = Habitat@prespawn_capacity,
    type = Habitat@prespawn_rel
  )
  HOS_stray <- calc_SRR(
    HOS_stray_init, sum(NOS_init, HOS_init, HOS_stray_init),
    p = Habitat@prespawn_prod, capacity = Habitat@prespawn_capacity,
    type = Habitat@prespawn_rel
  )
  HOS <- HOS_ + HOS_stray
  if (!is.null(hatchery_args$gamma)) {
    HOS_effective <- HOS * hatchery_args$gamma
  } else {
    HOS_effective <- HOS
  }

  # Spawners weighted by fecundity
  pHOSeff <- sum(fec * HOS_effective)/sum(fec * NOS, fec * HOS_effective)
  pHOScensus <- sum(fec * HOS)/sum(fec * NOS, fec * HOS)

  #### Egg production by parental age ----
  # Recorded as the egg production in the population
  # Potential future function for survival that is a function of parental age
  Egg_NOS <- NOS * p_female * fec
  Egg_HOS <- HOS_effective * p_female * fec

  total_egg <- sum(Egg_NOS, Egg_HOS)

  #### Calculate fitness ----
  Egg_NOB <- Brood_Calcs$broodtake$egg_NOB
  Egg_HOB <- Brood_Calcs$broodtake$egg_HOB_unmarked + Brood_Calcs$broodtake$egg_HOB_marked
  Egg_HOB[, 1] <- Egg_HOB[, 1] + Brood_Calcs$broodtake$egg_HOB_import
  fitness_calcs <- fitness_func(
    Egg_NOS = rowSums(Egg_NOS),
    Egg_HOS = rowSums(Egg_HOS),
    Egg_NOB = rowSums(Egg_NOB),
    Egg_HOB = rowSums(Egg_HOB),
    zbar_brood = zbar_brood,
    fitness_args
  )
  fitness_loss <- fitness_calcs$fitness_loss

  #### Egg Incubation survival ----
  # By parental life cycle group and age
  Egg_inc_NOS <- calc_SRR(
    Egg_NOS, sum(Egg_HOS, Egg_NOS),
    p = Habitat@egg_prod, capacity = Habitat@egg_capacity,
    type = Habitat@egg_rel
  )
  # By parental release strategy and age
  Egg_inc_HOS <- calc_SRR(
    Egg_HOS, sum(Egg_HOS, Egg_NOS),
    p = Habitat@egg_prod, capacity = Habitat@egg_capacity,
    type = Habitat@egg_rel
  )

  #### Egg-fry survival ----
  # By parental life cycle group
  fry_prod <- Habitat@fry_prod * fitness_loss[1, 1]
  fry_cap <- Habitat@fry_capacity * fitness_loss[1, 1]
  Fry_NOS <- calc_SRR(
    rowSums(Egg_inc_NOS), sum(Egg_inc_NOS, Egg_inc_HOS),
    p = fry_prod, capacity = fry_cap,
    type = Habitat@fry_rel
  ) * Habitat@fry_sdev[1, 1]

  # By parental release strategy
  Fry_HOS <- calc_SRR(
    rowSums(Egg_inc_HOS), sum(Egg_inc_NOS, Egg_inc_HOS),
    p = fry_prod, capacity = fry_cap,
    type = Habitat@fry_rel
  ) * Habitat@fry_sdev[1, 1]

  #### Fry-smolt survival ----
  # This can be density-dependence in competition with hatchery releases
  # Hatchery releases (after fitness loss) by release strategy (next generation)
  yearling <- Brood_Calcs$hatchery_production$yearling * fitness_loss[2, 1] * fitness_loss[2, 2]
  subyearling <- Brood_Calcs$hatchery_production$subyearling * fitness_loss[2, 1]

  # Juveniles (ostensibly fry) in freshwater environment that experience density-dependent survival together
  total_fry_DD <- sum(Fry_NOS, Fry_HOS)
  if (sum(yearling) && hatchery_args$yearling_DD) total_fry_DD <- sum(total_fry_DD, yearling)
  if (sum(subyearling) && hatchery_args$subyearling_DD) total_fry_DD <- sum(total_fry_DD, subyearling)

  # Calculate hatchery-origin density-dependent survival
  smolt_HO_prod <- Habitat@smolt_prod * fitness_loss[2, 2]
  smolt_HO_cap <- Habitat@smolt_capacity * fitness_loss[2, 2]
  if (hatchery_args$yearling_DD) {
    yearling_out <- calc_SRR(
      yearling, total_fry_DD,
      p = smolt_HO_prod, capacity = smolt_HO_cap,
      type = Habitat@smolt_rel
    )
  } else {
    yearling_out <- yearling
  }
  if (hatchery_args$subyearling_DD) {
    subyearling_out <- calc_SRR(
      subyearling, total_fry_DD,
      p = smolt_HO_prod, capacity = smolt_HO_cap,
      type = Habitat@smolt_rel
    )
  } else {
    subyearling_out <- subyearling
  }

  # Total outmigrating hatchery releases
  Smolt_RelOut <- yearling_out + subyearling_out

  # Natural-origin fry production: Calculate total number and re-distribute among life history group in next generation
  Fry_NOS_g <- sum(Fry_NOS) * p_LHG
  Fry_HOS_g <- sum(Fry_HOS) * p_LHG

  # Natural smolt production for life history group g (next generation)
  smolt_NO_prod <- Habitat@smolt_prod * fitness_loss[1, 2]
  smolt_NO_cap <- Habitat@smolt_capacity * fitness_loss[1, 2]
  Smolt_NOS <- calc_SRR(
    Fry_NOS_g, total_fry_DD,
    p = smolt_NO_prod, capacity = smolt_NO_cap,
    type = Habitat@smolt_rel
  ) * Habitat@smolt_sdev[1, 1]

  Smolt_HOS <- calc_SRR(
    Fry_HOS_g, total_fry_DD,
    p = smolt_NO_prod, capacity = smolt_NO_cap,
    type = Habitat@smolt_rel
  ) * Habitat@smolt_sdev[1, 1]

  output <- list(
    NOS = NOS,
    HOS = HOS,
    HOS_effective = HOS_effective,
    HOS_stray = HOS_stray,
    pHOSeff = pHOSeff,
    pHOScensus = pHOScensus,
    Egg_NOS = Egg_NOS,
    Egg_HOS = Egg_HOS,
    Smolt_RelOut = Smolt_RelOut,
    Fry_NOS = Fry_NOS_g,
    Fry_HOS = Fry_HOS_g,
    Smolt_NOS = Smolt_NOS,
    Smolt_HOS = Smolt_HOS
  )

  return(c(output, fitness_calcs))
}
