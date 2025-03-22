
#' Predict smolt production
#'
#' @description Internal function that predicts the natural origin or hatchery origin smolt production from the
#' escapement, and saves state variables to [salmonMSE::salmonMSE_env].
#' - `smolt_func()` is the population dynamics function
#' - `makeRel_smolt()` generates a list for openMSE to use in the simulations
#'
#' @param Nage_NOS Matrix `[n_age, n_g]` of natural escapement from openMSE. `n_g` is the number of life history groups (divisons within cohorts
#' that contribute to spawning.
#' @param Nage_HOS Matrix `[n_age, n_g]` of hatchery escapement from openMSE. `n_g` is the number of life history groups (divisons within cohorts
#' that contribute to spawning.
#' @param x Integer, simulation number from openMSE
#' @param y Integer, simulation year (including historical years)
#' @param output Character, whether to predic the natural origin or hatchery origin smolt production
#' @param s_enroute Numeric, en route survival of the escapement to the spawning grounds
#' @param p_female Numeric, proportion female for calculating the egg production.
#' @param fec Vector `[maxage]`. The fecundity at age schedule of spawners
#' @param s_egg_fry Numeric. Egg-fry survival
#' @param SRRpars Data frame containing stock recruit parameters for natural smolt production.
#' Column names include: SRrel, kappa, capacity_smolt, Smax, phi, kappa_improve, capacity_smolt_improve
#' @param hatchery_args Named list containing various arguments controlling broodtake and hatchery production. See details below.
#' @param fitness_args Named list containing various arguments controlling population fitness from hatchery production.
#' Names include: fitness_type, omega2, theta, fitness_variance, heritability, zbar_start, fitness_floor, rel_loss
#' @param s Integer, salmonMSE population index. Used to report variables to [salmonMSE::salmonMSE_env].
#' @section hatchery_args:
#' Hatchery control parameters are included in a named list with the following arguments:
#'
#' - `egg_target` Numeric, target egg production for hatchery. Set to zero for no hatchery production.
#' - `ptarget_NOB` Numeric, the target proportion of the natural origin broodtake relative to the overall broodtake
#' - `pmax_NOB` Numeric, the maximum proportion of the natural origin escapement to be used as broodtake
#' - `fec_brood` Vector `[maxage]`. the fecundity at age schedule of broodtake to calculate the total hatchery egg production
#' - `s_yearling` Numeric, the survival of eggs to the smolt life stage (for yearling release)
#' - `s_subyearling` Numeric. the survival of eggs to subyearling life stage (for subyearling release)
#' - `p_yearling` Numeric, the proportion of annual releases at the yearling life stage (vs. subyearling)
#' - `phatchery` Numeric, the proportion of the hatchery origin escapement that return to the hatchery,
#' for example, by removal from spawning grounds or swim-in facilities. These fish are available for broodtake.
#' - `premove_HOS` Numeric, the proportion of the hatchery origin escapement to be removed from the spawning
#' grounds (in order to ensure a high proportion of NOS). These fish are not available for broodtake.
#' For example, a value less than one can represent imperfect implementation of weir removal.
#' - `s_prespawn` Numeric, the survival of broodtake prior to egg production. `1 - s_prespawn` is the proportion of fish
#' not used for hatchery purposes, e.g., mortality or other resesarch purposes. Used to back-calculate the broodtake.
#' - `gamma` Numeric, the relative reproductive success of hatchery origin spawners (relative to natural origin spawners)
#' - `m` Numeric, mark rate for selective broodtake
#' @returns
#' - `smolt_func()` returns a numeric for the ratio of the realized smolt production vs. the hypothetical value if there were no
#' hatchery, en route mortality, or habitat improvement
#' - `makeRel_smolt()` returns a list that is passed to openMSE as a inter-population relationship
#' @keywords internal
smolt_func <- function(Nage_NOS, Nage_HOS, x = -1, y, output = c("natural", "hatchery"),
                       s_enroute, p_female, fec, s_egg_fry, SRRpars, # Spawning (natural production)
                       hatchery_args, fitness_args, s, g, prop_LHG, r) {

  output <- match.arg(output)
  Nage_NOS[is.na(Nage_NOS)] <- 0
  Nage_NOS_enroute <- Nage_NOS * s_enroute
  stray_external_enroute <- hatchery_args$stray_external * s_enroute

  # Hatchery
  if (missing(Nage_HOS)) Nage_HOS <- matrix(0, nrow(Nage_NOS), 1)
  Nage_HOS[is.na(Nage_HOS)] <- 0
  Nage_HOS_enroute <- Nage_HOS * s_enroute

  if (hatchery_args$egg_target > 0 && sum(Nage_NOS, Nage_HOS)) {

    Nage_NOS_avail_brood <- Nage_NOS_enroute * hatchery_args$pmax_esc
    Nage_HOS_avail_brood <- Nage_HOS_enroute * hatchery_args$pmax_esc
    broodtake <- calc_broodtake(
      NOR_escapement = Nage_NOS_avail_brood,
      HOR_escapement = Nage_HOS_avail_brood,
      stray_external_enroute,
      hatchery_args$brood_import,
      hatchery_args$ptarget_NOB,
      hatchery_args$pmax_NOB,
      hatchery_args$phatchery,
      hatchery_args$egg_target,
      hatchery_args$p_female,
      hatchery_args$fec_brood,
      hatchery_args$s_prespawn,
      hatchery_args$m
    )

    egg_NOB <- broodtake$egg_NOB
    egg_HOB <- broodtake$egg_HOB_unmarked + broodtake$egg_HOB_marked

    pNOB <- broodtake$pNOB

    hatchery_production <- calc_yearling(
      egg_NOB + egg_HOB,
      hatchery_args$s_yearling,
      hatchery_args$s_subyearling,
      hatchery_args$p_yearling,
      hatchery_args$p_subyearling
    )

  } else {

    # Check for same format as list returned by .broodtake_func
    broodtake <- list(
      egg_NOB = 0,
      egg_HOB_unmarked = 0,
      egg_HOB_marked = 0,
      ptake_unmarked = 0,
      ptake_marked = 0,
      pNOB = 0,
      NOB = array(0, dim(Nage_NOS)),
      HOB_unmarked = array(0, dim(Nage_HOS)),
      HOB_marked = array(0, dim(Nage_HOS)),
      HOB_import = rep(0, nrow(Nage_NOS)),
      HOB_stray = array(0, dim(Nage_HOS))
    )
    egg_NOB <- egg_HOB <- pNOB <- 0

    hatchery_production <- list(
      yearling = rep(0, ncol(Nage_HOS)),
      subyearling = rep(0, ncol(Nage_HOS))
    )

  }

  # Spawners
  if (sum(Nage_NOS, Nage_HOS)) {
    spawners <- calc_spawners(
      broodtake, Nage_NOS_enroute, Nage_HOS_enroute, stray_external_enroute,
      hatchery_args$phatchery, hatchery_args$premove_HOS, hatchery_args$m
    )
  } else {
    spawners <- list(
      NOS = array(0, dim(Nage_NOS)),
      HOS = array(0, dim(Nage_HOS))
    )
  }

  NOS <- spawners$NOS
  HOS <- spawners$HOS
  if (sum(HOS)) {
    HOS_effective <- spawners$HOS * hatchery_args$gamma
  } else {
    HOS_effective <- array(0, dim(HOS))
  }

  # Spawners weighted by fecundity
  pHOSeff <- sum(fec * HOS_effective)/sum(fec * NOS, fec * HOS_effective)
  pHOScensus <- sum(fec * HOS)/sum(fec * NOS, fec * HOS)

  # Natural egg production in the absence of fitness effects
  Egg_NOS <- colSums(NOS * p_female * fec)
  Egg_HOS <- colSums(HOS_effective * p_female * fec)

  total_egg <- sum(Egg_NOS, Egg_HOS)
  if (output == "natural" && !total_egg) return(0) # Perr_y = 0
  if (output == "hatchery" && !sum(unlist(hatchery_production))) return(0) # Perr_y = 0

  # Fitness
  if (x > 0 && any(fitness_args$fitness_type == "Ford")) {

    # Get zbar from salmonMSE_env
    zbar_prev <- filter(salmonMSE_env$Ford, .data$x == .env$x, .data$s == .env$s)

    if (nrow(zbar_prev)) {
      maxage <- nrow(NOS)
      zbar_brood <- matrix(0, maxage, 2) # Column 1 = natural environment, 2 = hatchery environment

      # Trait value by brood year
      for (a in 1:maxage) {
        if (sum(NOS[a, ]) || sum(broodtake$NOB[a, ])) {
          zbar1 <- dplyr::filter(zbar_prev, .data$t == .env$y - 2 * .env$a, .data$type == "natural") %>%
            pull(.data$zbar)
          if (!length(zbar1)) stop("Cannot find zbar") # zbar1 <- fitness_args$zbar_start[1]
          zbar_brood[a, 1] <- zbar1
        }
        if (sum(HOS_effective[a, ]) || sum(broodtake$HOB[a, ]) || sum(broodtake$HOB_import[a])) {
          zbar2 <- dplyr::filter(zbar_prev, .data$t == .env$y - 2 * .env$a, .data$type == "hatchery") %>%
            pull(.data$zbar)
          if (!length(zbar2)) stop("Cannot find zbar") # zbar2 <- fitness_args$zbar_start[2]
          zbar_brood[a, 2] <- zbar2
        }
      }

      zbar <- calc_zbar(
        rowSums(NOS), rowSums(HOS_effective), rowSums(broodtake$NOB),
        rowSums(broodtake$HOB_marked + broodtake$HOB_unmarked) + broodtake$HOB_import,
        fec, hatchery_args$fec_brood, zbar_brood,
        fitness_args$omega2, fitness_args$theta, fitness_args$fitness_variance, fitness_args$heritability
      )
    } else {
      zbar <- fitness_args$zbar_start
    }

    # Fitness in the natural and hatchery environments
    fitness <- rep(1, 2)
    for (i in 1:2) {
      if (fitness_args$fitness_type[i] == "Ford") {
        fitness[i] <- calc_fitness(
          zbar[i], fitness_args$theta[i], fitness_args$omega2,
          fitness_args$fitness_variance, fitness_args$fitness_floor
        )
      }
    }
    fitness_loss <- outer(fitness, fitness_args$rel_loss, "^")
  } else {
    zbar <- rep(100, 2)
    fitness <- c(1, 1)
    fitness_loss <- matrix(1, 2, 3)
  }

  # Fry production after fitness loss
  Fry_NOS <- Egg_NOS * s_egg_fry * fitness_loss[1, 1]
  Fry_HOS <- Egg_HOS * s_egg_fry * fitness_loss[1, 1]

  # Hatchery production after fitness loss
  yearling <- hatchery_production$yearling * fitness_loss[2, 1] * fitness_loss[2, 2]
  subyearling <- hatchery_production$subyearling * fitness_loss[2, 1]
  total_fry <- sum(Fry_NOS, Fry_HOS, subyearling)

  # Smolt production (natural, but in competition with hatchery subyearlings) after fitness loss
  xx <- max(x, 1)
  SRrel <- match.arg(SRRpars[xx, "SRrel"], choices = c("BH", "Ricker"))

  # Return total hatchery smolt releases after density-dependent survival of subyearlings for release strategy r
  if (output == "hatchery") {
    smolt_subyearling_r <- calc_smolt(
      subyearling[r], total_fry,
      SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
      SRRpars[xx, "kappa_improve"], SRRpars[xx, "capacity_smolt_improve"], fitness_loss[2, 2],
      SRrel
    )
    smolt_rel_r <- yearling[r] + smolt_subyearling_r

    # Save state variable for hatchery production
    if (x > 0) {
      df_H <- data.frame(
        x = x,
        s = s,
        r = r,
        t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
        a = 1:nrow(Nage_HOS),
        Esc_HOS = Nage_HOS[, r],
        HOB = broodtake$HOB_unmarked[, r] + broodtake$HOB_marked[, r] + broodtake$HOB_stray[, r],
        HOS = HOS[, r],
        HOS_effective = HOS_effective[, r]
      )
      df_H$HOB_import <- if (r == 1) broodtake$HOB_import else rep(NA, nrow(Nage_HOS))
      salmonMSE_env$H <- rbind(salmonMSE_env$H, df_H)

      df_stateH <- data.frame(
        x = x,
        s = s,
        r = r,
        t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
        Egg_HOS = Egg_HOS[r], # Spawning output by RS r of this generation
        yearling = yearling[r],
        subyearling = subyearling[r],
        smolt_rel = smolt_rel_r
      )
      salmonMSE_env$stateH <- rbind(salmonMSE_env$stateH, df_stateH)
    }

    return(smolt_rel_r)
  }

  # Fry production for life history group g
  Fry_NOS_g <- sum(Fry_NOS) * prop_LHG
  Fry_HOS_g <- sum(Fry_HOS) * prop_LHG

  smolt_NOS_g <- calc_smolt(
    Fry_NOS_g, total_fry,
    SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRRpars[xx, "kappa_improve"], SRRpars[xx, "capacity_smolt_improve"], fitness_loss[1, 2],
    SRrel
  )

  smolt_HOS_g <- calc_smolt(
    Fry_HOS_g, total_fry,
    SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRRpars[xx, "kappa_improve"], SRRpars[xx, "capacity_smolt_improve"], fitness_loss[1, 2],
    SRrel
  )

  # Predicted smolts from historical SRR parameters and openMSE setup
  # if there were no hatchery production, habitat improvement, enroute mortality, s_egg_fry = 1, or multiple LHG
  Egg_openMSE <- sum(Nage_NOS[, g] * p_female * fec)
  smolt_NOS_SRR <- calc_smolt(
    Egg_openMSE, Egg_openMSE, SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRrel = SRrel
  )

  # MICE parameter to be updated in operating model
  Perr_y <- (smolt_NOS_g + smolt_HOS_g)/smolt_NOS_SRR

  if (x > 0) {
    # Save state variables at age
    df_N <- data.frame(
      x = x,
      s = s,
      g = g,
      t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
      a = 1:nrow(Nage_NOS),
      Esc_NOS = Nage_NOS[, g],
      NOB = broodtake$NOB[, g],
      NOS = NOS[, g]
    )
    salmonMSE_env$N <- rbind(salmonMSE_env$N, df_N)

    df_stateN <- data.frame(
      x = x,
      s = s,
      g = g,
      t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
      Egg_NOS = Egg_NOS[g], # Egg production by LHG g of this generation
      Fry_NOS = Fry_NOS_g,  # Fry production assigned to LHG g for next generation
      Fry_HOS = Fry_HOS_g,
      smolt_NOS = smolt_NOS_g,
      smolt_HOS = smolt_HOS_g,
      pNOB = as.numeric(pNOB),
      pHOSeff = as.numeric(pHOSeff),
      pHOScensus = as.numeric(pHOScensus),
      Perr_y = Perr_y
    )
    salmonMSE_env$stateN <- rbind(salmonMSE_env$stateN, df_stateN)

    # Save zbar
    if (g == 1) {
      df_Ford <- data.frame(
        x = x,
        s = s,
        t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
        type = c("natural", "hatchery"),
        zbar = zbar,
        fitness = fitness
      )
      salmonMSE_env$Ford <- rbind(salmonMSE_env$Ford, df_Ford)

      # Report hatchery state variables for external strays if there is no hatchery production
      if (!sum(unlist(hatchery_production)) && !is.null(hatchery_args$stray_external)) {
        nr <- ncol(hatchery_args$stray_external)

        df_H <- lapply(1:nr, function(r) {
          d <- data.frame(
            x = x,
            s = s,
            r = r,
            t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
            a = 1:nrow(Nage_HOS),
            Esc_HOS = Nage_HOS[, r],
            HOB = broodtake$HOB_unmarked[, r] + broodtake$HOB_marked[, r] + broodtake$HOB_stray[, r],
            HOS = HOS[, r],
            HOS_effective = HOS_effective[, r]
          )
          d$HOB_import <- if (r == 1) broodtake$HOB_import else rep(NA, nrow(Nage_HOS))
          return(d)
        }) %>%
          bind_rows()

        salmonMSE_env$H <- rbind(salmonMSE_env$H, df_H)

        df_stateH <- lapply(1:nr, function(r) {
          data.frame(
            x = x,
            s = s,
            r = r,
            t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
            Egg_HOS = Egg_HOS[r], # Spawning output by RS r of this generation
            yearling = yearling[r],
            subyearling = subyearling[r],
            smolt_rel = 0
          )
        }) %>%
          bind_rows()
        salmonMSE_env$stateH <- rbind(salmonMSE_env$stateH, df_stateH)
      }

    }
  }

  return(Perr_y)
}

#' @rdname smolt_func
#' @param p_smolt Integer, the population index for the smolt production in the openMSE model, corresponding to `output`
#' @param p_natural Integer vector, the population index for the natural origin escapement in the openMSE model. Can be more than one
#' if spawning is from multiple life history groups
#' @param p_hatchery Integer vector, the population index for the hatchery origin escapement in the openMSE model. Can be more than one
#' if there are multiple release strategies. Set to `NULL` for no hatchery production
#' @param g Integer for the life history group of natural origin fish to pass the parameter back to openMSE (if `output = "natural"`)
#' @param prop_LHG Numeric, proportion of the egg production assign to life history group `g` corresponding to `p_smolt` for the next generation (only used if `output = "natural"`)
#' @param r Integer for the release strategy of hatchery origin fish to pass the parameter back to openMSE (if `output = "hatchery"`)
makeRel_smolt <- function(p_smolt = 1, s = 1, p_natural, p_hatchery = NULL,
                          output = c("natural", "hatchery"), s_enroute,
                          p_female, fec, s_egg_fry, SRRpars,  # Spawning (natural production)
                          hatchery_args, fitness_args, g, prop_LHG, r) {

  output <- match.arg(output)

  .smolt_func <- smolt_func

  formals(.smolt_func)$output <- output
  formals(.smolt_func)$s_enroute <- s_enroute
  formals(.smolt_func)$p_female <- p_female
  formals(.smolt_func)$fec <- fec
  formals(.smolt_func)$s_egg_fry <- s_egg_fry

  formals(.smolt_func)$SRRpars <- SRRpars

  formals(.smolt_func)$hatchery_args <- hatchery_args
  formals(.smolt_func)$fitness_args <- fitness_args
  formals(.smolt_func)$s <- s

  if (output == "natural") {
    formals(.smolt_func)$g <- g
    formals(.smolt_func)$prop_LHG <- prop_LHG
  }
  if (output == "hatchery") formals(.smolt_func)$r <- r

  maxage <- length(fec)
  N_natural <- matrix(1, maxage, length(p_natural))

  if (!is.null(p_hatchery)) {
    N_hatchery <- matrix(0, maxage, length(p_hatchery))

    Perr_y <- .smolt_func(N_natural, N_hatchery, x = -1, y = 1)
    model <- c(Perr_y = Perr_y, colSums(N_natural), colSums(N_hatchery), x = -1, y = 1)
    input <- paste0("Nage_", c(p_natural, p_hatchery))

    natural_origin <- c(rep(TRUE, length(p_natural)), rep(FALSE, length(p_hatchery)))
  } else {
    Perr_y <- .smolt_func(N_natural, x = -1, y = 1)
    model <- c(Perr_y = Perr_y, colSums(N_natural), x = -1, y = 1)
    input <- paste0("Nage_", p_natural)

    natural_origin <- rep(TRUE, length(p_natural))
  }

  response <- paste0("Perr_y_", p_smolt)
  terms <- c(response, input, "x", "y")

  out <- list(
    func = .smolt_func,
    model = structure(model, names = terms),
    fitted.values = Perr_y,
    CV = 0,
    terms = terms,
    lag = "next",
    type = "Recruitment deviation",
    Rel = switch(output,
                 "natural" = "Smolt natural production from escapement",
                 "hatchery" = "Smolt hatchery releases from escapement"),
    natural_origin = natural_origin
  )
  structure(out, class = "RelSmolt")
}

#' @importFrom stats predict
#' @export
predict.RelSmolt <- function(object, newdata, ...) {

  if (missing(newdata)) return(object$fitted.values)

  vars <- names(newdata)

  vars_Rel <- names(object$model)[-1]
  Nage <- vars_Rel[grepl("Nage", vars_Rel)]
  Nage_NOS <- Nage[object$natural_origin]
  Nage_HOS <- Nage[!object$natural_origin]

  do_hatchery <- length(Nage_HOS) > 0

  val <- sapply(1:nrow(newdata), function(i) {
    Esc_NOS <- sapply(Nage_NOS, function(j) as.numeric(newdata[i, grepl(j, vars)])) # matrix n_age x ng
    n_age <- nrow(Esc_NOS)
    a <- seq(3, n_age, 2)
    if (sum(Esc_NOS[!a, ])) stop("Internal salmonMSE error: there is natural escapement in even age classes")

    if (do_hatchery) {
      Esc_HOS <- sapply(Nage_HOS, function(j) as.numeric(newdata[i, grepl(j, vars)]))
      if (sum(Esc_HOS[!a, ])) stop("Internal salmonMSE error: there is hatchery escapement in even age classes")
    } else {
      Esc_HOS <- matrix(0, n_age, 1)
    }
    x <- newdata[i, "x"]
    y <- newdata[i, "y"]
    object$func(Nage_NOS = Esc_NOS[a, , drop = FALSE], Nage_HOS = Esc_HOS[a, , drop = FALSE], x = x, y = y)
  })

  return(val)

}

#' @importFrom stats simulate
#' @export
simulate.RelSmolt <- function(object, nsim = 1, seed = 1, ...) {

  set.seed(seed)
  stdev <- sqrt(log(1 + object$CV^2))
  do_sim <- object$fitted.values * rlnorm(nsim * length(object$fitted.values), -0.5 * stdev^2, stdev)
  val <- matrix(do_sim, length(object$fitted.values), nsim) %>% as.data.frame() %>%
    structure(names = paste0("sim_", 1:nsim), seed = seed)
  val
}

#' Update natural mortality of juveniles
#'
#' @description Internal function that updates juvenile natural mortality in the marine environement due
#' to fitness
#' - `SAR_fitness()` calculates the new natural mortality value
#' - `makeRel_SAR` generates a list for openMSE to use in the simulations
#' @param x Integer, simulation number from openMSE
#' @param y Integer, simulation year (including historical years)
#' @param envir Character, whether to obtain the fitness value for the natural or hatchery environment.
#' @param rel_loss Numeric, the loss exponent for the juveniles
#' @param s Integer, the salmonMSE population index. Used to search for the fitness value
#' @param nyears Integer, the number of historical years in the operating model
#' @param Mbase Array `[nsim, n_age, proyears]` the base natural mortality value in the openMSE operating model.
#' @returns
#' - `smolt_func()` returns a numeric for the ratio of the realized smolt production vs. the hypothetical value if there were no
#' hatchery, en route mortality, or habitat improvement
#' - `makeRel_smolt` returns a list that is passed to openMSE as a inter-population relationship
#' @keywords internal
SAR_fitness <- function(x = -1, y = 1,
                        envir = c("natural", "hatchery"),
                        rel_loss = 1, s = 1,
                        nyears, Mbase) {

  envir <- match.arg(envir)

  even_time_step <- !y %% 2
  if (even_time_step) {
    n_age <- dim(Mbase)[2]
    fitness_a <- rep(1, n_age)
    for (a in seq(2, n_age, 2)) {
      if (nrow(salmonMSE_env$Ford) && x > 0) {
        fitness_y <- dplyr::filter(
          salmonMSE_env$Ford,
          .data$x == .env$x, .data$s == .env$s, .data$t == y - a, .data$type == envir
        ) %>%
          pull(.data$fitness)

        if (length(fitness_y)) fitness_a[a] <- fitness_y
      }
    }

    fitness_loss <- fitness_a^rel_loss
    M <- Mbase[x, , y - nyears]
    SAR_base <- SAR_loss <- exp(-M)
    SAR_loss[M > .Machine$double.eps] <- SAR_base[M > .Machine$double.eps] * fitness_loss[M > .Machine$double.eps]
    return(-log(SAR_loss))
  } else {
    return(Mbase[x, , y - nyears])
  }
}

#' @rdname SAR_fitness
#' @param p_smolt Integer, the population index for the juvenile population in the openMSE model
#' @keywords internal
makeRel_SAR <- function(p_smolt = 1, s = 1, envir = c("natural", "hatchery"),
                        rel_loss, nyears, Mbase) {

  envir <- match.arg(envir)

  .SAR_fitness <- SAR_fitness

  formals(.SAR_fitness)$s <- s
  formals(.SAR_fitness)$rel_loss <- rel_loss
  formals(.SAR_fitness)$envir <- envir
  formals(.SAR_fitness)$nyears <- nyears
  formals(.SAR_fitness)$Mbase <- Mbase

  model <- data.frame(M = 1, N = seq(0, 10), x = -1, y = 1)

  response <- paste0("M_", p_smolt)
  input <- paste0("N_", p_smolt)
  terms <- c(response, input, "x", "y")

  out <- list(
    func = .SAR_fitness,
    model = structure(model, names = terms),
    fitted.values = model$M,
    CV = 0,
    terms = terms,
    type = "Natural mortality",
    Rel = "Marine survival reduced by fitness",
    age = seq(0, dim(Mbase)[2] - 1)
  )
  structure(out, class = "SARfitness")
}

#' @export
predict.SARfitness <- function(object, newdata, ...) {
  if (missing(newdata)) return(object$fitted.values)
  val <- sapply(1:nrow(newdata), function(i) object$func(x = newdata[i, "x"], y = newdata[i, "y"]))
  return(val)
}

#' @export
simulate.SARfitness <- simulate.RelSmolt
