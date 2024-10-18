
#' Predict smolt production
#'
#' @description Internal function that predicts the natural origin or hatchery origin smolt production from the
#' escapement, and saves state variables to [salmonMSE::salmonMSE_env].
#' - `smolt_func()` is the population dynamics function
#' - `makeRel_smolt()` generates a list for openMSE to use in the simulations
#'
#' @param Nage Matrix `[n_age, 2]` of natural and hatchery escapement from openMSE
#' @param x Integer, simulation number from openMSE
#' @param y Integer, simulation year (including historical years)
#' @param output Character, whether to predic the natural origin or hatchery origin smolt production
#' @param s_enroute Numeric, en route survival of the escapement to the spawning grounds
#' @param p_female Numeric, proportion female for calculating the egg production.
#' @param fec Vector `[maxage]`. The fecundity at age schedule of spawners
#' @param SRRpars Data frame containing stock recruit parameters for natural smolt production.
#' Column names include: SRrel, kappa, capacity_smolt, Smax, phi, kappa_improve, capacity_smolt_improve
#' @param hatchery_args Named list containing various arguments controlling broodtake and hatchery production. See details below.
#' @param fitness_args Named list containing various arguments controlling population fitness from hatchery production.
#' Names include: fitness_type, omega2, theta, fitness_variance, heritability, zbar_start, fitness_floor, rel_loss
#' @param p_naturalsmolt Integer, the population index for the natural smolts in the openMSE model.
#' Used to report variables to [salmonMSE::salmonMSE_env].
#' @section hatchery_args:
#' Hatchery control parameters are included in a named list with the following arguments:
#'
#' - `egg_local` Numeric, target egg production for hatchery. Set to zero for no hatchery production.
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
smolt_func <- function(Nage, x = -1, y, output = c("natural", "hatchery"),
                       s_enroute, p_female, fec, SRRpars, # Spawning (natural production)
                       hatchery_args, fitness_args, p_naturalsmolt) {

  output <- match.arg(output)
  Nage[is.na(Nage)] <- 0

  # Hatchery
  Nage_enroute <- Nage * s_enroute * hatchery_args$pmax_esc
  broodtake <- calc_broodtake(
    Nage_enroute,
    hatchery_args$ptarget_NOB,
    hatchery_args$pmax_NOB,
    hatchery_args$phatchery,
    hatchery_args$egg_local,
    hatchery_args$p_female,
    hatchery_args$fec_brood,
    hatchery_args$s_prespawn
  )

  egg_NOB <- sum(broodtake$NOB * hatchery_args$fec_brood * hatchery_args$s_prespawn * p_female)
  egg_HOB <- sum(broodtake$HOB * hatchery_args$fec_brood * hatchery_args$s_prespawn * p_female)

  hatchery_production <- calc_yearling(
    egg_NOB + egg_HOB,
    hatchery_args$s_yearling,
    hatchery_args$s_subyearling,
    hatchery_args$p_yearling
  )
  yearling <- hatchery_production[1]
  subyearling <- hatchery_production[2]

  # Spawners
  spawners <- calc_spawners(broodtake, Nage_enroute, hatchery_args$phatchery, hatchery_args$premove_HOS)

  NOS <- spawners$NOS
  HOS <- spawners$HOS
  if (sum(HOS)) {
    HOS_effective <- spawners$HOS * hatchery_args$gamma
  } else {
    HOS_effective <- rep(0, length(HOS))
  }

  # Spawners weighted by fecundity
  pNOB <- sum(hatchery_args$fec_brood * broodtake$NOB)/
    sum(hatchery_args$fec_brood * broodtake$NOB, hatchery_args$fec_brood * broodtake$HOB)
  pHOSeff <- sum(fec * HOS_effective)/sum(fec * NOS, fec * HOS_effective)
  pHOScensus <- sum(fec * HOS)/sum(fec * NOS, fec * HOS)

  # Natural egg production in the absence of fitness effects
  Egg_NOS <- sum(NOS * p_female * fec)
  Egg_HOS <- sum(HOS_effective * p_female * fec)

  total_egg <- Egg_NOS + Egg_HOS
  if (output == "natural" && !total_egg) return(0) # Perr_y = 0
  if (output == "hatchery" && !sum(hatchery_production)) return(0) # Perr_y = 0

  # Fitness
  if (x > 0 && any(fitness_args$fitness_type == "Ford")) {
    #browser(expr = x == 2)
    # Get zbar from salmonMSE_env
    zbar_prev <- filter(salmonMSE_env$Ford, x == .env$x, .data$p_smolt == .env$p_naturalsmolt)

    if (nrow(zbar_prev)) {
      maxage <- length(NOS)
      zbar_brood <- matrix(0, maxage, 2) # Column 1 = natural environment, 2 = hatchery environment

      # Trait value by brood year
      for (a in 1:maxage) {
        if (NOS[a] || broodtake$NOB[a]) {
          zbar1 <- dplyr::filter(zbar_prev, .data$t == .env$y - 2 * .env$a, .data$type == "natural") %>%
            pull(.data$zbar)
          if (!length(zbar1)) stop("Cannot find zbar") # zbar1 <- fitness_args$zbar_start[1]
          zbar_brood[a, 1] <- zbar1
        }
        if (HOS_effective[a] || broodtake$HOB[a]) {
          zbar2 <- dplyr::filter(zbar_prev, .data$t == .env$y - 2 * .env$a, .data$type == "hatchery") %>%
            pull(.data$zbar)
          if (!length(zbar2)) stop("Cannot find zbar") # zbar2 <- fitness_args$zbar_start[2]
          zbar_brood[a, 2] <- zbar2
        }
      }

      zbar <- calc_zbar(
        NOS, HOS_effective, broodtake$NOB, broodtake$HOB, fec, hatchery_args$fec_brood, zbar_brood,
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
    fitness <- c(1, 1)
    fitness_loss <- matrix(1, 2, 3)
  }

  # Natural egg production after fitness loss
  Egg_NOS_out <- Egg_NOS * fitness_loss[1, 1]
  Egg_HOS_out <- Egg_HOS * fitness_loss[1, 1]

  # Hatchery production after fitness loss
  yearling_out <- yearling * fitness_loss[2, 1] * fitness_loss[2, 2]
  subyearling_out <- subyearling * fitness_loss[2, 1]
  total_egg_out <- Egg_NOS_out + Egg_HOS_out + subyearling_out

  # Smolt production (natural, but in competition with hatchery subyearlings) after fitness loss
  xx <- max(x, 1)
  SRrel <- match.arg(SRRpars[xx, "SRrel"], choices = c("BH", "Ricker"))

  smolt_subyearling <- calc_smolt(
    subyearling_out, total_egg_out,
    SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRRpars[xx, "kappa_improve"], SRRpars[xx, "capacity_smolt_improve"], fitness_loss[2, 2],
    SRrel
  )
  smolt_rel <- yearling_out + smolt_subyearling

  # Return total hatchery smolt releases after density-dependent survival of subyearlings
  if (output == "hatchery") {
    return(smolt_rel)
  }

  smolt_NOS_proj <- calc_smolt(
    Egg_NOS_out, total_egg_out,
    SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRRpars[xx, "kappa_improve"], SRRpars[xx, "capacity_smolt_improve"], fitness_loss[1, 2],
    SRrel
  )

  smolt_HOS_proj <- calc_smolt(
    Egg_HOS_out, total_egg_out,
    SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRRpars[xx, "kappa_improve"], SRRpars[xx, "capacity_smolt_improve"], fitness_loss[1, 2],
    SRrel
  )

  # Predicted smolts from historical SRR parameters and openMSE setup
  # if there were no hatchery production, habitat improvement, or enroute mortality
  Egg_openMSE <- sum(Nage[, 1] * p_female * fec)
  smolt_NOS_SRR <- calc_smolt(
    Egg_openMSE, Egg_openMSE, SRRpars[xx, "kappa"], SRRpars[xx, "capacity_smolt"], SRRpars[xx, "Smax"], SRRpars[xx, "phi"],
    SRrel = SRrel
  )

  # MICE parameter to be updated in operating model
  Perr_y <- (smolt_NOS_proj + smolt_HOS_proj)/smolt_NOS_SRR

  if (x > 0) {
    # Save state variables
    df_N <- data.frame(
      x = x,
      p_smolt = p_naturalsmolt,
      t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
      a = 1:nrow(Nage),
      Esc_NOS = Nage[, 1],
      Esc_HOS = Nage[, 2],
      NOB = broodtake$NOB,
      HOB = broodtake$HOB,
      NOS = NOS,
      HOS = HOS,
      HOS_effective = HOS_effective
    )
    salmonMSE_env$N <- rbind(salmonMSE_env$N, df_N)

    df_state <- data.frame(
      x = x,
      p_smolt = p_naturalsmolt,
      t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
      Egg_NOS = Egg_NOS_out,
      Egg_HOS = Egg_HOS_out,
      smolt_NOS = smolt_NOS_proj,
      smolt_HOS = smolt_HOS_proj,
      fitness_natural = fitness[1],
      fitness_hatchery = fitness[2],
      pNOB = as.numeric(pNOB),
      pHOSeff = as.numeric(pHOSeff),
      pHOScensus = as.numeric(pHOScensus),
      Perr_y = Perr_y,
      #alpha = alpha_proj,
      #beta = beta_proj,
      yearling = yearling_out,
      subyearling = subyearling_out,
      smolt_rel = smolt_rel
    )
    salmonMSE_env$state <- rbind(salmonMSE_env$state, df_state)

    # Save zbar
    if (any(fitness_args$fitness_type == "Ford")) {
      df_Ford <- data.frame(
        x = x,
        p_smolt = p_naturalsmolt,
        t = y, # Even time steps (remember MICE predicts Perr_y for next time step)
        type = c("natural", "hatchery"),
        zbar = zbar
      )
      salmonMSE_env$Ford <- rbind(salmonMSE_env$Ford, df_Ford)
    }
  }

  return(Perr_y)
}

#' @rdname smolt_func
#' @param p_smolt Integer, the population index for the smolt production in the openMSE model, corresponding to `output`
#' @param p_natural Integer, the population index for the natural origin escapement in the openMSE model
#' @param p_hatchery Integer, the population index for the hatchery origin escapement in the openMSE model
makeRel_smolt <- function(p_smolt = 1, p_naturalsmolt = 1, p_natural, p_hatchery,
                          output = c("natural", "hatchery"), s_enroute,
                          p_female, fec, SRRpars,  # Spawning (natural production)
                          hatchery_args, fitness_args) {

  output <- match.arg(output)

  .smolt_func <- smolt_func

  formals(.smolt_func)$output <- output
  formals(.smolt_func)$s_enroute <- s_enroute
  formals(.smolt_func)$p_female <- p_female
  formals(.smolt_func)$fec <- fec

  formals(.smolt_func)$SRRpars <- SRRpars

  formals(.smolt_func)$hatchery_args <- hatchery_args
  formals(.smolt_func)$fitness_args <- fitness_args
  formals(.smolt_func)$p_naturalsmolt <- p_naturalsmolt

  maxage <- length(fec)
  N_natural <- rep(1, maxage) %>% structure(names = paste0("Nage_", p_natural, 1:maxage))

  if (!is.na(p_hatchery)) {
    N_hatchery <- rep(0, maxage) %>% structure(names = paste0("Nage_", p_hatchery, 1:maxage))
    Nage <- cbind(N_natural, N_hatchery)
    Perr_y <- .smolt_func(Nage, x = -1, y = 1)
    model <- c(Perr_y = Perr_y, sum(N_natural), sum(N_hatchery), x = -1, y = 1)
    input <- paste0("Nage_", c(p_natural, p_hatchery))
  } else {
    Nage <- matrix(N_natural, ncol = 1)
    Perr_y <- .smolt_func(Nage, x = -1, y = 1)
    model <- c(Perr_y = Perr_y, sum(N_natural), x = -1, y = 1)
    input <- paste0("Nage_", p_natural)
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
                 "hatchery" = "Smolt hatchery releases from escapement")
  )
  structure(out, class = "RelSmolt")
}

#' @importFrom stats predict
#' @export
predict.RelSmolt <- function(object, newdata, ...) {

  if (missing(newdata)) return(object$fitted.values)

  vars <- names(newdata)

  do_hatchery <- length(names(object$model)) == 5 # Perr_y, Nage, Nage, x, y
  vars_Rel <- names(object$model)[-1]

  val <- sapply(1:nrow(newdata), function(i) {
    Esc_NOS <- newdata[i, grepl(vars_Rel[1], vars)] %>% as.numeric()
    n_age <- length(Esc_NOS)
    a <- seq(3, n_age, 2)
    if (sum(Esc_NOS[!a])) stop("Internal salmonMSE error: there is natural escapement in even age classes")

    if (do_hatchery) {
      Esc_HOS <- newdata[i, grepl(vars_Rel[2], vars)] %>% as.numeric()
      if (sum(Esc_HOS[!a])) stop("Internal salmonMSE error: there is hatchery escapement in even age classes")
      Nage <- cbind(Esc_NOS[a], Esc_HOS[a])
    } else {
      Nage <- matrix(Esc_NOS[a], ncol = 1)
    }
    x <- newdata[i, "x"]
    y <- newdata[i, "y"]
    object$func(Nage = Nage, x = x, y = y)
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
#' @param p_naturalsmolt Integer, the population index for the natural smolts in the openMSE model. Used to search for the fitness value
#' @param nyears Integer, the number of historical years in the operating model
#' @param Mbase Array `[nsim, n_age, proyears]` the base natural mortality value in the openMSE operating model.
#' @returns
#' - `smolt_func()` returns a numeric for the ratio of the realized smolt production vs. the hypothetical value if there were no
#' hatchery, en route mortality, or habitat improvement
#' - `makeRel_smolt` returns a list that is passed to openMSE as a inter-population relationship
#' @keywords internal
SAR_fitness <- function(x = -1, y = 1,
                        envir = c("natural", "hatchery"),
                        rel_loss = 1, p_naturalsmolt = 1,
                        nyears, Mbase) {

  envir <- match.arg(envir)

  even_time_step <- !y %% 2
  if (even_time_step) {
    n_age <- dim(Mbase)[2]
    fitness_a <- rep(1, n_age)
    for (a in seq(2, n_age, 2)) {
      if (nrow(salmonMSE_env$state) && x > 0) {
        fitness_y <- dplyr::filter(
          salmonMSE_env$state,
          .data$x == .env$x, .data$p_smolt == .env$p_naturalsmolt, .data$t == y - a
        ) %>%
          pull(.data[[paste0("fitness_", envir)]])

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
makeRel_SAR <- function(p_smolt = 1, p_naturalsmolt = p_smolt, envir = c("natural", "hatchery"),
                        rel_loss, nyears, Mbase) {

  envir <- match.arg(envir)

  .SAR_fitness <- SAR_fitness

  formals(.SAR_fitness)$p_naturalsmolt <- p_naturalsmolt
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
