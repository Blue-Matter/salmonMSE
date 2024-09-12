

.smolt_func <- function(Nage, x = -1, y, output = c("natural", "hatchery"),
                        ptarget_NOB, pmax_NOB, egg_local, fec_brood,
                        s_yearling, s_subyearling, p_yearling,
                        phatchery, premove_HOS, s_prespawn, # Broodtake & hatchery production
                        p_female, fec, gamma, SRRpars, # Spawning (natural production)
                        fitness_args, p_naturalsmolt) {

  output <- match.arg(output)
  Nage[is.na(Nage)] <- 0

  # Hatchery
  broodtake <- calc_broodtake(Nage, ptarget_NOB, pmax_NOB, phatchery, egg_local, p_female, fec_brood, s_prespawn)

  egg_NOB <- sum(broodtake$NOB * fec_brood * s_prespawn * p_female)
  egg_HOB <- sum(broodtake$HOB * fec_brood * s_prespawn * p_female)

  hatchery_production <- calc_yearling(egg_NOB + egg_HOB, s_yearling, s_subyearling, p_yearling)
  yearling <- hatchery_production[1]
  subyearling <- hatchery_production[2]

  # Spawners
  spawners <- calc_spawners(broodtake, Nage, phatchery, premove_HOS)

  NOS <- spawners$NOS
  HOS <- spawners$HOS
  if (sum(HOS)) {
    HOS_effective <- spawners$HOS * gamma
  } else {
    HOS_effective <- rep(0, length(HOS))
  }

  # Spawners weighted by fecundity
  pNOB <- sum(fec_brood * broodtake$NOB)/sum(fec_brood * broodtake$NOB, fec_brood * broodtake$HOB)
  pHOSeff <- sum(fec * HOS_effective)/sum(fec * NOS, fec *HOS_effective)
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
    zbar_prev <- filter(salmonMSE_env$Ford, x == .env$x, p_smolt == .env$p_naturalsmolt)

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
        NOS, HOS_effective, broodtake$NOB, broodtake$HOB, fec, fec_brood, zbar_brood,
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
          zbar[i], fitness_args$theta[i], fitness_args$omega2, fitness_args$fitness_variance, fitness_args$fitness_floor
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

  # Predicted smolts from historical SRR parameters and openMSE setup (if there were no hatchery production or habitat improvement)
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

makeRel_smolt <- function(p_smolt = 1, p_naturalsmolt = 1, p_natural, p_hatchery,
                          output = c("natural", "hatchery"),
                          ptarget_NOB, pmax_NOB, egg_local, fec_brood,
                          s_yearling, s_subyearling, p_yearling, phatchery, premove_HOS, s_prespawn, # Broodtake & hatchery production
                          p_female, fec, gamma, SRRpars,  # Spawning (natural production)
                          fitness_args) {

  output <- match.arg(output)

  func <- .smolt_func

  formals(func)$ptarget_NOB <- ptarget_NOB
  formals(func)$pmax_NOB <- pmax_NOB
  formals(func)$egg_local <- egg_local
  formals(func)$fec_brood <- fec_brood

  formals(func)$s_yearling <- s_yearling
  formals(func)$s_subyearling <- s_subyearling
  formals(func)$p_yearling <- p_yearling

  formals(func)$phatchery <- phatchery
  formals(func)$s_prespawn <- s_prespawn
  formals(func)$p_female <- p_female
  formals(func)$output <- output
  formals(func)$gamma <- gamma

  formals(func)$premove_HOS <- premove_HOS
  formals(func)$fec <- fec

  formals(func)$SRRpars <- SRRpars

  formals(func)$fitness_args <- fitness_args
  formals(func)$p_naturalsmolt <- p_naturalsmolt

  maxage <- length(fec)
  N_natural <- rep(1, maxage) %>% structure(names = paste0("Nage_", p_natural, 1:maxage))

  if (!is.na(p_hatchery)) {
    N_hatchery <- rep(0, maxage) %>% structure(names = paste0("Nage_", p_hatchery, 1:maxage))
    Nage <- cbind(N_natural, N_hatchery)
    Perr_y <- func(Nage, x = -1, y = 1)
    model <- c(Perr_y = Perr_y, sum(N_natural), sum(N_hatchery), x = -1, y = 1)
    input <- paste0("Nage_", c(p_natural, p_hatchery))
  } else {
    Nage <- matrix(N_natural, ncol = 1)
    Perr_y <- func(Nage, x = -1, y = 1)
    model <- c(Perr_y = Perr_y, sum(N_natural), x = -1, y = 1)
    input <- paste0("Nage_", p_natural)
  }

  response <- paste0("Perr_y_", p_smolt)
  terms <- c(response, input, "x", "y")

  out <- list(
    func = func,
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

# Natural mortality in the marine environment of natural origin smolts, increased due to fitness loss
# Only applied at the start of even time step
.SAR_fitness <- function(x = -1, y = 1,
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


makeRel_SAR <- function(p_smolt = 1, p_naturalsmolt = p_smolt, envir = c("natural", "hatchery"), rel_loss, nyears, Mbase) {

  envir <- match.arg(envir)

  func <- .SAR_fitness

  #formals(func)$p_smolt <- p_smolt
  formals(func)$p_naturalsmolt <- p_naturalsmolt
  formals(func)$rel_loss <- rel_loss
  formals(func)$envir <- envir
  formals(func)$nyears <- nyears
  formals(func)$Mbase <- Mbase

  model <- data.frame(M = 1, N = seq(0, 10), x = -1, y = 1)

  response <- paste0("M_", p_smolt)
  input <- paste0("N_", p_smolt)
  terms <- c(response, input, "x", "y")

  out <- list(
    func = func,
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
