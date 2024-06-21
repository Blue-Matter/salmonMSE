
## Internal functions to predict either:
# (b) smolt hatchery production
# (c) smolt natural production
calc_broodtake <- function(Nage, ptarget_NOB, pmax_NOB, phatchery, egg_local, p_female,
                           fec_brood, gamma, s_prespawn) {

  NOR_escapement <- Nage[, 1]
  if (ncol(Nage) > 1) {
    HOR_escapement <- Nage[, 2]

    opt <- try(
      uniroot(.broodtake_func, c(1e-8, pmax_NOB),
              NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
              fec = fec_brood, gamma = gamma, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB),
      silent = TRUE
    )
    if (is.character(opt)) {
      ptake_NOB <- pmax_NOB
    } else {
      ptake_NOB <- opt$root
    }

    output <- .broodtake_func(
      ptake_NOB, NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
      fec = fec_brood, gamma = gamma, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, opt = FALSE
    )

    NOB <- output$NOB
    HOB <- output$HOB
  } else {

    NOBopt <- try(
      uniroot(.egg_func, c(1e-8, pmax_NOB), N = NOR_escapement, gamma = 1, fec = fec_brood, p_female = p_female,
              s_prespawn = s_prespawn, val = egg_local),
      silent = TRUE
    )
    if (is.character(NOBopt)) {
      ptake_NOB <- pmax_NOB
    } else {
      ptake_NOB <- NOBopt$root
    }
    NOB <- ptake_NOB * NOR_escapement
    HOB <- rep(0, length(NOB))
  }

  list(NOB = NOB, HOB = HOB)
}

calc_spawners <- function(broodtake, escapement, phatchery, premove_HOS) {
  spawners <- list()
  spawners$NOS <- escapement[, 1] - broodtake$NOB
  if (length(escapement) > 1) spawners$HOS <- escapement[, 2] * (1 - phatchery) * (1 - premove_HOS)
  return(spawners)
}


.egg_func <- function(ptake, N, gamma = 1, fec, p_female, s_prespawn, val = 0) sum(ptake * gamma * N * s_prespawn * fec * p_female) - val

.broodtake_func <- function(ptake_NOB, NOR_escapement, HOR_escapement, phatchery, p_female, fec, gamma,
                            egg_target, s_prespawn, ptarget_NOB, opt = TRUE) {

  NOB <- ptake_NOB * NOR_escapement
  egg_NOB <- sum(NOB * s_prespawn * fec * p_female)
  egg_HOB <- egg_target - egg_NOB

  if (egg_HOB < 0) {
    ptake_HOB <- 0
  } else {
    # solve for the required HOB
    HOBopt <- try(
      uniroot(.egg_func, c(1e-8, 1), N = HOR_escapement * phatchery, gamma = gamma, fec = fec, p_female = p_female,
              s_prespawn = s_prespawn, val = egg_HOB),
      silent = TRUE
    )
    if (is.character(HOBopt)) {
      ptake_HOB <- 1
    } else {
      ptake_HOB <- HOBopt$root
    }
  }

  HOB <- ptake_HOB * HOR_escapement * phatchery
  egg_HOB <- .egg_func(ptake_HOB, N = HOR_escapement * phatchery, gamma = gamma, fec = fec, p_female = p_female, s_prespawn = s_prespawn)

  if (opt) {
    obj <- log(sum(NOB)/sum(NOB + HOB)) - log(ptarget_NOB)   # Objective function, get close to zero, if not stay positive
    return(obj)
  } else {
    output <- list(egg_NOB = egg_NOB, egg_HOB = egg_HOB, ptake_HOB = ptake_HOB, NOB = NOB, HOB = HOB)
    return(output)
  }
}


.smolt_func <- function(Nage, x = -1, output = c("natural", "hatchery"),
                        ptarget_NOB, pmax_NOB, egg_local, fec_brood, s_egg, phatchery, premove_HOS, s_prespawn, # Broodtake & hatchery production
                        p_female, fec, gamma, # Spawning (natural production)
                        SRRpars_hist, SRRpars_proj, SRrel = c("BH", "Ricker"),
                        fitness_type = c("Ford", "none"), # Spawning (natural production)
                        fitness_args, p_smolt) {

  output <- match.arg(output)

  Nage[is.na(Nage)] <- 0
  broodtake <- calc_broodtake(Nage, ptarget_NOB, pmax_NOB, phatchery, egg_local, p_female, fec_brood, gamma, s_prespawn)

  if (output == "hatchery") {
    fry <- sum(broodtake$NOB + broodtake$HOB * fec_brood * s_prespawn * p_female)
    smolt <- fry * s_egg
    return(smolt)

  } else if (output == "natural") {
    SRrel <- match.arg(SRrel)
    fitness_type <- match.arg(fitness_type)
    spawners <- calc_spawners(broodtake, Nage, phatchery, premove_HOS)

    NOS <- spawners$NOS
    HOS <- spawners$HOS
    HOS_effective <- spawners$HOS * gamma

    # Fry production in the absence of fitness effects
    fry_NOS <- sum(NOS * p_female * fec)
    fry_HOS <- sum(HOS_effective * p_female * fec)

    total_fry <- fry_NOS + fry_HOS

    if (!total_fry) return(0) # Perr_y = 0

    alpha_hist <- SRRpars_hist[max(x, 1), "a"]
    beta_hist <- SRRpars_hist[max(x, 1), "b"]

    alpha_proj <- SRRpars_proj[max(x, 1), "a"]
    beta_proj <- SRRpars_proj[max(x, 1), "b"]

    # Predicted smolts from historical SRR parameters and openMSE setup (if there were no hatchery production)
    fry_openMSE <- sum(Nage[, 1] * p_female * fec)
    if (SRrel == "BH") {
      smolt_NOS_SRR <- alpha_hist * fry_openMSE/(1 + beta_hist * fry_openMSE)
    } else {
      smolt_NOS_SRR <- alpha_hist * fry_openMSE * exp(-beta_hist * fry_openMSE)
    }

    # Predicted fry and smolts from projected SRR parameters and fitness
    if (egg_local > 0 && fitness_type == "Ford" && x > 0) {
      pNOB <- sum(broodtake$NOB)/sum(broodtake$NOB, broodtake$HOB)
      pHOSeff <- sum(HOS_effective)/sum(NOS, HOS_effective)
      pHOScensus <- sum(HOS)/sum(NOS, HOS)

      # Get pbar from salmonMSE_env
      if (nrow(salmonMSE_env$Ford)) {
        Ford_x <- filter(salmonMSE_env$Ford, x == .env$x, p_smolt == .env$p_smolt)
        if (nrow(Ford_x)) {
          pbar_prev <- filter(Ford_x, t == max(.data$t)) %>% pull(.data$pbar)
          tprev <- max(Ford_x$t)
        } else {
          pbar_prev <- numeric(0)
          tprev <- 0
        }
      } else {
        tprev <- 0
      }
      if (nrow(salmonMSE_env$Ford) && length(pbar_prev)) {
        fitness_calcs <- calc_Ford_fitness(
          pNOB, pHOSeff, pbar_prev, fitness_args$theta, fitness_args$omega2, fitness_args$fitness_variance,
          fitness_args$A, fitness_args$fitness_floor, fitness_args$heritability
        )
        pbar <- fitness_calcs$pbar
        fitness <- fitness_calcs$fitness
      } else {
        pbar <- fitness_args$pbar_start
        fitness <- local({
          x <- exp(-0.5 * (pbar[1] - fitness_args$theta[1])^2/(fitness_args$omega2+fitness_args$fitness_variance))
          max(fitness_args$fitness_floor, x)
        })
      }
      fitness_loss <- fitness^fitness_args$rel_loss

      # Fry_NOS with fitness
      fry_NOS_out <- fry_NOS * fitness_loss[1]
      fry_HOS_out <- fry_HOS * fitness_loss[1]

      # Update smolt SRR with fitness
      alpha_fitness <- alpha_proj * fitness_loss[2]
      beta_fitness <- beta_proj / fitness_loss[2]

      # Save pbar for next generation
      df_Ford <- data.frame(
        x = x,
        p_smolt = p_smolt,
        t = tprev + 1,
        type = c("natural", "hatchery"),
        pbar = pbar
      )
      salmonMSE_env$Ford <- rbind(salmonMSE_env$Ford, df_Ford)

    } else { #if (fitness_type == "none") {
      alpha_fitness <- alpha_proj
      beta_fitness <- beta_proj

      fitness <- 1
      pNOB <- 1
      pHOSeff <- 0

      fry_NOS_out <- fry_NOS
      fry_HOS_out <- fry_HOS
    }

    total_fry_out <- fry_NOS_out + fry_HOS_out
    if (SRrel == "BH") {
      smolt_NOS_proj <- alpha_proj * fry_NOS_out/(1 + beta_proj * total_fry_out)
      smolt_HOS_proj <- alpha_proj * fry_HOS_out/(1 + beta_proj * total_fry_out)
    } else {
      smolt_NOS_proj <- alpha_proj * fry_NOS_out * exp(-beta_proj * total_fry_out)
      smolt_HOS_proj <- alpha_proj * fry_HOS_out * exp(-beta_proj * total_fry_out)
    }
    total_smolt <- smolt_NOS_proj + smolt_HOS_proj
    Perr_y <- total_smolt/smolt_NOS_SRR

    if (x > 0) {
      # Save state variables
      if (nrow(salmonMSE_env$N)) {
        tprevN <- local({
          dat <- filter(salmonMSE_env$N, x == .env$x, p_smolt == .env$p_smolt)
          ifelse(nrow(dat), max(dat$t), 0)
        })
      } else {
        tprevN <- 0
      }

      df_N <- data.frame(
        x = x,
        p_smolt = p_smolt,
        t = tprevN + 1,
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

      if (nrow(salmonMSE_env$state)) {
        tprevstate <- local({
          dat <- filter(salmonMSE_env$state, x == .env$x, p_smolt == .env$p_smolt)
          ifelse(nrow(dat), max(dat$t), 0)
        })
      } else {
        tprevstate <- 0
      }
      if (tprevstate != tprevN) stop("Problem with saving state variables")
      df_state <- data.frame(
        x = x,
        p_smolt = p_smolt,
        t = tprevstate + 1,
        fry_NOS = fry_NOS_out,
        fry_HOS = fry_HOS_out,
        smolt_NOS = smolt_NOS_proj,
        smolt_HOS = smolt_HOS_proj,
        fitness = fitness,
        pNOB = as.numeric(pNOB),
        pHOSeff = as.numeric(pHOSeff),
        pHOScensus = as.numeric(pHOScensus),
        Perr_y = Perr_y
      )
      salmonMSE_env$state <- rbind(salmonMSE_env$state, df_state)
    }

    return(Perr_y)
  }
}

makeRel_smolt <- function(p_smolt = 1, p_natural, p_hatchery,
                          output = c("natural", "hatchery"),
                          ptarget_NOB, pmax_NOB, egg_local, fec_brood, s_egg, phatchery, premove_HOS, s_prespawn, # Broodtake & hatchery production
                          p_female, fec, gamma, # Spawning (natural production)
                          SRRpars_hist, SRRpars_proj, SRrel = c("BH", "Ricker"),  # Spawning (natural production)
                          fitness_type = c("Ford", "none"), fitness_args) {

  output <- match.arg(output)

  func <- .smolt_func

  formals(func)$ptarget_NOB <- ptarget_NOB
  formals(func)$pmax_NOB <- pmax_NOB
  formals(func)$egg_local <- egg_local
  formals(func)$fec_brood <- fec_brood
  formals(func)$s_egg <- s_egg
  formals(func)$phatchery <- phatchery
  formals(func)$premove_HOS <- premove_HOS
  formals(func)$s_prespawn <- s_prespawn
  formals(func)$p_female <- p_female
  formals(func)$output <- output

  maxage <- length(fec)

  if (output == "natural") {

    SRrel <- match.arg(SRrel)

    formals(func)$p_female <- p_female
    formals(func)$fec <- fec
    formals(func)$gamma <- gamma

    formals(func)$SRRpars_hist <- SRRpars_hist
    formals(func)$SRRpars_proj <- SRRpars_proj
    formals(func)$SRrel <- SRrel
    formals(func)$fitness_type <- fitness_type

    if (fitness_type != "none") {
      formals(func)$fitness_args <- fitness_args
      formals(func)$p_smolt <- p_smolt
    }
  }

  N_natural <- rep(1, maxage) %>% structure(names = paste0("Nage_", p_natural, 1:maxage))
  N_hatchery <- rep(0, maxage) %>% structure(names = paste0("Nage_", p_hatchery, 1:maxage))

  Perr_y <- func(cbind(N_natural, N_hatchery))

  model <- c(Perr_y = Perr_y, N_natural, N_hatchery, x = -1) %>% t() %>% as.data.frame()

  response <- paste0("Perr_y_", p_smolt)
  input <- paste0("Nage_", c(p_natural, p_hatchery))
  terms <- c(response, input, "x")

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
  vars_Rel <- names(object$model)[-1]

  val <- sapply(1:nrow(newdata), function(i) {
    Esc_NOS <- newdata[i, grepl(vars_Rel[1], vars)] %>% as.numeric()
    Esc_HOS <- newdata[i, grepl(vars_Rel[2], vars)] %>% as.numeric()
    Nage <- cbind(Esc_NOS, Esc_HOS)
    x <- newdata[i, "x"]
    object$func(Nage = Nage, x = x)
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

# Marine survival of natural origin smolts, as reduced by fitness
# N is a dummy variable
.SAR_fitness <- function(N, x = -1,
                         fitness_type = c("Ford", "none"), # Spawning (natural production)
                         rel_loss = 1, p_naturalsmolt = 1) {

  fitness_type <- match.arg(fitness_type)
  stopifnot(fitness_type == "Ford")

  if (nrow(salmonMSE_env$state) && x > 0) {
    fitness_df <- filter(salmonMSE_env$state, x == .env$x, .data$p_smolt == .env$p_naturalsmolt)

    if (nrow(fitness_df)) {
      tmax <- ifelse(nrow(fitness_df), max(fitness_df$t), 1)
      fitness <- filter(fitness_df, t == .env$tmax) %>%
        pull(.data$fitness) %>% unique()
      if (!length(fitness)) fitness <- 1
    } else {
      fitness <- 1
    }

  } else {
    fitness <- 1
  }

  fitness_loss <- fitness^rel_loss
  return(fitness_loss)
}


makeRel_SAR <- function(p_smolt = 1, p_naturalsmolt = p_smolt, fitness_type = c("Ford", "none"), rel_loss, maxage) {

  fitness_type <- match.arg(fitness_type)
  stopifnot(fitness_type == "Ford")

  func <- .SAR_fitness

  #formals(func)$p_smolt <- p_smolt
  formals(func)$p_naturalsmolt <- p_naturalsmolt
  formals(func)$rel_loss <- rel_loss
  formals(func)$fitness_type <- fitness_type

  model <- data.frame(M = 1, N = seq(0, 10), x = -1)

  response <- paste0("M_", p_smolt)
  input <- paste0("N_", p_smolt)
  terms <- c(response, input, "x")

  out <- list(
    func = func,
    model = structure(model, names = terms),
    fitted.values = model$M,
    CV = 0,
    terms = terms,
    type = "Natural mortality",
    Rel = "Marine survival reduced by fitness",
    mult = TRUE,
    age = seq(0, maxage)
  )
  structure(out, class = "SARfitness")
}

#' @export
predict.SARfitness <- predict.RelSmolt

#' @export
simulate.SARfitness <- simulate.RelSmolt
