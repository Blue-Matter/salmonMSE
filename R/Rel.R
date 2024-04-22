

# Only used if there is a hatchery

## Internal functions to predict either:
# (b) smolt hatchery production
# (c) smolt natural production
calc_broodtake <- function(N, ptarget_NOB, pmax_NOB, phatchery, brood_local,
                           fec_brood, s_egg) {

  NOR_escapement <- N[1]
  if (length(N) > 1) {
    HOR_escapement <- N[2]

    # NOB is the minimum of the target proportion of natural broodtake (relative to hatchery broodtake) and
    # maximum allowable proportion of the natural escapement
    # brood_local can be stochastic (to be added in the future)
    target_NOB <- ptarget_NOB * brood_local
    max_NOB <- pmax_NOB * NOR_escapement
    NOB <- min(target_NOB, max_NOB)

    # HOB is the minimum of brood_local - NOB or the hatchery origin escapement
    HOB <- min(brood_local - NOB, 0.999 * HOR_escapement * phatchery)
  } else {
    NOB <- NOR_escapement
    HOB <- 0
  }
  c("NOB" = NOB, "HOB" = HOB)
}

calc_spawners <- function(broodtake, escapement, phatchery, premove_HOS) {
  spawners <- structure(rep(0, length(escapement)), names = c("NOS", "HOS"))
  spawners[1] <- escapement[1] - broodtake[1]
  if (length(escapement) > 1) spawners[2] <- escapement[2] * (1 - phatchery) * (1 - premove_HOS)
  return(spawners)
}


.smolt_func <- function(N, x = -1, output = c("natural", "hatchery"),
                        ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg, phatchery, premove_HOS, s_prespawn, # Broodtake & hatchery production
                        p_female, fec, gamma, # Spawning (natural production)
                        SRRpars_hist, SRRpars_proj, SRrel = c("BH", "Ricker"),
                        fitness_type = c("Ford", "none"), # Spawning (natural production)
                        fitness_args, p_smolt) {

  output <- match.arg(output)

  N[is.na(N)] <- 0
  broodtake <- calc_broodtake(N, ptarget_NOB, pmax_NOB, phatchery, brood_local, fec_brood, s_egg)

  if (output == "hatchery") {
    fry <- sum(broodtake) * fec_brood * s_prespawn * p_female
    smolt <- fry * s_egg
    return(smolt)

  } else if (output == "natural") {
    SRrel <- match.arg(SRrel)
    fitness_type <- match.arg(fitness_type)
    spawners <- calc_spawners(broodtake, N, phatchery, premove_HOS)

    NOS <- spawners[1]
    HOS <- spawners[2]
    HOS_effective <- spawners[2] * gamma

    # Fry production in the absence of fitness effects
    fry_NOS <- NOS * p_female * fec
    fry_HOS <- HOS_effective * p_female * fec

    total_fry <- fry_NOS + fry_HOS

    if (!total_fry) return(0) # Perr_y = 0

    alpha_hist <- SRRpars_hist[1, max(x, 1)]
    beta_hist <- SRRpars_hist[2, max(x, 1)]

    alpha_proj <- SRRpars_proj[1, max(x, 1)]
    beta_proj <- SRRpars_proj[2, max(x, 1)]

    # Predicted smolts from historical SRR parameters and openMSE setup (if there were no hatchery production)
    fry_openMSE <- N[1] * p_female * fec
    capacity_SRR <- ifelse(SRrel == "BH", alpha_hist/beta_hist, exp(-1) * alpha_hist/beta_hist)
    smolt_NOS_SRR <- calc_SRR(fry_openMSE, fry_openMSE, p = alpha_hist, capacity = capacity_SRR, SRrel)

    # Predicted fry and smolts from projected SRR parameters and fitness
    if (brood_local > 0 && fitness_type == "Ford" && x > 0) {
      pNOB <- broodtake[1]/sum(broodtake)
      pHOSeff <- HOS_effective/(NOS + HOS_effective)
      pHOScensus <- HOS/(NOS + HOS)

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
      prod_smolt <- alpha_proj * fitness_loss[2]
      capacity_smolt <- fitness_loss[2] * ifelse(SRrel == "BH", alpha_proj/beta_proj, exp(-1) * alpha_proj/beta_proj)

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
      prod_smolt <- alpha_proj
      capacity_smolt <- ifelse(SRrel == "BH", alpha_proj/beta_proj, exp(-1) * alpha_proj/beta_proj)

      fitness <- 1
      pNOB <- 1
      pHOSeff <- 0

      fry_NOS_out <- fry_NOS
      fry_HOS_out <- fry_HOS
    }

    total_fry_out <- fry_NOS_out + fry_HOS_out

    smolt_NOS_proj <- calc_SRR(fry_NOS_out, total_fry_out, prod_smolt, capacity_smolt, SRrel)
    smolt_HOS_proj <- calc_SRR(fry_HOS_out, total_fry_out, prod_smolt, capacity_smolt, SRrel)
    total_smolt <- smolt_NOS_proj + smolt_HOS_proj

    if (x > 0) {
      # Save fry, smolts, and spawners for next generation
      if (nrow(salmonMSE_env$N)) {
        tprev <- local({
          dat <- filter(salmonMSE_env$N, x == .env$x, p_smolt == .env$p_smolt)
          ifelse(nrow(dat), max(dat$t), 0)
        })
      } else {
        tprev <- 0
      }

      df_N <- data.frame(
        x = x,
        p_smolt = p_smolt,
        t = tprev + 1,
        Esc_NOS = N[1],
        Esc_HOS = N[2],
        NOB = broodtake[1],
        HOB = broodtake[2],
        NOS = NOS,
        HOS = HOS,
        HOS_effective = HOS_effective,
        fry_NOS = fry_NOS_out,
        fry_HOS = fry_HOS_out,
        smolt_NOS = smolt_NOS_proj,
        smolt_HOS = smolt_HOS_proj,
        fitness = fitness,
        pNOB = as.numeric(pNOB),
        pHOSeff = as.numeric(pHOSeff),
        pHOScensus = as.numeric(pHOScensus),
        Perr_y = total_smolt/smolt_NOS_SRR
      )
      salmonMSE_env$N <- rbind(salmonMSE_env$N, df_N)
    }

    Perr_y <- total_smolt/smolt_NOS_SRR
    return(Perr_y)
  }
}

makeRel_smolt <- function(p_smolt = 1, p_natural = 2, p_hatchery = 4,
                          output = c("natural", "hatchery"),
                          ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg, phatchery, premove_HOS, s_prespawn, # Broodtake & hatchery production
                          p_female, fec, gamma, # Spawning (natural production)
                          SRRpars_hist, SRRpars_proj, SRrel = c("BH", "Ricker"),  # Spawning (natural production)
                          fitness_type = c("Ford", "none"), fitness_args) {

  output <- match.arg(output)

  func <- .smolt_func

  formals(func)$ptarget_NOB <- ptarget_NOB
  formals(func)$pmax_NOB <- pmax_NOB
  formals(func)$brood_local <- brood_local
  formals(func)$fec_brood <- fec_brood
  formals(func)$s_egg <- s_egg
  formals(func)$phatchery <- phatchery
  formals(func)$premove_HOS <- premove_HOS
  formals(func)$s_prespawn <- s_prespawn
  formals(func)$p_female <- p_female
  formals(func)$output <- output

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
    N_natural <- seq(0, 1/SRRpars_proj[2, 1], length.out = 10)
    N_hatchery <- seq(0, 1/SRRpars_proj[2, 1], length.out = 10)
  } else {
    N_natural <- seq(0, 1000, length.out = 10)
    N_hatchery <- seq(0, 1000, length.out = 10)
  }

  Perr_y <- sapply(1:length(N_natural), function(x) func(c(N_natural[x], N_hatchery[x])))

  model <- data.frame(Perr_y = Perr_y, N_natural = N_natural, N_hatchery = N_hatchery, x = -1)

  response <- paste0("Perr_y_", p_smolt)
  input <- paste0("N_", c(p_natural, p_hatchery))
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

  if(missing(newdata)) {
    return(object$fitted.values)
  }

  vars <- names(newdata)
  vars_Rel <- names(object$model)[-1]
  var_match <- vars_Rel %in% vars

  if (any(!var_match)) {
    stop(paste("Provide a data frame with column names:", paste(vars_Rel, collapse = ", "), "in newdata."))
  }

  vars_data <- sapply(vars, function(x) strsplit(x, "_")[[1]][1])
  vars_call <- unique(vars_data)
  args <- lapply(vars_call, function(x) newdata[vars_data == x] %>% as.numeric()) %>%
    structure(names = vars_call)
  do.call(object$func, args)
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
                         SAR = 0.01, rel_loss = 1, p_naturalsmolt = 1) {

  fitness_type <- match.arg(fitness_type)
  stopifnot(fitness_type == "Ford")

  if (nrow(salmonMSE_env$N) && x > 0) {
    fitness_df <- filter(salmonMSE_env$N, x == .env$x, .data$p_smolt == .env$p_naturalsmolt)

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
  SAR_loss <- SAR[x] * fitness_loss
  M_loss <- -log(SAR_loss)
  return(M_loss)
}


makeRel_SAR <- function(p_smolt = 1, p_naturalsmolt = p_smolt, fitness_type = c("Ford", "none"), SAR, rel_loss, age_mat) {

  fitness_type <- match.arg(fitness_type)
  stopifnot(fitness_type == "Ford")

  func <- .SAR_fitness

  #formals(func)$p_smolt <- p_smolt
  formals(func)$p_naturalsmolt <- p_naturalsmolt
  formals(func)$rel_loss <- rel_loss
  formals(func)$fitness_type <- fitness_type
  formals(func)$SAR <- SAR

  model <- data.frame(M = -log(SAR[1]), N = seq(0, 10), x = -1)

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
    Rel = "Marine survival reduced by fitness"
  )
  if (!missing(age_mat)) out$age <- age_mat - 2
  structure(out, class = "SARfitness")
}

#' @export
predict.SARfitness <- predict.RelSmolt

#' @export
simulate.SARfitness <- simulate.RelSmolt
