

# Only used if there is a hatchery

## Internal functions to predict either:
# (b) smolt hatchery production
# (c) smolt natural production
calc_broodtake <- function(N, ptarget_NOB, pmax_NOB, brood_local,
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
    HOB <- min(brood_local - NOB, 0.999 * HOR_escapement)
  } else {
    NOB <- NOR_escapement
    HOB <- 0
  }
  c("NOB" = NOB, "HOB" = HOB)
}

calc_spawners <- function(N, ptarget_NOB, pmax_NOB, brood_local,
                          fec_brood, s_egg, premove_HOS) {

  broodtake <- calc_broodtake(N, ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg)
  spawners <- structure(N - broodtake, names = c("NOS", "HOS"))
  spawners[2] <- spawners[2] * premove_HOS

  return(spawners)
}


.smolt_func <- function(N, x, output = c("natural", "hatchery"),
                        ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg, premove_HOS, s_prespawn, # Broodtake & hatchery production
                        p_female, fec, gamma, # Spawning (natural production)
                        SRRpars_hist, SRRpars_proj,
                        fitness_type = c("Ford", "none"), # Spawning (natural production)
                        fitness_args, p_r) {

  output <- match.arg(output)
  broodtake <- calc_broodtake(N, ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg)

  if (output == "hatchery") {
    fry <- sum(broodtake) * fec_brood * s_prespawn * p_female
    smolt <- fry * s_egg
    return(smolt)

  } else if (output == "natural") {
    fitness_type <- match.arg(fitness_type)
    spawners <- structure(N - broodtake, names = c("NOS", "HOS"))
    spawners[2] <- spawners[2] * (1 - premove_HOS)

    NOS <- spawners[1]
    HOS <- spawners[2]
    HOS_effective <- spawners[2] * gamma

    fry_NOS <- NOS * p_female * fec
    fry_HOS <- HOS_effective * p_female * fec

    total_fry <- fry_NOS + fry_HOS

    if (!total_fry) return(0)

    if (missing(x)) {
      alpha_hist <- SRRpars_hist["alpha", 1]
      beta_hist <- SRRpars_hist["beta", 1]

      alpha_proj <- SRRpars_proj["alpha", 1]
      beta_proj <- SRRpars_proj["beta", 1]
    } else {
      alpha_hist <- SRRpars_hist["alpha", x]
      beta_hist <- SRRpars_hist["beta", x]

      alpha_proj <- SRRpars_proj["alpha", x]
      beta_proj <- SRRpars_proj["beta", x]
    }

    # Predicted smolts from historical SRR parameters and openMSE setup (if there were no hatchery production)
    fry_openMSE <- N[1] * p_female * fec
    smolt_NOS_SRR <- .AHA_SRR(fry_openMSE, fry_openMSE, p = alpha_hist, capacity = alpha_hist/beta_hist)

    # Predicted smolts from projected SRR parameters and fitness
    if (fitness_type == "Ford" && !missing(x)) {
      pNOB <- broodtake[1]/sum(broodtake)
      pHOS <- HOS_effective/(NOS + HOS_effective)

      # Get pbar from salmonMSE_env
      if (nrow(salmonMSE_env$Ford)) {
        pbar_prev <- filter(salmonMSE_env$Ford, x == .env$x, p_r == .env$p_r, t == max(.data$t)) %>%
          pull(.data$pbar)
        if (!length(pbar_prev)) pbar_prev <- fitness_args$pbar_start
      } else {
        pbar_prev <- fitness_args$pbar_start
      }

      fitness_calcs <- calc_Ford_fitness(
        pNOB, pHOS, pbar_prev, fitness_args$theta, fitness_args$omega2, fitness_args$fitness_variance,
        fitness_args$A, fitness_args$fitness_floor, fitness_args$heritability
      )
      pbar <- fitness_calcs$pbar
      fitness <- fitness_calcs$fitness

      fitness_loss <- fitness^fitness_args$rel_loss

      prod_smolt <- alpha_proj * prod(fitness_loss[1:2])
      capacity_smolt <- alpha_proj/beta_proj * prod(fitness_loss[1:2])

      # Save pbar and fitness for next generation
      if (nrow(salmonMSE_env$Ford)) {
        t <- filter(salmonMSE_env$Ford, x == .env$x, p_r == .env$p_r, t == max(.data$t)) %>%
          pull(.data$t) %>% unique()
        if (!length(t)) t <- 1
      } else {
        t <- 1
      }
      df_Ford <- data.frame(
        x = x,
        p_r = p_r,
        t = t + 1,
        pbar = pbar,
        fitness = fitness,
        pNOB = as.numeric(pNOB),
        pHOS = as.numeric(pHOS)
      )

      salmonMSE_env$Ford <- rbind(salmonMSE_env$Ford, df_Ford)

    } else { #if (fitness_type == "none") {
      prod_smolt <- alpha_proj
      capacity_smolt <- alpha_proj/beta_proj
    }

    smolt_NOS_proj <- .AHA_SRR(fry_NOS, total_fry, p = prod_smolt, capacity = capacity_smolt)
    smolt_HOS_proj <- .AHA_SRR(fry_HOS, total_fry, p = prod_smolt, capacity = capacity_smolt)
    total_smolt <- smolt_NOS_proj + smolt_HOS_proj

    if (!missing(x)) {
      # Save fry, smolts, and spawners for next generation
      if (nrow(salmonMSE_env$N)) {
        t <- filter(salmonMSE_env$N, x == .env$x, p_r == .env$p_r, t == max(.data$t)) %>%
          pull(.data$t) %>% unique()
        if (!length(t)) t <- 1
      } else {
        t <- 1
      }

      df_N <- data.frame(
        x = x,
        p_r = p_r,
        t = t + 1,
        NOS = NOS,
        HOS = HOS,
        HOS_effective = HOS_effective,
        fry_NOS = fry_NOS,
        fry_HOS = fry_HOS,
        smolt_NOS = smolt_NOS_proj,
        smolt_HOS = smolt_HOS_proj
      )
      salmonMSE_env$N <- rbind(salmonMSE_env$N, df_N)
    }

    Perr_y <- total_smolt/smolt_NOS_SRR
    return(Perr_y)
  }
}

makeRel_smolt <- function(p_r = 1, p_natural = 2, p_hatchery = 4,
                          output = c("natural", "hatchery"),
                          ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg, premove_HOS, s_prespawn, # Broodtake & hatchery production
                          p_female, fec, gamma, # Spawning (natural production)
                          SRRpars_hist, SRRpars_proj,  # Spawning (natural production)
                          fitness_type = c("Ford", "none"), fitness_args) {

  output <- match.arg(output)

  func <- .smolt_func

  formals(func)$ptarget_NOB <- ptarget_NOB
  formals(func)$pmax_NOB <- pmax_NOB
  formals(func)$brood_local <- brood_local
  formals(func)$fec_brood <- fec_brood
  formals(func)$s_egg <- s_egg
  formals(func)$premove_HOS <- premove_HOS
  formals(func)$s_prespawn <- s_prespawn
  formals(func)$p_female <- p_female
  formals(func)$output <- output

  if (output == "natural") {

    formals(func)$p_female <- p_female
    formals(func)$fec <- fec
    formals(func)$gamma <- gamma

    formals(func)$SRRpars_hist <- SRRpars_hist
    formals(func)$SRRpars_proj <- SRRpars_proj
    formals(func)$fitness_type <- fitness_type

    if (fitness_type != "none") {
      formals(func)$fitness_args <- fitness_args
      formals(func)$p_r <- p_r
    }
    N_natural <- seq(0, 1/SRRpars_proj["beta", 1], length.out = 10)
    N_hatchery <- seq(0, 1/SRRpars_proj["beta", 1], length.out = 10)
  } else {
    N_natural <- seq(0, 1000, length.out = 10)
    N_hatchery <- seq(0, 1000, length.out = 10)
  }

  Perr_y <- sapply(1:length(N_natural), function(x) func(c(N_natural[x], N_hatchery[x])))

  model <- data.frame(Perr_y = Perr_y, N_natural = N_natural, N_hatchery = N_hatchery, x = 1)

  response <- paste0("Perr_y_", p_r)
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
.SAR_fitness <- function(N, x = 1,
                         fitness_type = c("Ford", "none"), # Spawning (natural production)
                         SAR = 0.01, rel_loss = 1, p_r = 1) {

  fitness_type <- match.arg(fitness_type)
  stopifnot(fitness_type == "Ford")

  if (nrow(salmonMSE_env$Ford)) {
    fitness <- filter(salmonMSE_env$Ford, x == .env$x, p_r == .env$p_r, t == max(.data$t)) %>%
      pull(.data$fitness) %>% unique()
    if (!length(fitness)) fitness <- 1
  } else {
    fitness <- 1
  }

  fitness_loss <- fitness^rel_loss
  SAR_loss <- SAR[x] * fitness_loss
  M_loss <- -log(SAR_loss)
  return(M_loss)
}


makeRel_SAR <- function(p_r = 1, fitness_type = c("Ford", "none"), SAR, rel_loss, age_mat) {

  fitness_type <- match.arg(fitness_type)
  stopifnot(fitness_type == "Ford")

  func <- .SAR_fitness

  formals(func)$p_r <- p_r
  formals(func)$rel_loss <- rel_loss
  formals(func)$fitness_type <- fitness_type
  formals(func)$SAR <- SAR

  model <- data.frame(M = -log(SAR[1]), N = seq(0, 10), x = 1)

  response <- paste0("M_", p_r)
  input <- paste0("N_", p_r)
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
