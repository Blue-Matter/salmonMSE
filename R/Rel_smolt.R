

# Only used if there is a hatchery

## Internal function to predict either:
# (a) "survival" from broodtake of escapement to spawning
# (b) smolt hatchery production
# (c) smolt natural production
# ii = Stock number of adult NOS
# jj = Stock number of adult HOS
# ... = Broodtake hatchery parameters
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


.smolt_func <- function(N, output = c("natural", "hatchery"),
                        ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg, premove_HOS, # Broodtake & hatchery production
                        p_female, fec, gamma, # Spawning (natural production)
                        SRRpars_hist, SRRpars_proj, fitness = "none" # Spawning (natural production)
                        ) {

  output <- match.arg(output)

  if (output == "hatchery") {
    broodtake <- calc_broodtake(N, ptarget_NOB, pmax_NOB, brood_local,
                                fec_brood, s_egg)
    fry <- sum(broodtake) * fec_brood
    smolt <- fry * s_egg
    return(smolt)

  } else if (output == "natural") {

    spawners <- calc_spawners(N, ptarget_NOB, pmax_NOB, brood_local,
                              fec_brood, s_egg, premove_HOS)

    fry_NOS <- spawners[1] * p_female * fec
    fry_HOS <- spawners[2] * p_female * fec * gamma

    total_fry <- fry_NOS + fry_HOS

    if (!total_fry) return(0)

    alpha_hist <- SRRpars_hist["alpha"]
    beta_hist <- SRRpars_hist["beta"]

    alpha_proj <- SRRpars_proj["alpha"]
    beta_proj <- SRRpars_proj["beta"]

    # Predicted smolts from historical SRR parameters and openMSE setup (if there were no hatchery production)
    fry_openMSE <- N[1] * p_female * fec
    smolt_NOS_SRR <- .AHA_SRR(fry_openMSE, fry_openMSE, p = alpha_hist, capacity = alpha_hist/beta_hist)

    # Predicted smolts from projected SRR parameters

    # Still need to add fitness !
    smolt_NOS_proj <- .AHA_SRR(fry_NOS, total_fry, p = alpha_proj, capacity = alpha_proj/beta_proj)

    Perr_y <- smolt_NOS_proj/smolt_NOS_SRR
    return(Perr_y)
  }

}

makeRel_smolt <- function(p_r = 1, p_natural = 2, p_hatchery = 4,
                          output = c("natural", "hatchery"),
                          ptarget_NOB, pmax_NOB, brood_local, fec_brood, s_egg, premove_HOS, # Broodtake & hatchery production
                          p_female, fec, gamma, # Spawning (natural production)
                          SRRpars_hist, SRRpars_proj, fitness = "none" # Spawning (natural production)
                          ) {

  output <- match.arg(output)

  func <- .smolt_func

  formals(func)$ptarget_NOB <- ptarget_NOB
  formals(func)$pmax_NOB <- pmax_NOB
  formals(func)$brood_local <- brood_local
  formals(func)$fec_brood <- fec_brood
  formals(func)$s_egg <- s_egg
  formals(func)$premove_HOS <- premove_HOS
  formals(func)$output <- output

  if (output == "natural") {

    formals(func)$p_female <- p_female
    formals(func)$fec <- fec
    formals(func)$gamma <- gamma

    formals(func)$SRRpars_hist <- SRRpars_hist
    formals(func)$SRRpars_proj <- SRRpars_proj
    formals(func)$fitness <- fitness

    N_natural <- seq(0, 1/SRRpars_proj["beta"], length.out = 10)
    N_hatchery <- seq(0, 1/SRRpars_proj["beta"], length.out = 10)
  } else {
    N_natural <- seq(0, 1000, length.out = 10)
    N_hatchery <- seq(0, 1000, length.out = 10)
  }


  Perr_y <- sapply(1:length(N_natural), function(x) func(c(N_natural[x], N_hatchery[x])))

  model <- data.frame(Perr_y = Perr_y, N_natural = N_natural, N_hatchery = N_hatchery)

  response <- paste0("Perr_y_", p_r)
  input <- paste0("N_", c(p_natural, p_hatchery))
  terms <- c(response, input)

  out <- list(
    func = func,
    model = structure(model, names = terms),
    fitted.values = Perr_y,
    CV = 0,
    terms = terms,
    lag = "next",
    type = "Recruitment deviation",
    Rel = switch(output,
                 "natural" = "Smolt production from escapement",
                 "hatchery" = "Smolt releases from escapement")
  )
  structure(out, class = "RelSmolt")
}


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

  object$func(newdata[vars_Rel] %>% as.numeric())
}

simulate.RelSmolt <- function(object, nsim = 1, seed = 1, ...) {

  set.seed(seed)
  stdev <- sqrt(log(1 + object$CV^2))
  do_sim <- object$fitted.values * rlnorm(nsim * length(object$fitted.values), -0.5 * stdev^2, stdev)
  val <- matrix(do_sim, length(object$fitted.values), nsim) %>% as.data.frame() %>%
    structure(names = paste0("sim_", 1:nsim), seed = seed)
  val

}
