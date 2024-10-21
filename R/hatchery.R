
## Internal functions to predict either:
# (b) smolt hatchery production
# (c) smolt natural production
calc_broodtake <- function(Nage, ptarget_NOB, pmax_NOB, phatchery, egg_local, p_female,
                           fec_brood, s_prespawn, m) {

  NOR_escapement <- Nage[, 1]
  if (ncol(Nage) > 1) {
    HOR_escapement <- Nage[, 2]

    # Optimization can fail for two reasons:
    # (1) not enough unmarked fish
    # (2) more than enough but too many unmarked hatchery fish but can't solve for ptarget_NOB
    opt <- try(
      uniroot(.broodtake_func, c(1e-8, pmax_NOB),
              NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
              fec = fec_brood, gamma = 1, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m),
      silent = TRUE
    )

    # Assume fail for (1), then set ptake_marked = pmax_NOB
    if (is.character(opt)) {
      ptake_unmarked <- pmax_NOB
    } else {
      ptake_unmarked <- opt$root
    }

    output <- .broodtake_func(
      ptake_unmarked, NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
      fec = fec_brood, gamma = 1, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
    )

    egg_total <- output$egg_NOB + output$egg_HOB_unmarked + output$egg_HOB_marked

    # Ensure hatchery egg production doesn't exceed target
    # If it does, then optimization failed for reason (2), in which case, we use unmarked fish only to reach production target
    if (egg_total/egg_local > 1.01) {
      unmarked_opt <- try(
        uniroot(.egg_func, c(1e-8, pmax_NOB), N = NOR_escapement + (1 - m) * HOR_escapement,
                gamma = 1, fec = fec_brood, p_female = p_female,
                s_prespawn = s_prespawn, val = egg_local),
        silent = TRUE
      )
      output <- .broodtake_func(
        unmarked_opt$root, NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
        fec = fec_brood, gamma = 1, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
      )
    }

    NOB <- output$NOB
    HOB <- output$HOB_unmarked + output$HOB_marked

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
  if (ncol(escapement) > 1) {
    spawners$HOS <- escapement[, 2] * (1 - phatchery) * (1 - premove_HOS)
  } else {
    spawners$HOS <- rep(0, length(spawners$NOS))
  }
  return(spawners)
}


.egg_func <- function(ptake, N, gamma = 1, fec, p_female, s_prespawn, val = 0) {
  sum(ptake * gamma * N * s_prespawn * fec * p_female) - val
}

.broodtake_func <- function(ptake_NOB, NOR_escapement, HOR_escapement, phatchery, p_female, fec, gamma,
                            egg_target, s_prespawn, ptarget_NOB, m = 1, opt = TRUE) {

  NOB <- ptake_NOB * NOR_escapement
  HOB_unmarked <- ptake_NOB * (1 - m) * HOR_escapement

  egg_NOB <- sum(NOB * s_prespawn * fec * p_female)
  egg_HOB_unmarked <- sum(HOB_unmarked * s_prespawn * fec * p_female)

  egg_HOB_marked <- egg_target - egg_NOB - egg_HOB_unmarked

  if (egg_HOB_marked < 0 || !sum(HOR_escapement)) {
    ptake_marked <- 0
  } else {
    # solve for the required HOB marked
    HOBopt <- try(
      uniroot(.egg_func, c(1e-8, 1), N = m * phatchery * HOR_escapement, gamma = gamma, fec = fec, p_female = p_female,
              s_prespawn = s_prespawn, val = egg_HOB_marked),
      silent = TRUE
    )
    if (is.character(HOBopt)) {
      ptake_marked <- 1
    } else {
      ptake_marked <- HOBopt$root
    }
  }

  HOB_marked <- ptake_marked * m * HOR_escapement * phatchery
  egg_HOB_marked_actual <- .egg_func(
    ptake_marked, N = m * phatchery * HOR_escapement, gamma = gamma, fec = fec, p_female = p_female, s_prespawn = s_prespawn
  )

  if (opt) {
    p <- sum(NOB + HOB_unmarked)/sum(NOB + HOB_unmarked + HOB_marked)
    obj <- log(p) - log(ptarget_NOB)   # Objective function, get close to zero, if not stay positive
    return(obj)
  } else {
    pNOB <- sum(NOB)/sum(NOB + HOB_unmarked + HOB_marked)
    output <- list(
      egg_NOB = egg_NOB,
      egg_HOB_unmarked = egg_HOB_unmarked, egg_HOB_marked = egg_HOB_marked_actual,
      ptake_marked = ptake_marked,
      pNOB = pNOB,
      NOB = NOB, HOB_unmarked = HOB_unmarked, HOB_marked = HOB_marked
    )
    return(output)
  }
}

calc_yearling <- function(egg_local, s_yearling, s_subyearling, p_yearling) {

  if (p_yearling == 1 || p_yearling == 0) {
    ans <- .yearling_func(
      p_egg_yearling = p_yearling, egg_local = egg_local, s_yearling = s_yearling, s_subyearling = s_subyearling, opt = FALSE
    )
  } else {
    opt <- uniroot(
      .yearling_func, interval = c(1e-8, 0.9999),
      egg_local = egg_local, s_yearling = s_yearling, s_subyearling = s_subyearling, p_yearling = p_yearling
    )
    ans <- .yearling_func(
      p_egg_yearling = opt$root, egg_local = egg_local, s_yearling = s_yearling, s_subyearling = s_subyearling, opt = FALSE
    )
  }

  return(ans)
}

.yearling_func <- function(p_egg_yearling, egg_local, s_yearling, s_subyearling, p_yearling, opt = TRUE) {
  egg_yearling <- p_egg_yearling * egg_local
  egg_subyearling <- egg_local - egg_yearling
  yearling <- egg_yearling * s_yearling
  subyearling <- egg_subyearling * s_subyearling

  if (opt) {
    frac <- yearling/(yearling + subyearling)
    return(frac - p_yearling)
  } else {
    return(c(yearling, subyearling))
  }
}

#' Smolt production
#'
#' Calculate smolt production from base stock-recruit parameters, habitat improvement, and fitness loss
#'
#' @param N1 Egg production for the density-independent component of the stock-recruit relationship. Can be the number of spawners if `phi = 1`
#' and `Smax` is in units of spawners.
#' @param N2 Egg production for the density-dependent component of the stock-recruit relationship (only used if `per_recruit = FALSE`)
#' @param kappa Base compensation ratio parameter
#' @param capacity Base capacity parameter if `SRrel = "BH"`
#' @param Smax Base Smax parameter if `SRrel = "Ricker"`
#' @param phi Unfished egg per smolt (`1/phi` is the replacement line)
#' @param kappa_improve Proportional change in kappa from habitat improvement
#' @param capacity_improve Proportional change in capacity or Smax from habitat improvement
#' @param fitness_loss Survival term to reduce smolt production due to fitness, between 0-1
#' @param SRrel Character for the stock-recruit function
#' @param per_recruit Logical, whether N1 is a per recruit quantity (TRUE) or in absolute numbers (FALSE)
#' @return Numeric
#' @export
calc_smolt <- function(N1, N2 = N1, kappa, capacity, Smax, phi = 1, kappa_improve = 1, capacity_improve = 1, fitness_loss = 1,
                       SRrel = c("BH", "Ricker"), per_recruit = FALSE) {
  SRrel <- match.arg(SRrel)

  a <- kappa/phi * kappa_improve * fitness_loss
  if (SRrel == "BH") {
    b <- a/(capacity * capacity_improve * fitness_loss)

    if (per_recruit) {
      smolt <- (a * N1 - 1)/b/N1
    } else {
      smolt <- a * N1 / (1 + b * N2)
    }
  } else {
    b <- 1/(Smax * capacity_improve * fitness_loss)

    if (per_recruit) {
      smolt <- log(a * N1)/b/N1
    } else {
      smolt <- a * N1 * exp(-b * N2)
    }
  }
  return(smolt)
}

