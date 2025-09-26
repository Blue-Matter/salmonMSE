
## Internal functions to predict either:
# (b) smolt hatchery production
# (c) smolt natural production
calc_broodtake <- function(NOR_escapement, HOR_escapement, stray, brood_import, ptarget_NOB, pmax_NOB, phatchery, egg_target, p_female,
                           fec_brood, s_prespawn, m) {

  if (is.na(phatchery)) {
    HO_avail <- cbind(
      HOR_escapement,
      stray,
      brood_import
    )
  } else {
    HO_avail <- cbind(
      phatchery * HOR_escapement,
      stray,
      brood_import
    )
  }

  if (sum(HO_avail)) {

    # Optimization can fail for two reasons:
    # (1) not enough unmarked fish
    # (2) more than enough but too many unmarked hatchery fish but can't solve for ptarget_NOB
    opt <- try(
      uniroot(.broodtake_func, c(1e-8, pmax_NOB),
              NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, stray = stray, brood_import = brood_import,
              phatchery = phatchery, p_female = p_female,
              fec = fec_brood, gamma = 1, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m),
      silent = TRUE
    )

    # Assume fail for (1), then set ptake_marked = pmax_NOB
    if (is.character(opt)) {
      ptake_unmarked <- pmax_NOB
    } else {
      ptake_unmarked <- opt$root
    }

    output <- .broodtake_func(
      ptake_unmarked, NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, stray = stray, brood_import = brood_import,
      phatchery = phatchery, p_female = p_female,
      fec = fec_brood, gamma = 1, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
    )

    egg_total <- output$egg_NOB + output$egg_HOB_unmarked + output$egg_HOB_marked

    # Ensure hatchery egg production doesn't exceed target
    # If it does, then optimization failed for reason (2), in which case, we use unmarked fish only to reach production target
    if (egg_total/egg_target > 1.01) {
      unmarked_opt <- try(
        uniroot(.egg_func, c(1e-8, pmax_NOB), N = NOR_escapement + (1 - m) * HOR_escapement,
                gamma = 1, fec = fec_brood, p_female = p_female,
                s_prespawn = s_prespawn, val = egg_target),
        silent = TRUE
      )
      if (is.character(unmarked_opt)) {
        ptake_unmarked <- pmax_NOB
      } else {
        ptake_unmarked <- unmarked_opt$root
      }
      output <- .broodtake_func(
        ptake_unmarked, NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, stray = stray, brood_import = brood_import,
        phatchery = phatchery, p_female = p_female,
        fec = fec_brood, gamma = 1, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
      )
    }

    NOB <- output$NOB
    HOB <- output$HOB_unmarked + output$HOB_marked
    HOB_import <- output$HOB_import
    pNOB <- output$pNOB

  } else {

    NOBopt <- try(
      uniroot(.egg_func, c(1e-8, pmax_NOB), N = NOR_escapement, gamma = 1, fec = fec_brood, p_female = p_female,
              s_prespawn = s_prespawn, val = egg_target),
      silent = TRUE
    )
    if (is.character(NOBopt)) {
      ptake_NOB <- pmax_NOB
    } else {
      ptake_NOB <- NOBopt$root
    }
    NOB <- ptake_NOB * NOR_escapement
    egg_NOB <- sum(NOB * s_prespawn * fec_brood * p_female)
    HOB_zero <- array(0, dim(HOR_escapement))
    HOB_import <- rep(0, nrow(NOB))
    pNOB <- 1

    # Check that output is the same as .broodtake_func()
    output <- list(
      egg_NOB = egg_NOB,
      egg_HOB_unmarked = 0,
      egg_HOB_marked = 0,
      ptake_unmarked = ptake_NOB,
      ptake_marked = 0,
      pNOB = 1,
      NOB = NOB,
      HOB_unmarked = HOB_zero,
      HOB_marked = HOB_zero,
      HOB_import = rep(0, nrow(NOB)),
      HOB_stray = HOB_zero
    )
  }

  return(output)
}

calc_spawners <- function(broodtake, escapement_NOS, escapement_HOS, stray, phatchery, premove_HOS, m) {
  spawners <- list()
  spawners$NOS <- escapement_NOS - broodtake$NOB

  if (sum(stray)) {
    spawners$HOS_stray <- stray - broodtake$HOB_stray # Assume unmarked so not subject to premove_HOS!
  } else {
    spawners$HOS_stray <- array(0, dim(escapement_HOS))
  }
  if (sum(escapement_HOS)) {
    if (is.na(phatchery)) {
      HO <- escapement_HOS - broodtake$HOB_marked - broodtake$HOB_unmarked
    } else {
      HO <- escapement_HOS * (1 - phatchery)
    }
    if (is.function(premove_HOS)) {
      p <- premove_HOS(spawners$NOS, HO, m)
    } else if (is.numeric(premove_HOS)) {
      p <- premove_HOS * m
    } else {
      stop("Error with premove_HOS")
    }
    HOS_local <- HO * (1 - p)
    spawners$HO_remove <- HO * p
  } else {
    HOS_local <- spawners$HO_remove <- array(0, dim(escapement_HOS))
  }
  spawners$HOS <- HOS_local + spawners$HOS_stray
  return(spawners)
}


.egg_func <- function(ptake, N, gamma = 1, fec, p_female, s_prespawn, val = 0) {
  sum(ptake * gamma * N * s_prespawn * fec * p_female) - val
}

calc_broodtake_custom <- function(f_brood, NOR_escapement, HOR_escapement, stray, p_female, fec, s_prespawn, m,
                                  egg_target) {

  brood <- f_brood(NOR_escapement, HOR_escapement, stray, m)
  brood$egg_NOB <- sum(brood$NOB * s_prespawn * fec * p_female)
  brood$egg_HOB_unmarked <- sum((brood$HOB_unmarked + brood$HOB_stray) * s_prespawn * fec * p_female)
  brood$egg_HOB_marked <- sum(brood$HOB_marked * s_prespawn * fec * p_female)
  brood$total_egg <- brood$egg_NOB + brood$egg_HOB_unmarked + brood$egg_HOB_marked

  if (brood$total_egg/egg_target > 1.001) {
    brood_init <- cbind(brood$NOB, brood$HOB_marked, brood$HOB_unmarked, brood$HOB_stray)
    opt <- uniroot(.egg_func, c(1e-8, 1), N = brood_init, fec = fec,
                   p_female = p_female, s_prespawn = s_prespawn, val = egg_target)
    p_adjust <- opt$root

  } else {
    p_adjust <- 1
  }

  NOB <- p_adjust * brood$NOB
  HOB_marked <- p_adjust * brood$HOB_marked
  HOB_unmarked <- p_adjust * brood$HOB_unmarked
  HOB_stray <- p_adjust * brood$HOB_stray

  egg_NOB <- sum(NOB * s_prespawn * fec * p_female)
  egg_HOB_unmarked <- sum((HOB_unmarked + HOB_stray) * s_prespawn * fec * p_female)
  egg_HOB_marked <- sum(HOB_marked * s_prespawn * fec * p_female)

  pNOB <- sum(fec * NOB)/sum(fec * (NOB + HOB_unmarked + HOB_marked + HOB_stray))

  # Check output is the same as .broodtake_func()
  output <- list(
    egg_NOB = egg_NOB,
    egg_HOB_unmarked = egg_HOB_unmarked,
    egg_HOB_marked = egg_HOB_marked,
    ptake_unmarked = sum(NOB, HOB_unmarked, HOB_stray)/sum(NOR_escapement, HOR_escapement, stray),
    ptake_marked = sum(HOB_marked)/sum(NOR_escapement, HOR_escapement, stray),
    pNOB = pNOB,
    NOB = NOB,
    HOB_unmarked = HOB_unmarked,
    HOB_marked = HOB_marked,
    HOB_import = rep(0, nrow(HOR_escapement)),
    HOB_stray = HOB_stray
  )
  return(output)
}

.broodtake_func <- function(ptake_unmarked, NOR_escapement, HOR_escapement, stray, brood_import, phatchery, p_female, fec, gamma,
                            egg_target, s_prespawn, ptarget_NOB, m = 1, opt = TRUE) {

  NOB <- ptake_unmarked * NOR_escapement

  if (is.na(phatchery)) {
    HOB_unmarked <- ptake_unmarked * (1 - m) * HOR_escapement
    HO_avail_marked <- cbind(
      m * HOR_escapement,
      brood_import
    )
  } else {
    HOB_unmarked <- ptake_unmarked * phatchery * (1 - m) * HOR_escapement
    HO_avail_marked <- cbind(
      m * phatchery * HOR_escapement,
      brood_import
    )
  }
  HOB_stray <- ptake_unmarked * stray

  egg_NOB <- sum(NOB * s_prespawn * fec * p_female)
  egg_HOB_unmarked <- sum((HOB_unmarked + HOB_stray) * s_prespawn * fec * p_female)

  egg_HOB_marked <- egg_target - egg_NOB - egg_HOB_unmarked

  if (egg_HOB_marked < 0 || !sum(HO_avail_marked)) {
    ptake_marked <- 0
  } else {
    # solve for the required HOB marked
    HOBopt <- try(
      uniroot(.egg_func, c(1e-8, 1), N = HO_avail_marked, gamma = gamma, fec = fec, p_female = p_female,
              s_prespawn = s_prespawn, val = egg_HOB_marked),
      silent = TRUE
    )
    if (is.character(HOBopt)) {
      ptake_marked <- 1
    } else {
      ptake_marked <- HOBopt$root
    }
  }

  egg_HOB_marked_actual <- .egg_func(
    ptake_marked, N = HO_avail_marked, gamma = gamma, fec = fec, p_female = p_female, s_prespawn = s_prespawn
  )
  HOB_marked <- ptake_marked * HO_avail_marked[, 1:ncol(HOR_escapement), drop = FALSE]
  HOB_import <- ptake_marked * HO_avail_marked[, ncol(HO_avail_marked)]

  if (opt) {
    p_unmarked <- sum(NOB, HOB_unmarked, HOB_stray)/sum(NOB, HOB_unmarked, HOB_marked, HOB_import, HOB_stray)
    obj <- log(p_unmarked) - log(ptarget_NOB)   # Objective function, get close to zero, if not stay positive
    return(obj)
  } else {
    pNOB <- sum(fec * NOB)/sum(fec * (NOB + HOB_unmarked + HOB_marked + HOB_stray), fec * HOB_import)
    output <- list(
      egg_NOB = egg_NOB,
      egg_HOB_unmarked = egg_HOB_unmarked,
      egg_HOB_marked = egg_HOB_marked_actual,
      ptake_unmarked = ptake_unmarked,
      ptake_marked = ptake_marked,
      pNOB = pNOB,
      NOB = NOB,
      HOB_unmarked = HOB_unmarked,
      HOB_marked = HOB_marked,
      HOB_import = HOB_import,
      HOB_stray = HOB_stray
    )
    return(output)
  }
}

calc_yearling <- function(egg_target, s_yearling, s_subyearling, p_yearling, p_subyearling) {

  if (sum(p_yearling) == 1 || sum(p_yearling) == 0) {
    ans <- .yearling_func(
      p_egg_yearling = sum(p_yearling), egg_target = egg_target, s_yearling = s_yearling, s_subyearling = s_subyearling, opt = FALSE
    )
  } else {
    opt <- uniroot(
      .yearling_func, interval = c(1e-8, 0.9999),
      egg_target = egg_target, s_yearling = s_yearling, s_subyearling = s_subyearling, p_yearling = sum(p_yearling)
    )
    ans <- .yearling_func(
      p_yearling = opt$root, egg_target = egg_target, s_yearling = s_yearling, s_subyearling = s_subyearling, opt = FALSE
    )
  }

  yearling <- ans[1] * p_yearling/sum(p_yearling)
  subyearling <- ans[2] * p_subyearling/sum(p_subyearling)

  yearling[is.na(yearling)] <- 0
  subyearling[is.na(subyearling)] <- 0

  list(yearling = yearling, subyearling = subyearling)
}

.yearling_func <- function(p_egg_yearling, egg_target, s_yearling, s_subyearling, p_yearling, opt = TRUE) {
  egg_yearling <- p_egg_yearling * egg_target
  egg_subyearling <- egg_target - egg_yearling
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
#' Calculate smolt production from base stock-recruit parameters and fitness loss
#'
#' @param N1 Egg production for the density-independent component of the stock-recruit relationship. Can be the number of spawners if `phi = 1`
#' and `Smax` is in units of spawners.
#' @param N2 Egg production for the density-dependent component of the stock-recruit relationship (only used if `per_recruit = FALSE`)
#' @param kappa Base productivity parameter
#' @param capacity Base capacity parameter if `SRrel = "BH"`
#' @param Smax Base Smax parameter if `SRrel = "Ricker"`
#' @param phi Unfished egg per smolt (`1/phi` is the replacement line)
#' @param fitness_loss Survival term to reduce smolt production due to fitness, between 0-1
#' @param SRrel Character for the stock-recruit function
#' @param per_recruit Logical, whether N1 is a per recruit quantity (TRUE) or in absolute numbers (FALSE)
#' @return Numeric
#' @export
calc_smolt <- function(N1, N2 = N1, kappa, capacity, Smax, phi = 1, fitness_loss = 1,
                       SRrel = c("BH", "Ricker"), per_recruit = FALSE) {
  SRrel <- match.arg(SRrel)

  a <- kappa/phi * fitness_loss
  if (SRrel == "BH") {
    b <- a/(capacity * fitness_loss)

    if (per_recruit) {
      smolt <- (a * N1 - 1)/b/N1
    } else {
      smolt <- a * N1 / (1 + b * N2)
    }
  } else {
    b <- 1/(Smax * fitness_loss)

    if (per_recruit) {
      smolt <- log(a * N1)/b/N1
    } else {
      smolt <- a * N1 * exp(-b * N2)
    }
  }
  return(smolt)
}

