
#' Solve for broodtake numbers
#'
#' Internal functions to calculate broodtake with various constraints, called by [brood_func()].
#' `calc_broodtake()` is a wrapper function.
#'
#' @param NO Matrix `[nage, n_g]`, natural-origin fish available for broodtake
#' @param HO Matrix `[nage, n_r]`, hatchery-origin fish, potentially available for broodtake
#' @param stray Matrix `[nage, n_r]`, hatchery-origin strays available for broodtake
#' @param brood_import Vector `[nage]` of imported brood
#' @param ptarget_NOB Numeric, target proportion of `NOB/(NOB + HOB)`. If `m = 1`, then the realized pNOB should be ptarget_NOB.
#' If `m < 1`, then the system achieves `unmarked brood/(NOB + HOB)  = ptarget_NOB`.
#' @param pmax_NOB Numeric, maximum proportion of `NOB/NO`
#' @param phatchery Numeric, proportion of `HO` that return to hatchery instead of the spawning ground. HOB is taken from this subset. Set to `NA` to
#' obtain HOB on the way to spawning ground.
#' @param egg_target Numeric, target egg production from which to back-calculate brood numbers
#' @param p_female Vector `[nage]` of proportion female to calculate egg production
#' @param fec Vector `[nage]`, egg production per female
#' @param s_prespawn Numeric, survival of brood prior to egg production
#' @param m Numeric, mark rate of `HO` used. Discounts `ptarget_NOB` and `pmax_NOB` based on marked fish
#' @returns
#' `calc_broodtake()` returns a list from `.broodtake_func()`
calc_broodtake <- function(NO, HO, stray, brood_import, ptarget_NOB, pmax_NOB,
                           phatchery, egg_target, p_female,
                           fec, s_prespawn, m) {

  if (is.null(phatchery)) phatchery <- NA

  if (is.na(phatchery)) { # Hatchery-origin broodtake on the way to spawning ground
    HO_avail <- cbind(
      HO,
      stray,
      brood_import
    )
  } else {  # Hatchery-origin broodtake taken from fish that return to hatchery (seining, swim-in facilities, etc.)
    HO_avail <- cbind(
      phatchery * HO,
      stray,
      brood_import
    )
  }

  if (sum(HO_avail)) {

    # Optimization for pNOB target, but can fail for two reasons:
    # (1) not enough unmarked fish
    # (2) more than enough but too many unmarked hatchery fish but can't solve for ptarget_NOB
    opt <- try(
      uniroot(.broodtake_func, c(.Machine$double.eps, pmax_NOB),
              tol = .Machine$double.eps,
              NO = NO, HO = HO, stray = stray, brood_import = brood_import,
              phatchery = phatchery, p_female = p_female,
              fec = fec, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m),
      silent = TRUE
    )

    # Assume fail for (1), then set ptake_marked = pmax_NOB
    if (is.character(opt)) {
      ptake_unmarked <- pmax_NOB
    } else {
      ptake_unmarked <- opt$root
    }

    output <- .broodtake_func(
      ptake_unmarked, NO = NO, HO = HO, stray = stray, brood_import = brood_import,
      phatchery = phatchery, p_female = p_female,
      fec = fec, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
    )

    egg_total <- sum(output$egg_NOB, output$egg_HOB_unmarked, output$egg_HOB_marked)

    # Ensure hatchery egg production doesn't exceed target
    # If it does, then optimization failed for reason (2), in which case, we use unmarked fish only to reach production target
    if (egg_total/egg_target > 1.01) {
      unmarked_opt <- try(
        uniroot(.egg_func, c(.Machine$double.eps, pmax_NOB),
                tol = .Machine$double.eps,
                N = NO + (1 - m) * HO,
                gamma = 1, fec = fec, p_female = p_female,
                s_prespawn = s_prespawn, val = egg_target),
        silent = TRUE
      )
      if (is.character(unmarked_opt)) {
        ptake_unmarked <- pmax_NOB
      } else {
        ptake_unmarked <- unmarked_opt$root
      }
      output <- .broodtake_func(
        ptake_unmarked, NO = NO, HO = HO, stray = stray, brood_import = brood_import,
        phatchery = phatchery, p_female = p_female,
        fec = fec, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
      )
    }

  } else {

    if (sum(NO) && egg_target > 0) {

      NOBopt <- try(
        uniroot(.egg_func, c(1e-8, pmax_NOB), N = NO, gamma = 1, fec = fec, p_female = p_female,
                s_prespawn = s_prespawn, val = egg_target),
        silent = TRUE
      )
      if (is.character(NOBopt)) {
        ptake_NOB <- pmax_NOB
      } else {
        ptake_NOB <- NOBopt$root
      }
    } else {
      ptake_NOB <- 0
    }

    output <- .broodtake_func(
      ptake_NOB, NO = NO, HO = HO, stray = stray, brood_import = brood_import,
      phatchery = phatchery, p_female = p_female,
      fec = fec, egg_target = egg_target, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, m = m, opt = FALSE
    )
  }

  return(output)
}

#' Calculate spawners after broodtake
#'
#' Internal function that calculates remaining spawners after broodtake and additional removals.
#'
#' @param broodtake List, returned by [calc_broodtake()]
#' @param NO Matrix `[nage, n_g]`, natural-origin in-river return
#' @param HO Matrix `[nage, n_r]`, hatchery-origin in-river return
#' @param stray Matrix `[nage, n_r]`, hatchery-origin strays in-river return
#' @param phatchery Numeric, proportion of `HO` that return to hatchery instead of the spawning ground. HOB is taken from this subset. Set to `NA` to
#' obtain HOB on the way to spawning ground.
#' @param premove_HOS Numeric, proportion of HO fish to remove (after brood). Can also be a function
#' @param premove_NOS Numeric, proportion of NO fish to remove (after brood). Can also be a function
#' @param m Numeric, mark rate
#' @returns
#' Named list:
#' - NOS `[nage, n_g]` natural-origin spawners
#' - HOS `[nage, n_r]` hatchery-origin spawners
#' - HOS_stray `[nage, n_r]` spawners that are strays
#' - NO_remove `[nage, n_g]` natural-origin fish removed before spawning
#' - HO_remove `[nage, n_r]` hatchery-origin fish removed before spawning
#' @keywords internal
calc_spawners <- function(broodtake, NO, HO, stray, phatchery, premove_HOS, premove_NOS, m) {
  if (is.null(premove_HOS)) premove_HOS <- 0
  if (is.null(premove_NOS)) premove_NOS <- 0
  if (is.null(m)) m <- 0

  spawners <- list()

  NO_after_brood <- NO - broodtake$NOB

  if (sum(stray)) {
    spawners$HOS_stray <- stray - broodtake$HOB_stray # Assume unmarked so not subject to premove_HOS!
  } else {
    spawners$HOS_stray <- array(0, dim(HO))
  }

  if (sum(HO)) {
    if (is.na(phatchery)) {
      HO_after_brood <- HO - broodtake$HOB_marked - broodtake$HOB_unmarked
    } else {
      HO_after_brood <- HO * (1 - phatchery)
    }
    if (is.function(premove_HOS)) {
      pHO <- premove_HOS(NO_after_brood, HO_after_brood, m)
    } else if (is.numeric(premove_HOS)) {
      pHO <- premove_HOS * m
    } else {
      stop("Error with premove_HOS")
    }
    HOS_local <- HO_after_brood * (1 - pHO)
    spawners$HO_remove <- HO_after_brood * pHO
  } else {
    HOS_local <- spawners$HO_remove <- HO_after_brood <- array(0, dim(HO))
  }
  spawners$HOS <- HOS_local + spawners$HOS_stray

  if (sum(NO)) {

    if (is.function(premove_NOS)) {
      pNO <- premove_NOS(NO_after_brood, HO_after_brood, m)
    } else if (is.numeric(premove_NOS)) {
      pNO <- premove_NOS * (1 - m)
    } else {
      stop("Error with premove_NOS")
    }

    spawners$NOS <- NO_after_brood * (1 - pNO)
    spawners$NO_remove <- NO_after_brood * pNO
  } else {
    spawners$NOS <- spawners$NO_remove <- array(0, dim(NO))
  }

  return(spawners)
}

#' Egg production function
#'
#' Simple wrapper function to calculate hatchery egg production
#'
#' @param ptake Numeric, proportion of spawners that spawn
#' @param N Numeric, spawners
#' @param gamma Numeric, relative reproductive success of spawners
#' @param p_female Numeric, proportion female
#' @param s_prespawn Numeric, survival of spawners prior to egg production
#' @param val Numeric, target egg production. Used to optimize for `ptake` if `opt = TRUE`
#' @param opt Logical, whether the function is used to optimize for `ptake`
#' @returns Numeric
#' @keywords internal
.egg_func <- function(ptake = 1, N, gamma = 1, fec, p_female, s_prespawn, val = 0, opt = TRUE) {
  egg <- ptake * gamma * N * s_prespawn * fec * p_female
  if (opt) {
    return(sum(egg) - val)
  } else {
    return(egg)
  }
}

#' @name calc_broodtake
#' @description `calc_broodtake_custom()` uses a user-provided function to generate brood numbers and adjusts
#' downwards if egg production exceeds the target.
#' @param f_brood Function that calculates the brood numbers
#' @returns
#' `calc_broodtake_custom()` returns a named list, same format as `calc_broodtake()`
#' @keywords internal
calc_broodtake_custom <- function(f_brood, NO, HO, stray, p_female, fec, s_prespawn, m,
                                  egg_target) {

  brood <- f_brood(NO, HO, stray, m)
  brood$egg_NOB <- brood$NOB * s_prespawn * fec * p_female
  brood$egg_HOB_unmarked <- (brood$HOB_unmarked + brood$HOB_stray) * s_prespawn * fec * p_female
  brood$egg_HOB_marked <- brood$HOB_marked * s_prespawn * fec * p_female

  total_egg <- sum(brood$egg_NOB, brood$egg_HOB_unmarked, brood$egg_HOB_marked)

  if (total_egg/egg_target > 1.001) {
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

  egg_NOB <- NOB * s_prespawn * fec * p_female
  egg_HOB_unmarked <- (HOB_unmarked + HOB_stray) * s_prespawn * fec * p_female
  egg_HOB_marked <- HOB_marked * s_prespawn * fec * p_female

  pNOB <- sum(fec * NOB)/sum(fec * NOB, fec * HOB_unmarked, fec * HOB_marked, fec * HOB_stray)

  # Check output is the same as .broodtake_func()
  output <- list(
    egg_NOB = egg_NOB,
    egg_HOB_unmarked = egg_HOB_unmarked,
    egg_HOB_marked = egg_HOB_marked,
    ptake_unmarked = sum(NOB, HOB_unmarked, HOB_stray)/sum(NO, HO, stray),
    ptake_marked = sum(HOB_marked)/sum(NO, HO, stray),
    pNOB = pNOB,
    NOB = NOB,
    HOB_unmarked = HOB_unmarked,
    HOB_marked = HOB_marked,
    HOB_import = rep(0, nrow(HO)),
    HOB_stray = HOB_stray
  )
  return(output)
}


#' @name calc_broodtake
#' @param ptake_unmarked Numeric, proportion of unmarked fish used for brood
#' @param opt Logical, whether the function is used for optimization (TRUE) or reporting (FALSE)
#' @description `.broodtake_func()` is the optimization function used to calculate broodtake
#' @returns
#' `.broodtake_func()` returns a numeric if `opt = TRUE`: `log(p_unmarked) - log(ptarget_NOB)`, otherwise
#' a named list of brood and egg production.
#'
#' - `egg_NOB` Matrix `[nage, n_g]`
#' - `egg_HOB_unmarked` Matrix `[nage, n_r]` (including strays)
#' - `egg_HOB_marked` Matrix `[nage, n_r]`
#' - `egg_HOB_import` Vector `[nage]`
#' - `ptake_unmarked` Numeric
#' - `ptake_marked` Numeric
#' - `pNOB` Numeric
#' - `NOB` Numeric
#' - `HOB_unmarked` Matrix `[nage, n_r]`
#' - `HOB_marked` Matrix `[nage, n_r]`
#' - `HOB_import` Matrix `[nage, n_r]`
#' - `HOB_stray` Vector `[nage]`
#' @keywords internal
.broodtake_func <- function(ptake_unmarked, NO, HO, stray, brood_import, phatchery, p_female, fec,
                            egg_target, s_prespawn, ptarget_NOB, m = 1, opt = TRUE) {

  NOB <- ptake_unmarked * NO
  if (is.na(phatchery)) { # Hatchery-origin broodtake on the way to spawning ground
    HOB_unmarked <- ptake_unmarked * (1 - m) * HO
    HO_avail_marked <- cbind(
      m * HO,
      brood_import
    )
  } else { # Hatchery-origin broodtake taken from fish that return to hatchery (seining, swim-in facilities, etc.)
    HOB_unmarked <- ptake_unmarked * phatchery * (1 - m) * HO
    HO_avail_marked <- cbind(
      m * phatchery * HO,
      brood_import
    )
  }
  HOB_stray <- ptake_unmarked * stray

  if (sum(NOB)) {
    egg_NOB <- NOB * s_prespawn * fec * p_female
  } else {
    egg_NOB <- array(0, dim(NO))
  }

  if (sum(HOB_unmarked, HOB_stray)) {
    egg_HOB_unmarked <- (HOB_unmarked + HOB_stray) * s_prespawn * fec * p_female
  } else {
    egg_HOB_unmarked <- array(0, dim(HO))
  }

  egg_HOB_marked <- egg_target - sum(egg_NOB, egg_HOB_unmarked)

  if (egg_HOB_marked <= 0 || !sum(HO_avail_marked)) {
    ptake_marked <- 0
    HOB_marked <- array(0, dim(HO))
    HOB_import <- rep(0, nrow(HO))

    egg_HOB_marked_actual <- cbind(
      HOB_marked,
      HOB_import
    )
  } else {
    # solve for the required HOB marked
    HOBopt <- try(
      uniroot(.egg_func, c(1e-8, 1), N = HO_avail_marked, gamma = 1, fec = fec, p_female = p_female,
              s_prespawn = s_prespawn, val = egg_HOB_marked),
      silent = TRUE
    )
    if (is.character(HOBopt)) {
      ptake_marked <- 1
    } else {
      ptake_marked <- HOBopt$root
    }

    egg_HOB_marked_actual <- .egg_func(
      ptake_marked, N = HO_avail_marked, gamma = 1, fec = fec, p_female = p_female, s_prespawn = s_prespawn,
      opt = FALSE
    )
    HOB_marked <- ptake_marked * HO_avail_marked[, 1:ncol(HO), drop = FALSE]
    HOB_import <- ptake_marked * HO_avail_marked[, ncol(HO_avail_marked)]
  }

  if (opt) {
    p_unmarked <- sum(NOB, HOB_unmarked, HOB_stray)/sum(NOB, HOB_unmarked, HOB_marked, HOB_import, HOB_stray)
    if (is.na(p_unmarked)) p_unmarked <- 0
    obj <- log(p_unmarked) - log(ptarget_NOB)   # Objective function, get close to zero, if not stay positive
    return(obj)
  } else {

    pNOB <- sum(fec * NOB)/sum(fec * NOB, fec * HOB_unmarked, fec * HOB_marked, fec * HOB_stray, fec * HOB_import)
    output <- list(
      egg_NOB = egg_NOB,
      egg_HOB_unmarked = egg_HOB_unmarked, # Include strays
      egg_HOB_marked = egg_HOB_marked_actual[, -ncol(egg_HOB_marked_actual), drop = FALSE],
      egg_HOB_import = egg_HOB_marked_actual[, ncol(egg_HOB_marked_actual)],
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

#' Calculate hatchery releases
#'
#' From egg target, calculate production of yearlings and subyearlings conditional on survival and
#' proportion of releases.
#'
#' @param egg_target Numeric
#' @param s_yearling Numeric, survival from egg to yearling
#' @param s_subyearling Numeric, survival from egg to subyearling
#' @param p_yearling Vector `[n_r]`, proportion of releases at yearling stage, where `sum(p_yearling, p_subyearling) = 1`
#' @param p_yearling Vector `[n_r]`, proportion of releases at subyearling stage, where `sum(p_yearling, p_subyearling) = 1`
#' @returns
#' `calc_yearling()` returns a named list returned by `.yearling_func()`
#' @keywords internal
calc_yearling <- function(egg_target, s_yearling, s_subyearling, p_yearling, p_subyearling) {

  if (egg_target > 0) {
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

    yearling <- sum(ans$yearling) * p_yearling/sum(p_yearling)
    subyearling <- sum(ans$subyearling) * p_subyearling/sum(p_subyearling)

    yearling[is.na(yearling)] <- 0
    subyearling[is.na(subyearling)] <- 0
  } else {
    yearling <- subyearling <- 0
  }

  list(yearling = yearling, subyearling = subyearling)
}

#' @name calc_yearling
#' @description `.yearling_func()` is the optimization function that determines proportions.
#' @param p_egg_yearling Numeric, proportion to eggs to become yearlings
#' @param opt Logical, whether the function is used for optimization (TRUE) or reporting (FALSE)
#' @keywords internal
#' @returns
#' `.yearling_func()` returns a numeric if `opt = TRUE`: `yearling/(yearling + subyearling) - p_yearling`.
#' Otherwise, log(p_unmarked) - log(ptarget_NOB)`, otherwise a named list:
#' - `yearling` Vector `[n_r]` of yearling releases
#' - `subyearling` Vector `[n_r]` of subyearling releases
.yearling_func <- function(p_egg_yearling, egg_target, s_yearling, s_subyearling, p_yearling, opt = TRUE) {
  egg_yearling <- p_egg_yearling * egg_target
  egg_subyearling <- egg_target - egg_yearling
  yearling <- egg_yearling * s_yearling
  subyearling <- egg_subyearling * s_subyearling

  if (opt) {
    frac <- yearling/(yearling + subyearling)
    return(frac - p_yearling)
  } else {
    return(list(yearling = yearling, subyearling = subyearling))
  }
}


#' Proportion wild spawners
#'
#' @description Calculate the proportion of wild spawners from a time series of spawners
#' - `calc_pwild()` is the simple calculation based on the proportion of hatchery spawners
#' - `calc_pwild_age()` performs the calculation weighted by age class fecundity
#'
#' @param NOS_a Array `[nsim, maxage, years]` for natural spawners
#' @param HOS_a Array `[nsim, maxage, years]` for hatchery spawners
#' @param fec Array `[nsim, maxage, years]` for age class fecundity
#' @param gamma Numeric, reduced reproductive success of hatchery spawners
#' @param pHOS_cur Numeric, proportion of hatchery spawners in current generation
#' @param pHOS_prev Numeric, proportion of hatchery spawners in previous generation
#' @return `calc_pwild_age()` a matrix of pWILD by simulation and year. `calc_pwild()` returns a numeric
#' @keywords internal
calc_pwild_age <- function(NOS_a, HOS_a, fec, gamma) {

  nsim <- dim(NOS_a)[1]
  maxage <- dim(NOS_a)[2]
  proyears <- dim(NOS_a)[3]

  prob <- array(NA_real_, c(nsim, proyears))

  for (y in 1:proyears) {
    if (sum(NOS_a[, , y], HOS_a[, , y])) {
      pNOS_y <- matrix(1, nsim, maxage)
      pnatural_b <- phetero_b <- phatch_b <- matrix(0, nsim, maxage)

      pNOS_y[] <- NOS_a[, , y]/rowSums(NOS_a[, , y] + HOS_a[, , y])

      for (a in 1:maxage) {
        pNOS_b <- matrix(1, nsim, maxage)
        pHOS_b <- matrix(0, nsim, maxage)

        b <- y - a # brood_year
        if (b > 0) {
          pNOS_b[] <- NOS_a[, , b]/rowSums(fec[, , b] * NOS_a[, , b] + fec[, , b] * HOS_a[, , b])
          pHOS_b[] <- HOS_a[, , b]/rowSums(fec[, , b] * NOS_a[, , b] + fec[, , b] * HOS_a[, , b])
        }
        pnatural_b[, a] <- rowSums(pNOS_b, na.rm = TRUE)^2
        phetero_b[, a] <- 2 * gamma * rowSums(pNOS_b, na.rm = TRUE) * rowSums(pHOS_b, na.rm = TRUE)
        phatch_b[, a] <- gamma^2 * rowSums(pHOS_b, na.rm = TRUE)^2
      }
      denom <- pnatural_b + phetero_b + phatch_b

      prob[, y] <- rowSums(pNOS_y * pnatural_b/denom, na.rm = TRUE)
    }
  }

  return(prob)
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
#' @param phi Unfished egg per smolt (`1/phi` corresponds to the one-to-one adult/spawner replacement line)
#' @param tau Unfished spawner per smolt
#' @param fitness_loss Survival term to reduce smolt production due to fitness, between 0-1
#' @param SRrel Character for the stock-recruit function
#' @param per_recruit Logical, whether N1 is a per recruit quantity (TRUE) or in absolute numbers (FALSE)
#' @return Numeric
#' @export
calc_smolt <- function(N1, N2 = N1, kappa, capacity, Smax, phi = 1, tau = 1, fitness_loss = 1,
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
    Smax_prime <- Smax * fitness_loss
    Emax <- Smax_prime * phi/tau
    b <- 1/Emax

    if (per_recruit) {
      smolt <- log(a * N1)/b/N1
    } else {
      smolt <- a * N1 * exp(-b * N2)
    }
  }
  return(smolt)
}

