
## Internal functions to predict either:
# (b) smolt hatchery production
# (c) smolt natural production
calc_broodtake <- function(Nage, ptarget_NOB, pmax_NOB, phatchery, egg_local, p_female,
                           fec_brood, s_prespawn) {

  NOR_escapement <- Nage[, 1]
  if (ncol(Nage) > 1) {
    HOR_escapement <- Nage[, 2]

    opt <- try(
      uniroot(.broodtake_func, c(1e-8, pmax_NOB),
              NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
              fec = fec_brood, gamma = 1, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB),
      silent = TRUE
    )
    if (is.character(opt)) {
      ptake_NOB <- pmax_NOB
    } else {
      ptake_NOB <- opt$root
    }

    output <- .broodtake_func(
      ptake_NOB, NOR_escapement = NOR_escapement, HOR_escapement = HOR_escapement, phatchery = phatchery, p_female = p_female,
      fec = fec_brood, gamma = 1, egg_target = egg_local, s_prespawn = s_prespawn, ptarget_NOB = ptarget_NOB, opt = FALSE
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
  egg_HOB_actual <- .egg_func(
    ptake_HOB, N = HOR_escapement * phatchery, gamma = gamma, fec = fec, p_female = p_female, s_prespawn = s_prespawn
  )

  if (opt) {
    obj <- log(sum(NOB)/sum(NOB + HOB)) - log(ptarget_NOB)   # Objective function, get close to zero, if not stay positive
    return(obj)
  } else {
    output <- list(egg_NOB = egg_NOB, egg_HOB = egg_HOB_actual, ptake_HOB = ptake_HOB, NOB = NOB, HOB = HOB)
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
