

calc_zbar <- function(NOS, HOS_effective, NOB, HOB, fec, fec_brood,
                      zbar_prev, zbar_start, y, omega2, theta, fitness_variance, heritability) {

  pNOS <- NOS * fec/sum(fec * (NOS + HOS_effective))
  pHOSeff <- HOS_effective * fec/sum(fec * (NOS + HOS_effective))

  pNOB <- NOB * fec_brood/sum(fec_brood * (NOB + HOB))
  pHOB <- HOB * fec_brood/sum(fec_brood * (NOB + HOB))

  maxage <- length(NOS)
  zbar <- matrix(0, maxage, 2) # Column 1 = natural environment, 2 = hatchery environment

  # Trait value by brood year
  for (a in 1:maxage) {
    if (NOS[a] || NOB[a]) {
      zbar1 <- dplyr::filter(zbar_prev, .data$t == .env$y - 2 * .env$a, .data$type == "natural") %>%
        pull(.data$zbar)
      if (!length(zbar1)) zbar1 <- zbar_start[1]
      zbar[a, 1] <- zbar1
    }
    if (HOS_effective[a] || HOB[a]) {
      zbar2 <- dplyr::filter(zbar_prev, .data$t == .env$y - 2 * .env$a, .data$type == "hatchery") %>%
        pull(.data$zbar)
      if (!length(zbar2)) zbar2 <- zbar_start[2]
      zbar[a, 2] <- zbar2
    }
  }

  # Trait value after selection
  zprime_natural <- (zbar * omega2 + theta[1] * fitness_variance)/(omega2 + fitness_variance)
  zprime_hatchery <- (zbar * omega2 + theta[2] * fitness_variance)/(omega2 + fitness_variance)

  # Change in trait value in the next generation, indexed by environment and brood year
  znext_natural <- zbar + (zprime_natural - zbar) * heritability
  znext_hatchery <- zbar + (zprime_hatchery - zbar) * heritability

  # Realized change after weighting by environment, brood year, and age class fecundity
  z_natural <- sum(znext_natural[, 1] * pNOS, znext_natural[, 2] * pHOSeff)
  z_hatchery <- sum(znext_hatchery[, 1] * pNOB, znext_hatchery[, 2] * pHOB)

  c(z_natural, z_hatchery)
}


calc_fitness <- function(zbar, theta, omega2, fitness_variance, fitness_floor) {
  W <- exp(-0.5 * (zbar - theta)^2/(omega2 + fitness_variance))
  pmax(W, fitness_floor)
}

