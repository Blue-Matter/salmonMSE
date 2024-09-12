

calc_zbar <- function(NOS, HOS_effective, NOB, HOB, fec, fec_brood,
                      zbar_prev, omega2, theta, fitness_variance, heritability) {

  pNOS <- NOS * fec/sum(fec * (NOS + HOS_effective))
  pHOSeff <- HOS_effective * fec/sum(fec * (NOS + HOS_effective))

  pNOB <- NOB * fec_brood/sum(fec_brood * (NOB + HOB))
  pHOB <- HOB * fec_brood/sum(fec_brood * (NOB + HOB))

  maxage <- length(NOS)

  # Trait value after selection
  zprime_natural <- (zbar_prev * omega2 + theta[1] * fitness_variance)/(omega2 + fitness_variance)
  zprime_hatchery <- (zbar_prev * omega2 + theta[2] * fitness_variance)/(omega2 + fitness_variance)

  # Change in trait value in the next generation, indexed by environment and brood year
  znext_natural <- zbar_prev + (zprime_natural - zbar_prev) * heritability
  znext_hatchery <- zbar_prev + (zprime_hatchery - zbar_prev) * heritability

  # Realized change after weighting by environment, brood year, and age class fecundity
  z_natural <- sum(znext_natural[, 1] * pNOS, znext_natural[, 2] * pHOSeff)
  z_hatchery <- sum(znext_hatchery[, 1] * pNOB, znext_hatchery[, 2] * pHOB)

  c(z_natural, z_hatchery)
}


calc_fitness <- function(zbar, theta, omega2, fitness_variance, fitness_floor) {
  W <- exp(-0.5 * (zbar - theta)^2/(omega2 + fitness_variance))
  pmax(W, fitness_floor)
}

