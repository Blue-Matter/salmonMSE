

calc_zbar <- function(NOS, HOS_effective, NOB, HOB, fec, fec_brood,
                      zbar_prev, fitness_variance, theta, phenotype_variance, heritability) {

  tiny <- 1e-8

  pNOS <- NOS * fec/sum(fec * (NOS + HOS_effective))
  pHOSeff <- HOS_effective * fec/sum(fec * (NOS + HOS_effective))

  pNOB <- NOB * fec_brood/sum(fec_brood * (NOB + HOB) + tiny)
  pHOB <- HOB * fec_brood/sum(fec_brood * (NOB + HOB) + tiny)

  # Trait value after selection
  zprime_natural <- (zbar_prev * fitness_variance + theta[1] * phenotype_variance)/(fitness_variance + phenotype_variance)
  zprime_hatchery <- (zbar_prev * fitness_variance + theta[2] * phenotype_variance)/(fitness_variance + phenotype_variance)

  # Change in trait value in the next generation, indexed by environment and brood year
  znext_natural <- zbar_prev + (zprime_natural - zbar_prev) * heritability
  znext_hatchery <- zbar_prev + (zprime_hatchery - zbar_prev) * heritability

  # Realized change after weighting by environment, brood year, and age class fecundity
  z_natural <- sum(znext_natural[, 1] * pNOS, znext_natural[, 2] * pHOSeff)
  z_hatchery <- sum(znext_hatchery[, 1] * pNOB, znext_hatchery[, 2] * pHOB)

  c(z_natural, z_hatchery)
}


calc_fitness <- function(zbar, theta, fitness_variance, phenotype_variance, fitness_floor) {
  W <- exp(-0.5 * (zbar - theta)^2/(fitness_variance + phenotype_variance))
  pmax(W, fitness_floor)
}

