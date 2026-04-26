
#' Calculate mean phenotype and fitness
#'
#' Quantitative genetics model of mean phenotypic trait and fitness of next brood year.
#'
#' @param Egg_NOS Numeric, egg production by natural-origin spawners
#' @param Egg_HOS Numeric, egg production by hatchery-origin spawners
#' @param Egg_NOB Numeric, egg production by natural-origin brood
#' @param Egg_HOS Numeric, egg production by hatchery-origin brood
#' @param zbar_prev Numeric, mean phenotype of parents
#' @param fitness_args List, containing `fitness_variance`, `theta`, `phenotype_variance`, `heritability`, `rel_loss`, and
#' `fitness_floor`
#' @returns `fitness_func` returns a list:
#' - `zbar` Numeric length 2, mean phenotype of next generation to natural `zbar[1]` and hatchery `zbar[2]` environments
#' - `fitness` Numeric length 2, fitness of next generation to natural and hatchery environments
#' - `fitness_loss` Matrix `[2, 3]` penalty in survival of next generation due to fitness effects
#' @keywords internal
fitness_func <- function(Egg_NOS, Egg_HOS, Egg_NOB, Egg_HOB, zbar_brood, fitness_args = list()) {

  maxage <- length(Egg_NOS)

  # Column 1 = natural environment, 2 = hatchery environment
  zbar <- rep(NA_real_, 2)
  fitness <- rep(NA_real_, 2)

  # Row 1 = natural environment, 2 = hatchery environment
  # Column = life stages
  fitness_loss <- matrix(1, 2, 3)

  do_fitness <- any(fitness_args$fitness_type == "Ford") && !is.null(fitness_args$heritability) && sum(Egg_NOS, Egg_HOS)

  if (do_fitness) {
    zbar[] <- calc_zbar(
      Egg_NOS, Egg_HOS, Egg_NOB, Egg_HOB, zbar_brood,
      fitness_args$fitness_variance, fitness_args$theta, fitness_args$phenotype_variance, fitness_args$heritability
    )

    for (i in 1:2) {
      if (fitness_args$fitness_type[i] == "Ford") {
        fitness[i] <- calc_fitness(
          zbar[i], fitness_args$theta[i], fitness_args$fitness_variance,
          fitness_args$phenotype_variance, fitness_args$fitness_floor
        )
        fitness_loss[i, ] <- fitness[i]^fitness_args$rel_loss
      }
    }
  }

  output <- list(
    zbar = zbar,
    fitness = fitness,
    fitness_loss = fitness_loss
  )
  return(output)
}

#' @name fitness_func
#' @description `calc_zbar()` is the internal function that calculates mean phenotype.
#' @param fitness_variance Numeric, variance (omega-squared) of the fitness function
#' @param theta Numeric length 2, optimum phenotype value for the natural and hatchery environments, respectively.
#' @param phenotype_variance Numeric, variance (sigma-squared) of the phenotypic trait (theta)
#' @param heritability Numeric, heritability (h-squared) of the phenotypic trait
#' @returns `calc_zbar` returns numeric length 2, mean phenotype of next generation
#' @keywords internal
calc_zbar <- function(Egg_NOS, Egg_HOS, Egg_NOB, Egg_HOB,
                      zbar_prev, fitness_variance, theta, phenotype_variance, heritability) {

  pNOS <- Egg_NOS/sum(Egg_NOS, Egg_HOS)
  if (!sum(Egg_NOS)) pNOS[] <- 0

  pHOSeff <- Egg_HOS/sum(Egg_NOS, Egg_HOS)
  if (!sum(Egg_HOS)) pHOSeff[] <- 0

  pNOB <- Egg_NOB/sum(Egg_NOB, Egg_HOB)
  if (!sum(Egg_NOB)) pNOB[] <- 0
  pHOB <- Egg_HOB/sum(Egg_NOB, Egg_HOB)
  if (!sum(Egg_HOB)) pHOB[] <- 0

  # Trait value after selection
  zprime_natural <- (zbar_prev * fitness_variance + theta[1] * phenotype_variance)/(fitness_variance + phenotype_variance)
  zprime_hatchery <- (zbar_prev * fitness_variance + theta[2] * phenotype_variance)/(fitness_variance + phenotype_variance)

  # Change in trait value in the next generation, indexed by environment and brood year
  znext_natural <- zbar_prev + (zprime_natural - zbar_prev) * heritability
  znext_hatchery <- zbar_prev + (zprime_hatchery - zbar_prev) * heritability

  # Realized change after weighting by environment, brood year, and age class fecundity
  znext_natural[!pNOS, 1] <- znext_natural[!pHOSeff, 2] <- 0
  znext_hatchery[!pNOB, 1] <- znext_hatchery[!pHOB, 2] <- 0

  z_natural <- sum(znext_natural[, 1] * pNOS, znext_natural[, 2] * pHOSeff)
  z_hatchery <- sum(znext_hatchery[, 1] * pNOB, znext_hatchery[, 2] * pHOB)

  c(z_natural, z_hatchery)
}


calc_fitness <- function(zbar, theta, fitness_variance, phenotype_variance, fitness_floor) {
  W <- exp(-0.5 * (zbar - theta)^2/(fitness_variance + phenotype_variance))
  pmax(W, fitness_floor)
}

