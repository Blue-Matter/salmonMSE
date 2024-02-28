
# pbar Vector length 2. Pbar value from the previous generation.
# First value is for the natural population and the second value for the hatchery
calc_Ford_fitness <- function(pNOB, pHOS, pbar_prev, theta, omega2, fitness_variance, A, fitness_floor, heritability) {

  pbar <- numeric(2)

  pbar[1] <- local({
    qHOS <- 1 - pHOS

    aa <- (pbar_prev[1] * omega2 + theta[1] * fitness_variance)/(omega2+fitness_variance)
    bb <- qHOS * (pbar_prev[1] + (aa - pbar_prev[1]) * heritability)
    cc <- (pbar_prev[2] * omega2 + theta[1] * fitness_variance)/(omega2+fitness_variance)
    dd <- pHOS * (pbar_prev[2] + (cc - pbar_prev[2]) * heritability)

    bb + dd
  })

  pbar[2] <- local({
    qNOB <- 1 - pNOB

    aa <- (pbar_prev[2] * omega2 + theta[2] * fitness_variance)/(omega2+fitness_variance)
    bb <- qNOB * (pbar_prev[2] + (aa - pbar_prev[2]) * heritability)
    cc <- (pbar_prev[1] * omega2 + theta[2] * fitness_variance)/(omega2+fitness_variance)
    dd <- pNOB * (pbar_prev[1] + (cc - pbar_prev[1]) * heritability)

    bb + dd
  })

  fitness <- local({
    x <- exp(-0.5 * (pbar[1] - theta[1])^2/(omega2+fitness_variance))
    max(fitness_floor, x)
  })

  return(list(pbar = pbar, fitness = fitness))
}

## INCOMPLETE
calc_Busack_fitness <- function(pNOS, Zpop_prev, theta) {

  # Zpop1
  pNOS <- NOS/(NOS + RRS_HOS * HOS_total)

  HOS_seg <- SAR_ratio * strays_total
  HOS_int <- HOS_total - HOS_seg

  pHOS_seg <- RRS_HOS * HOS_seg/(NOS + RRS_HOS * HOS_total)
  pHOS_int <- RRS_HOS * HOS_int/(NOS + RRS_HOS * HOS_total)

  Zpop[1, g] <- A * pNOS * (Zpop[1, g-1] - theta[1]) +
    pHOS_int * (Zpop[2, g-1] - theta[1]) + pHOS_seg * (Zpop[3, g-1] - theta[1]) + theta[1]

  pNOB1 <- brood_NOB/(brood_NOB + brood_HOB + brood_import)
  pHOBint1 <- brood_HOB/(brood_NOB + brood_HOB + brood_import)
  pHOBseg1 <- brood_import/(brood_NOB + brood_HOB + brood_import)

  Zpop[2, g] <- A * pNOB1 * (Zpop[1, g-1] - theta[2]) +
    pHOBint1 * (Zpop[2, g-1] - theta[2]) + pHOBseg1 * (Zpop[3, g-1] - theta[2]) + theta[2]

  pNOB2 <- 0 #BV60/(BV60+BW60+BX60)
  pHOBint2 <- max(1, export) #export #
  pHOBseg2 <- 0 #brood_import/(brood_NOB + brood_HOB + brood_import)

  #$I$49*(CE59*(CH59-$F$42)+CF59*(CI59-$F$42)+CG59*(CJ59-$F$42))+$F$42
  Zpop[3, g] <- A * pNOB2 * (Zpop[1, g-1] - theta[3]) +
    pHOBint2 * (Zpop[2, g-1] - theta[3]) + pHOBseg2 * (Zpop[3, g-1] - theta[3]) + theta[3]

  fitness <- abs(Zpop[1, g] - mean(theta[-1]))/abs(theta[1] - mean(theta[-1]))

  list(Zpop = Zpop, fitness = fitness)
}


catch_fn <- function(return_size, u, surv_passage = c(1, 1, 1, 1), nfishery = 4) {

  catch <- vector(length = nfishery)
  N <- vector(length = nfishery + 1)
  N[1] <- return_size

  for(i in 1:nfishery) {
    catch[i] <- N[i] * surv_passage[i] * u[i]
    N[i+1] <- N[i] * surv_passage[i] * (1 - u[i])
  }

  return(list(catch = catch, N = N))
}

.AHA_SRR <- function(Nnum, Nden = Nnum, p, capacity) {
  denom <- 1 + Nden * p /capacity
  Nnum * p/denom
}

.AHA_SRRpars <- function(prod, capacity, f, p_female) {
  alpha <- prod/f/p_female
  beta <- alpha/capacity
  c("alpha" = alpha, "beta" = beta)
}


