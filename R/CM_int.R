
#' Estimation function for conditioning model
#'
#' Population dynamics model of an age structured salmon population. Used with RTMB to estimate historical reconstruction from data.
#'
#' @param p List of parameter variables. See [fit_CM()].
#' @param d List of data variables. See [fit_CM()].
#' @keywords internal
#' @return Numeric, objective function value (log-posterior)
#' @author Q. Huynh with Stan code provided by J. Korman and C. Walters
#' @export
CM_int <- function(p, d) {

  # Dispatch method for AD variables ----
  is_ad <- any(sapply(p, inherits, "advector"))
  if (is_ad) {
    `[<-` <- RTMB::ADoverload("[<-")
  }

  # Transformed data ----
  tiny <- 1e-6
  moaddcwt <- -log(d$hatchsurv)
  spawnhist <- d$obsescape[1]
  logobsesc <- log(d$obsescape)

  # Vulnerability (= 0 for age 1, = 1 for oldest age)
  if (is_ad) {
    vulPT <- vulT <- advector(rep(0, d$Nages))
  } else {
    vulPT <- vulT <- numeric(d$Nages)
  }
  vulPT[seq(2, d$Nages-1)] <- plogis(p$logit_vulPT)
  vulPT[d$Nages] <- 1
  vulT[seq(2, d$Nages-1)] <- plogis(p$logit_vulT)
  vulT[d$Nages] <- 1

  lo <- numeric(d$Nages)
  if (is_ad) {
    lhist <- advector(rep(0, d$Nages))
  } else {
    lhist <- numeric(d$Nages)
  }

  lo[1] <- 1
  lhist[1] <- 1
  for (a in 2:d$Nages) {
    lo[a] <- lo[a-1] * exp(-d$mobase[a-1]) * (1 - d$bmatt[a-1]) # unfished juvenile survival
    lhist[a] <- lhist[a-1] * exp(-d$mobase[a-1] - vulPT[a] * d$finitPT) * (1 - d$bmatt[a-1]) # historical juvenile survival
  }

  epro <- d$ssum * sum(lo * d$fec * d$bmatt)                     # unfished egg production per smolt (recruit, pr)
  spro <- d$ssum * sum(lo * d$bmatt)                             # unfished spawners per smolt
  sprhist <- d$ssum * sum(lhist * exp(-vulPT * d$finitPT) * d$bmatt * exp(vulT * d$finitT))   # historcial spawners per recruit (pr)

  memax <- -log(1.0/epro) # unfished M from egg to smolt
  rhist <- spawnhist/sprhist # historical recruitment

  # Transformed parameters ----
  so <- exp(p$log_so)
  ro <- so/spro                 # unfished recruitment ro
  eo <- ro * epro               # unfished total egg production
  mden <- p$log_cr/eo           # Ricker b parameter for egg-smolt relationship
  memin <- memax - p$log_cr     # minimum egg-smolt M at low population density

  alpha <- exp(-memin)
  beta <- mden

  mo <- matrix(0, d$Ldyr, d$Nages-1) # annual ocean M by age

  # predicted survival rate from fishing
  # predicted catch at age vector
  # predicted spawners
  # predicted cwt catches for year
  # predicted cwt spawners for the year
  survPT <- ccwtPT <- survT <- ccwtT <- ecwt <- matrix(NA_real_, d$Ldyr, d$Nages)
  cyearPT <- cyearT <- recr <- syear <- escyear <- array(NA_real_, c(d$Ldyr, d$Nages, 2))

  spawners <- catchPT <- catchT <- egg <- megg <- logpredesc <- moplot <- FPT <- FT <- numeric(d$Ldyr)
  matt <- matrix(0, d$Ldyr, d$Nages)
  matt[, d$Nages] <- 1

  N <- array(0, c(d$Ldyr+1, d$Nages, 2)) # Array slice 1 = natural origin fish, 2 = hatchery fish
  Ncwt <- matrix(0, d$Ldyr+1, d$Nages)

  # Initialize mean phenotype, fitness, etc.
  brood <- array(NA_real_, c(d$Ldyr, d$Nages, 2))
  pNOB <- pHOSeff <- pHOScensus <- numeric(d$Ldyr)

  if (d$fitness) {
    zbar <- matrix(0, d$Ldyr, 2)
    fitness <- matrix(1, d$Ldyr, 2)
    omega <- sqrt(d$fitness_variance) * d$selection_strength
    omega2 <- omega * omega
  }

  # Initialize N ----
  N[1, , 1] <- rhist * lhist               # initial numbers at age year 1
  N[1, 1, 2] <- d$hatchsurv * d$hatchrelease[1] # initial age 1 numbers for hatchery release in year 1
  if (d$lht==2) {  #in case spring run type where age of ocean entry=2, not 1
    N[2, 1, ] <- N[1, 1, ]
  }

  # Initialize F
  FPT[] <- exp(p$log_FbasePT) * d$RelRegFPT * exp(p$log_fanomalyPT)
  FT[] <- exp(p$log_FbaseT) * d$RelRegFT * exp(p$log_fanomalyT)

  # Loop over years ----
  for (t in 1:d$Ldyr) {
    # First year ocean M
    if (length(d$covariate1)) {
      mo[t, 1] <- d$mobase[1] + p$moadd + p$wto[t] + sum(d$covariate1[t, ] * p$b1)
    } else {
      mo[t, 1] <- d$mobase[1] + p$moadd + p$wto[t]
    }

    # M's for older ages
    if (length(d$covariate)) {
      mo[t, 2:(d$Nages-1)] <- d$mobase[2:(d$Nages-1)] + sum(d$covariate[t, ] * p$b)
    } else {
      mo[t, 2:(d$Nages-1)] <- d$mobase[2:(d$Nages-1)]
    }

    Ncwt[t, 1] <- d$cwtrelease[t] * d$hatchsurv

    matt[t, 2:(d$Nages-1)] <- plogis(p$logit_matt[t, ])

    # set age-specific survival rates through fishing
    survPT[t, ] <- exp(-vulPT * FPT[t])
    survT[t, ] <- exp(-vulT * FT[t])

    # predict preterminal catch at age for the year
    cyearPT[t, , ] <- N[t, , ] * (1 - survPT[t, ])
    ccwtPT[t, ] <- Ncwt[t, ] * (1 - survPT[t, ])

    # predict recruitment at age for the year
    recr[t, , ] <- N[t, , ] * survPT[t, ] * matt[t, ]

    # predict terminal catch at age for the year
    cyearT[t, , ] <- recr[t, , ] * (1 - survT[t, ])
    ccwtT[t, ] <- Ncwt[t, ] * survPT[t, ] * matt[t, ] * (1 - survT[t, ])

    # predict escapement at age for the year
    escyear[t, , ] <- d$ssum * recr[t, , ] * survT[t, ]
    ecwt[t, ] <- d$ssum * Ncwt[t, ] * survPT[t, ] * matt[t, ] * survT[t, ]

    # predict spawners at age for the year
    syear[t, , ] <- d$propwildspawn[t] * escyear[t, , ] * d$s_enroute

    # predict egg production for the year
    egg[t] <- sum(d$fec * (syear[t, , 1] + d$gamma * syear[t, , 2]))

    # Assume broodtake is equal between NOB and HOB
    # predict pHOS, pNOB
    brood[t, , ] <- (1 - d$propwildspawn[t]) * escyear[t, , ]

    pNOB[t] <- sum(d$fec * brood[t, , 1])/sum(d$fec * (brood[t, , 1] + brood[t, , 2]))
    pHOSeff[t] <- sum(d$fec * d$gamma * syear[t, , 2])/
      sum(d$fec * (syear[t, , 1] + d$gamma * syear[t, , 2]))
    pHOScensus[t] <- sum(d$fec * syear[t, , 2])/sum(d$fec * (syear[t, , 1] + syear[t, , 2]))

    if (d$fitness) {
      # Retrieve zbar and fitness by brood year
      zbar_brood <- matrix(0, d$Nages, 2)
      for (a in 1:d$Nages) {
        tt <- t - a
        if (tt > 0) {
          zbar_brood[a, ] <- zbar[tt, ]
          if (a < d$Nages) {
            surv_fitness <- exp(-mo[t, a]) * fitness[tt, 1]^d$rel_loss[3]
            mo[t, a] <- -log(surv_fitness)
          }
        } else {
          zbar_brood[a, ] <- d$zbar_start
        }
      }
      zbar[t, ] <- calc_zbar(
        syear[t, , 1], d$gamma * syear[t, , 2], brood[t, , 1], brood[t, , 2],
        d$fec, d$fec, zbar_brood,
        omega2, d$theta, d$fitness_variance, d$heritability
      )
      fitness[t, ] <- exp(-0.5 * (zbar[t, ] - d$theta)^2/(omega2 + d$fitness_variance)) # fitness_floor not used here!
    }

    # survive fish over the year, removing maturing fish that will spawn
    N[t+1, 2:d$Nages, ] <- N[t, 2:d$Nages - 1, ] * survPT[t, 2:d$Nages - 1] * exp(-mo[t, 2:d$Nages - 1]) * (1 - matt[t, 2:d$Nages - 1])
    Ncwt[t+1, 2:d$Nages] <- Ncwt[t, 2:d$Nages - 1] * survPT[t, 2:d$Nages - 1] * exp(-mo[t, 2:d$Nages - 1]) * (1 - matt[t, 2:d$Nages - 1])

    if (d$fitness) {
      egg_fitness <- egg[t] * fitness[t, 1]^d$rel_loss[1]
      alpha_fitness <- alpha * fitness[t, 1]^d$rel_loss[2]
      beta_fitness <- beta/fitness[t, 1]^d$rel_loss[2]

      megg[t] <- memin + mden * egg_fitness + p$wt[t]
      N[t + d$lht, 1, 1] <- alpha_fitness * egg_fitness * exp(-beta_fitness * egg_fitness) * exp(-p$wt[t])
    } else {
      megg[t] <- memin + mden * egg[t] + p$wt[t]
      N[t + d$lht, 1, 1] <- alpha * egg[t] * exp(-beta * egg[t]) * exp(-p$wt[t])
    }
    N[t + d$lht, 1, 2] <- d$hatchsurv * d$hatchrelease[t+1]

    # total spawners and escapement
    spawners[t] <- sum(syear[t, , ])
    logpredesc[t] <- log(sum(escyear[t, , ]) + tiny)

    catchPT[t] <- sum(cyearPT[t, , ]) # add up total catch for the year
    catchT[t] <- sum(cyearT[t, , ]) # add up total catch for the year
    moplot[t] <- mo[t, 1] + moaddcwt # for plotting first year ocean M, including hatchery survival
  }

  # Convert predicted cwt catches and escapement from calendar year to brood year for likelihoods
  # data are in brood year format
  cbroodPT <- cbroodT <- ebrood <- matrix(tiny, d$Ldyr, d$Nages)
  for (t in 1:(d$Ldyr)) {
    for (a in 1:d$Nages) {
      if (t+a-1 <= d$Ldyr) {
        cbroodPT[t, a] <- d$cwtExp * ccwtPT[t+a-1, a] + tiny
        cbroodT[t, a] <- d$cwtExp * ccwtT[t+a-1, a] + tiny
        ebrood[t, a] <- d$cwtExp * ecwt[t+a-1, a] + tiny
      }
    }
  }

  # Log prior for parameters
  logprior_so <- dnorm(p$log_so, d$so_mu, d$so_sd, log = TRUE) # prior on so in log space as it is poorly determined from data
  logprior_wt <- dnorm(p$wt, 0, p$wt_sd, log = TRUE)
  logprior_wto <- dnorm(p$wto, 0, p$wto_sd, log = TRUE)

  if (sum(d$cwtcatPT)) {
    logprior_fanomPT <- dnorm(p$log_fanomalyPT, 0, p$fanomalyPT_sd, log = TRUE)
    logprior_fanomPT_sd <- dgamma(p$fanomalyPT_sd, 2, scale = 0.2, log = TRUE)
  } else if (is_ad) {
    logprior_fanomPT <- logprior_fanomPT_sd <- advector(0)
  } else {
    logprior_fanomPT <- logprior_fanomPT_sd <- 0
  }
  if (sum(d$cwtcatT)) {
    logprior_fanomT <- dnorm(p$log_fanomalyT, 0, p$fanomalyT_sd, log = TRUE)
    logprior_fanomT_sd <- dgamma(p$fanomalyT_sd, 2, scale = 0.2, log = TRUE)
  } else if (is_ad) {
    logprior_fanomT <- logprior_fanomT_sd <- advector(0)
  } else {
    logprior_fanomT <- logprior_fanomT_sd <- 0
  }

  # Log prior for vulnerability
  logit_bvulPT <- qlogis(squeeze(d$bvulPT[2:(d$Nages-1)]))
  logit_bvulT <- qlogis(squeeze(d$bvulT[2:(d$Nages-1)]))

  if (sum(d$bvulPT)) {
    logprior_vulPT <- dnorm(p$logit_vulPT, logit_bvulPT, 1.75, log = TRUE)
  } else if (is_ad) {
    logprior_vulPT <- advector(0)
  } else {
    logprior_vulPT <- 0
  }

  if (sum(d$bvulT)) {
    logprior_vulT <- dnorm(p$logit_vulT, logit_bvulT, 1.75, log = TRUE)
  } else if (is_ad) {
    logprior_vulT <- advector(0)
  } else {
    logprior_vulT <- 0
  }

  #convert base maturity rates to logit for calculation of time-varying rates
  #no time variation for first and last age
  logit_bmatt <- qlogis(d$bmatt[2:(d$Nages-1)])
  logprior_matt <- sapply(2:(d$Nages - 1), function(a) {
    dnorm(p$logit_matt[, a-1], logit_bmatt[a-1], p$sd_matt[a-1], log = TRUE)
  })

  logprior <- sum(
    logprior_so, logprior_wt, logprior_wto,
    logprior_fanomPT, logprior_fanomT, logprior_matt,
    logprior_vulPT, logprior_vulT
  )

  # Log prior for hyperparameters ----
  # Gamma mean=shape/rate, variance=shape/rate^2 (scale = 1/rate)
  # 2, 5 has mean of 0.4 and cv of 0.71 (95% CI = 0.03-1.0)
  logprior <- logprior +
    sum(
      dgamma(p$sd_matt, 2, scale = 0.2, log = TRUE),
      dgamma(p$wt_sd, 2, scale = 0.2, log = TRUE),
      dgamma(p$wto_sd, 2, scale = 0.2, log = TRUE),
      dgamma(p$lnE_sd, 2, scale = 0.2, log = TRUE)
    ) +
    logprior_fanomPT_sd + logprior_fanomT_sd

  # Log likelihood
  loglike_esc <- dnorm(logobsesc, logpredesc, p$lnE_sd, log = TRUE)
  loglike_cwtesc <- dpois(d$cwtesc, ebrood, log = TRUE)

  if (sum(d$cwtcatPT)) {
    loglike_cwtcatPT <- dpois(d$cwtcatPT, cbroodPT, log = TRUE)
  } else if (is_ad) {
    loglike_cwtcatPT <- advector(0)
  } else {
    loglike_cwtcatPT <- 0
  }
  if (sum(d$cwtcatT)) {
    loglike_cwtcatT <- dpois(d$cwtcatT, cbroodT, log = TRUE)
  } else if (is_ad) {
    loglike_cwtcatT <- advector(0)
  } else {
    loglike_cwtcatT <- 0
  }

  loglike <- sum(loglike_esc, loglike_cwtcatPT, loglike_cwtcatT, loglike_cwtesc)

  # Objective function
  fn <- -1 * (logprior + loglike)

  # Report variables
  REPORT(so)
  REPORT(eo)
  REPORT(ro)
  REPORT(mden)
  REPORT(memax)
  REPORT(memin)
  REPORT(alpha)
  REPORT(beta)

  REPORT(survPT)
  REPORT(survT)
  REPORT(cyearPT)
  REPORT(cyearT)
  REPORT(recr)
  REPORT(syear)
  REPORT(escyear)
  REPORT(ccwtPT)
  REPORT(ccwtT)
  REPORT(ecwt)

  REPORT(spawners)
  REPORT(logpredesc)
  REPORT(catchPT)
  REPORT(catchT)
  REPORT(egg)
  REPORT(megg)
  REPORT(mo)
  REPORT(moplot)
  REPORT(FPT)
  REPORT(FT)
  REPORT(vulPT)
  REPORT(vulT)

  REPORT(matt)
  REPORT(N)
  REPORT(Ncwt)

  REPORT(cbroodPT)
  REPORT(cbroodT)
  REPORT(ebrood)

  REPORT(epro)

  REPORT(brood)
  REPORT(pNOB)
  REPORT(pHOSeff)
  REPORT(pHOScensus)

  if (d$fitness) {
    REPORT(fitness)
    REPORT(zbar)
  }

  REPORT(loglike_esc)
  REPORT(loglike_cwtcatPT)
  REPORT(loglike_cwtcatT)
  REPORT(loglike_cwtesc)

  REPORT(loglike)
  REPORT(logprior)
  REPORT(fn)

  return(fn)
}


# Make list of starting values
make_CMpars <- function(p, d) {
  if (is.null(p$log_cr)) p$log_cr <- 3
  if (is.null(p$log_so)) p$log_so <- log(3 * max(d$obsescape))
  if (is.null(p$moadd)) p$moadd <- 0
  if (is.null(p$wt)) p$wt <- rep(0, d$Ldyr)
  if (is.null(p$wto)) p$wto <- rep(0, d$Ldyr)
  if (is.null(p$log_fanomalyPT)) p$log_fanomalyPT <- rep(0, d$Ldyr)
  if (is.null(p$log_fanomalyT)) p$log_fanomalyT <- rep(0, d$Ldyr)
  if (is.null(p$lnE_sd)) p$lnE_sd <- 0.1

  if (is.null(p$log_FbasePT)) p$log_FbasePT <- log(0.1)
  if (is.null(p$log_FbaseT)) p$log_FbaseT <- log(0.1)
  if (is.null(p$logit_vulPT)) p$logit_vulPT <- qlogis(squeeze(d$bvulPT[-c(1, d$Nages)]))
  if (is.null(p$logit_vulT)) p$logit_vulT <- qlogis(squeeze(d$bvulT[-c(1, d$Nages)]))

  if (is.null(p$logit_matt)) p$logit_matt <- matrix(qlogis(d$bmatt[-c(1, d$Nages)]), d$Ldyr, d$Nages - 2, byrow = TRUE)
  if (is.null(p$sd_matt)) p$sd_matt <- rep(0.5, d$Nages - 2)

  if (is.null(p$wt_sd)) p$wt_sd <- 1
  if (is.null(p$wto_sd)) p$wto_sd <- 1
  if (is.null(p$fanomalyPT_sd)) p$fanomalyPT_sd <- 1
  if (is.null(p$fanomalyT_sd)) p$fanomalyT_sd <- 1

  if (is.null(p[["b1"]])) {
    if (length(d[["covariate1"]])) {
      p[["b1"]] <- rep(0, ncol(d[["covariate1"]]))
    } else {
      p[["b1"]] <- 0
    }
  }

  if (is.null(p[["b"]])) {
    if (length(d[["covariate"]])) {
      p[["b"]] <- rep(0, ncol(d[["covariate"]]))
    } else {
      p[["b"]] <- 0
    }
  }

  return(p)
}

check_data <- function(data) {

  if (is.null(data$Nages)) stop("data$Nages not found")
  if (is.null(data$Ldyr)) stop("data$Ldyr not found")
  if (is.null(data$lht)) data$lht <- 1

  if (is.null(data$hatchsurv)) stop("data$hatchsurv should be between 0-1")
  if (is.null(data$gamma)) data$gamma <- 1
  if (is.null(data$ssum)) stop("data$ssum should be between 0-1")

  if (is.null(data$finitPT)) data$finitPT <- 0
  if (is.null(data$finitT)) data$finitT <- 0

  if (is.null(data$bmatt) || length(data$bmatt) != data$Nages) {
    stop("data$bmatt should be a vector length Nages")
  }

  if (is.null(data$fec) || length(data$fec) != data$Nages) {
    stop("data$fec should be a vector length Nages")
  }

  if (is.null(data$cwtcatPT) || !is.matrix(data$cwtcatPT)) {
    stop("data$cwtcatPT should be a matrix: Ldyr (by brood year) x Nages")
  }
  if (is.null(data$cwtcatT) || !is.matrix(data$cwtcatT)) {
    stop("data$cwtcatT should be a matrix: Ldyr (by brood year) x Nages")
  }

  FPT <- sum(data$cwtcatPT)
  FT <- sum(data$cwtcatT)
  if (!FPT && !FT) stop("No CWT preterminal or terminal catch found")

  if (FPT) {
    if (is.null(data$bvulPT) || length(data$bvulPT) != data$Nages) {
      stop("data$bvulPT should be a vector length Nages")
    }
    if (is.null(data$RelRegFPT) || length(data$RelRegFPT) != data$Ldyr) {
      stop("data$RelRegFPT should be a vector length Ldyr")
    }
  } else {
    data$bvulPT <- rep(0, data$Nages)
    data$RelRegFPT <- rep(0, data$Ldyr)
  }

  if (FT) {
    if (is.null(data$bvulT) || length(data$bvulT) != data$Nages) {
      stop("data$bvulT should be a vector length Nages")
    }
    if (is.null(data$RelRegFT) || length(data$RelRegFT) != data$Ldyr) {
      stop("data$RelRegFT should be a vector length Ldyr")
    }
  } else {
    data$bvulT <- rep(0, data$Nages)
    data$RelRegFT <- rep(0, data$Ldyr)
  }

  if (is.null(data$mobase) || length(data$mobase) != data$Nages) {
    stop("data$mobase should be a vector length Nages")
  }

  if (is.null(data$cwtrelease) || length(data$cwtrelease) != data$Ldyr) {
    stop("data$cwtrelease should be a vector length Ldyr")
  }

  if (is.null(data$obsescape) || length(data$obsescape) != data$Ldyr) {
    stop("data$obsescape should be a vector length Ldyr")
  }

  if (is.null(data$propwildspawn) || length(data$propwildspawn) != data$Ldyr) {
    stop("data$propwildspawn should be a vector length Ldyr")
  }

  if (is.null(data$hatchrelease) || length(data$hatchrelease) != data$Ldyr + 1) {
    stop("data$hatchrelease should be a vector length Ldyr+1")
  }

  if (is.null(data$cwtExp)) data$cwtExp <- 1

  if (!is.null(data$covariate1)) {
    if (!is.matrix(data$covariate1)) stop("data$covariate1 should be a matrix. See help('fit-CM') for dimensions")
  }

  if (!is.null(data$covariate)) {
    if (!is.matrix(data$covariate)) stop("data$covariate should be a matrix. See help('fit-CM') for dimensions")
  }

  if (is.null(data$s_enroute)) data$s_enroute <- 1

  if (is.null(data$so_mu)) data$so_mu <- log(3 * max(data$obsescape))
  if (is.null(data$so_sd)) data$so_sd <- 0.5

  if (is.null(data$fitness)) data$fitness <- FALSE

  if (data$fitness) {
    if (is.null(data$theta)) data$theta <- c(100, 80)
    if (is.null(data$rel_loss)) data$rel_loss <- c(0.3, 0.2, 0.5)
    if (is.null(data$zbar_start)) data$zbar_start <- c(100, 100)
    if (is.null(data$fitness_variance)) data$fitness_variance <- 10
    if (is.null(data$selection_strength)) data$selection_strength <- 3
    if (is.null(data$heritability)) data$heritability <- 0.5
    if (is.null(data$fitness_floor)) data$fitness_floor <- 0.5
  }

  return(data)
}


make_map <- function(map = list(), p, d) {

  FPT <- sum(d$cwtcatPT)
  FT <- sum(d$cwtcatT)

  if (!FPT) {
    map$FbasePT <- factor(NA)
    map$fanomalyPT <- factor(rep(NA, d$Ldyr))
    map$fanomalyPT_sd <- factor(NA)
  }

  if (!FT) {
    map$FbaseT <- factor(NA)
    map$fanomalyT <- factor(rep(NA, d$Ldyr))
    map$fanomalyT_sd <- factor(NA)
  }

  if (!length(d[["covariate1"]])) map[["b1"]] <- factor(NA)
  if (!length(d[["covariate"]])) map[["b"]] <- factor(NA)
  return(map)
}

make_bounds <- function(par_names, data, lower = list(), upper = list()) {

  .lower <- structure(rep(-Inf, length(par_names)), names = par_names)
  .upper <- structure(rep(Inf, length(par_names)), names = par_names)

  # Add important defaults first
  if ("log_cr" %in% names(.lower)) .lower["log_cr"] <- 1e-8
  if ("moadd" %in% names(.lower)) .lower["moadd"] <- 1e-8
  #if ("FbasePT" %in% names(.lower)) .lower["FbasePT"] <- 1e-8
  #if ("FbaseT" %in% names(.lower)) .lower["FbaseT"] <- 1e-8

  if ("lnE_sd" %in% names(.lower)) .lower["lnE_sd"] <- 1e-8
  if ("wt_sd" %in% names(.lower)) .lower["wt_sd"] <- 1e-8
  if ("wto_sd" %in% names(.lower)) .lower["wto_sd"] <- 1e-8
  if ("fanomalyPT_sd" %in% names(.lower)) .lower["fanomalyPT_sd"] <- 1e-8
  if ("fanomalyT_sd" %in% names(.lower)) .lower["fanomalyT_sd"] <- 1e-8

  if ("sd_matt" %in% names(.lower)) .lower[names(.lower) == "sd_matt"] <- 1e-8

  #if ("wt" %in% names(.lower)) {
  #  .lower[names(.lower) == "wt"] <- -3
  #  .upper[names(.upper) == "wt"] <- 3
  #}
  #if ("wto" %in% names(.lower)) {
  #  .lower[names(.lower) == "wto"] <- -3
  #  .upper[names(.upper) == "wto"] <- 3
  #}
  #if ("fanomalyPT" %in% names(.lower)) {
  #  .lower[names(.lower) == "fanomalyPT"] <- -3
  #  .upper[names(.upper) == "fanomalyPT"] <- 3
  #}
  #if ("fanomalyT" %in% names(.lower)) {
  #  .lower[names(.lower) == "fanomalyT"] <- -3
  #  .upper[names(.upper) == "fanomalyT"] <- 3
  #}

  # Override with user-defined bounds
  for (i in unique(par_names)) {
    if (!is.null(lower[[i]])) .lower[names(.lower) == i] <- lower[[i]]
    if (!is.null(upper[[i]])) .upper[names(.upper) == i] <- upper[[i]]
  }

  list(lower = .lower, upper = .upper)
}

get_report <- function(stanfit, sims) {
  if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan package is needed.")

  pars <- rstan::extract(stanfit)
  if (missing(sims)) sims <- seq_len(length(pars[["lp__"]]))

  pars_samp <- lapply(pars[names(pars) != "lp__"], function(x) {
    if (is.matrix(x)) x[sims, , drop = FALSE] else x[sims]
  })

  fit <- stanfit@.MISC$CMfit
  if (is.null(fit)) stop("CM fitted object not found in stanfit@.MISC$CMfit")
  report <- lapply(1:length(sims), function(i) {
    par_x <- lapply(pars_samp, function(x) {
      if (is.matrix(x)) x[i, ] else x[i]
    })
    fit$obj$report(unlist(par_x))
  })

  return(report)
}

get_CMdata <- function(fit) {
  func <- attr(fit$obj$env$data, "func")
  get("data", envir = environment(func), inherits = FALSE)
}


#' @importFrom abind abind
expand_array <- function(x, proyears) {
  xlast <- x[, , dim(x)[3]]
  xproj <- array(xlast, c(dim(x)[1:2], proyears))
  abind::abind(x, xproj, along = 3) %>%
    structure(dimnames = NULL)
}

squeeze <- function(x) (1 - .Machine$double.eps) * (x - 0.5) + 0.5
