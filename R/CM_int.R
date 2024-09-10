
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

  lo <- lhist <- numeric(d$Nages)
  lo[1] <- lhist[1] <- 1
  for (a in 2:d$Nages) {
    lo[a] <- lo[a-1] * exp(-d$mobase[a-1]) * (1 - d$bmatt[a-1]) # unfished survivorship recursion
    lhist[a] <- lhist[a-1] * exp(-d$mobase[a-1] - d$vul[a-1] * d$fhist) * (1 - d$bmatt[a-1]) # historical survivorship
  }

  epro <- d$ssum * sum(lo * d$fec * d$bmatt)                     # unfished egg production per smolt (recruit, pr)
  spro <- d$ssum * sum(lo * d$bmatt)                             # unfished spawners per smolt
  sprhist <- d$ssum * sum(lhist * exp(-d$vul * d$fhist) * d$bmatt)   # historcial spawners per recruit (pr)

  memax <- -log(1.0/epro) # unfished M from egg to smolt
  rhist <- spawnhist/sprhist # historical recruitment

  # Transformed parameters ----
  so <- exp(p$log_so)
  ro <- so/spro            # unfished recruitment ro
  eo <- ro * epro          # unfished total egg production
  mden <- p$cr/eo          # Ricker b parameter for egg-smolt relationship
  memin <- memax - p$cr    # minimum egg-smolt M at low population density

  alpha <- exp(-memin)
  beta <- mden

  mo <- matrix(0, d$Ldyr, d$Nages-1) # annual ocean M by age

  # predicted survival rate from fishing
  # predicted catch at age vector
  # predicted spawners
  # predicted cwt catches for year
  # predicted cwt spawners for the year
  sfish <- cyear <- syear <- ccwt <- scwt <- matrix(NA_real_, d$Ldyr, d$Nages)
  cyear <- syear <- array(NA_real_, c(d$Ldyr, d$Nages, 2))
  spawners <- tcatch <- egg <- megg <- logpredesc <- moplot <- ft <- numeric(d$Ldyr)
  matt <- matrix(0, d$Ldyr, d$Nages)
  matt[, d$Nages] <- 1

  N <- array(0, c(d$Ldyr+1, d$Nages, 2)) # Array slice 1 = natural origin fish, 2 = hatchery fish
  Ncwt <- matrix(0, d$Ldyr+1, d$Nages)

  # Initialize N ----
  N[1, , 1] <- rhist * lhist               # initial numbers at age year 1
  N[1, 1, 2] <- d$hatchsurv * d$hatchrelease[1] # initial age 1 numbers for hatchery release in year 1
  if (d$lht==2) {  #in case spring run type where age of ocean entry=2, not 1
    N[2, 1, ] <- N[1, 1, ]
  }

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

    ft[t] <- p$Fbase * d$RelRegF[t] + p$fanomaly[t]

    matt[t, 2:(d$Nages-1)] <- plogis(p$logit_matt[t, ])
    sfish[t, ] <- exp(-d$vul * ft[t])  # set age-specific survival rates through fishing

    cyear[t, , ] <- N[t, , ] * (1 - sfish[t, ]) # predict catch at age for the year
    syear[t, , ] <- d$ssum * N[t, , ] * sfish[t, ] * matt[t, ] # predict spawners at age for the year

    ccwt[t, ] <- Ncwt[t, ] * (1 - sfish[t, ])
    scwt[t, ] <- d$ssum * Ncwt[t, ] * sfish[t, ] * matt[t, ]

    sp_t <- syear[t, , 1] + d$gamma * syear[t, , 2]

    egg[t] <- sum(sp_t * d$fec * d$propwildspawn[t])

    # survive fish over the year, removing maturing fish that will spawn
    N[t+1, 2:d$Nages, ] <- N[t, 2:d$Nages - 1, ] * sfish[t, 2:d$Nages - 1] * exp(-mo[t, 2:d$Nages - 1]) * (1 - matt[t, 2:d$Nages - 1])
    Ncwt[t+1, 2:d$Nages] <- Ncwt[t, 2:d$Nages - 1] * sfish[t, 2:d$Nages - 1] * exp(-mo[t, 2:d$Nages - 1]) * (1 - matt[t, 2:d$Nages - 1])

    megg[t] <- memin + mden * egg[t] + p$wt[t]
    N[t + d$lht, 1, 1] <- alpha * egg[t] * exp(-beta * egg[t]) * exp(-p$wt[t])
    N[t + d$lht, 1, 2] <- d$hatchsurv * d$hatchrelease[t+1]

    spawners[t] <- sum(syear[t, , ])
    logpredesc[t] <- log(spawners[t] + tiny)

    tcatch[t] <- sum(cyear[t, , ]) # add up total catch for the year
    moplot[t] <- mo[t, 1] + moaddcwt # for plotting first year ocean M, including hatchery survival
  }

  # Convert predicted cwt catches and escapement from calendar year to brood year for likelihoods
  # data are in brood year format
  cbrood <- ebrood <- matrix(tiny, d$Ldyr, d$Nages)
  for (t in 1:(d$Ldyr)) {
    for (a in 1:d$Nages) {
      if (t+a-1 <= d$Ldyr) {
        cbrood[t, a] <- ccwt[t+a-1, a] + tiny
        ebrood[t, a] <- scwt[t+a-1, a] + tiny
      }
    }
  }

  cbrood2 <- cbrood * d$cwtExp
  ebrood2 <- ebrood * d$cwtExp

  # Log prior for parameters
  logprior_so <- dnorm(p$log_so, d$so_mu, d$so_sd, log = TRUE) # prior on so in log space as it is poorly determined from data
  logprior_wt <- dnorm(p$wt, 0, p$wt_sd, log = TRUE)
  logprior_wto <- dnorm(p$wto, 0, p$wto_sd, log = TRUE)
  logprior_f <- dnorm(p$fanomaly, 0, p$fanomaly_sd, log = TRUE)

  #convert base maturity rates to logit for calculation of time-varying rates
  #no time variation for first and last age
  logit_bmatt <- qlogis(d$bmatt[2:(d$Nages-1)])
  logprior_matt <- sapply(2:(d$Nages - 1), function(a) {
    dnorm(p$logit_matt[, a-1], logit_bmatt[a-1], p$sd_matt[a-1], log = TRUE)
  })

  logprior <- sum(logprior_so, logprior_wt, logprior_wto, logprior_f, logprior_matt)

  # Log prior for hyperparameters ----
  # Gamma mean=shape/rate, variance=shape/rate^2 (scale = 1/rate)
  # 2, 5 has mean of 0.4 and cv of 0.71 (95% CI = 0.03-1.0)
  logprior <- logprior +
    sum(
      dgamma(p$sd_matt, 2, scale = 0.2, log = TRUE),
      dgamma(p$wt_sd, 2, scale = 0.2, log = TRUE),
      dgamma(p$wto_sd, 2, scale = 0.2, log = TRUE),
      dgamma(p$fanomaly_sd, 2, scale = 0.2, log = TRUE),
      dgamma(p$lnS_sd, 2, scale = 0.2, log = TRUE)
    )

  # Log likelihood
  loglike_esc <- dnorm(logobsesc, logpredesc, p$lnS_sd, log = TRUE) # likelihood for spawner predictions
  loglike_cwtcat <- dpois(d$cwtcat, cbrood2, log = TRUE)
  loglike_cwtesc <- dpois(d$cwtesc, ebrood2, log = TRUE)

  loglike <- sum(loglike_esc, loglike_cwtcat, loglike_cwtesc)

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

  REPORT(sfish)
  REPORT(cyear)
  REPORT(syear)
  REPORT(ccwt)
  REPORT(scwt)

  REPORT(spawners)
  REPORT(tcatch)
  REPORT(egg)
  REPORT(megg)
  REPORT(mo)
  REPORT(moplot)
  REPORT(ft)

  REPORT(matt)
  REPORT(N)
  REPORT(Ncwt)

  REPORT(cbrood)
  REPORT(ebrood)

  REPORT(epro)

  REPORT(loglike_esc)
  REPORT(loglike_cwtcat)
  REPORT(loglike_cwtesc)

  REPORT(loglike)
  REPORT(logprior)
  REPORT(fn)

  return(fn)
}


# Make list of starting values
make_CMpars <- function(p, d) {
  if (is.null(p$cr)) p$cr <- 3
  if (is.null(p$log_so)) p$log_so <- log(3 * max(d$obsescape))
  if (is.null(p$moadd)) p$moadd <- 0
  if (is.null(p$wt)) p$wt <- rep(0, d$Ldyr)
  if (is.null(p$wto)) p$wto <- rep(0, d$Ldyr)
  if (is.null(p$fanomaly)) p$fanomaly <- rep(0, d$Ldyr)
  if (is.null(p$lnS_sd)) p$lnS_sd <- 0.1

  if (is.null(p$Fbase)) p$Fbase <- 1
  if (is.null(p$logit_matt)) p$logit_matt <- matrix(0.1, d$Ldyr, d$Nages - 2)
  if (is.null(p$sd_matt)) p$sd_matt <- rep(0.5, d$Nages - 2)

  if (is.null(p$wt_sd)) p$wt_sd <- 1
  if (is.null(p$wto_sd)) p$wto_sd <- 1
  if (is.null(p$fanomaly_sd)) p$fanomaly_sd <- 1

  if (is.null(p$b1)) {
    if (length(d$covariate1)) {
      p$b1 <- rep(0, ncol(d$covariate1))
    } else {
      p$b1 <- 0
    }
  }
  if (is.null(p$b)) {
    if (length(d$covariate)) {
      p$b <- rep(0, ncol(d$covariate))
    } else {
      p$b <- 0
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

  if (is.null(data$fhist)) stop("data$fhist should be a numeric")

  if (is.null(data$bmatt) || length(data$bmatt) != data$Nages) {
    stop("data$bmatt should be a vector length Nages")
  }

  if (is.null(data$fec) || length(data$fec) != data$Nages) {
    stop("data$fec should be a vector length Nages")
  }

  if (is.null(data$vul) || length(data$vul) != data$Nages) {
    stop("data$vul should be a vector length Nages")
  }

  if (is.null(data$mobase) || length(data$mobase) != data$Nages) {
    stop("data$mobase should be a vector length Nages")
  }

  if (is.null(data$cwtrelease) || length(data$cwtrelease) != data$Ldyr) {
    stop("data$cwtrelease should be a vector length Ldyr")
  }

  if (is.null(data$cwtesc) || !is.matrix(data$cwtesc)) {
    stop("data$cwtesc should be a matrix: Ldyr (by brood year) x Nages")
  }

  if (is.null(data$cwtcat) || !is.matrix(data$cwtcat)) {
    stop("data$cwtcat should be a matrix: Ldyr (by brood year) x Nages")
  }

  if (is.null(data$RelRegF) || length(data$RelRegF) != data$Ldyr) {
    stop("data$RelRegF should be a vector length Ldyr")
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
    if (!is.matrix(data$covariate1)) stop("data$covariate1 should be a matrix")
  }

  if (!is.null(data$covariate1)) {
    if (!is.matrix(data$covariate1)) stop("data$covariate1 should be a matrix. See help('fit-CM') for dimensions")
  }

  if (!is.null(data$covariate)) {
    if (!is.matrix(data$covariate)) stop("data$covariate should be a matrix. See help('fit-CM') for dimensions")
  }


  if (is.null(data$so_mu)) data$so_mu <- log(3 * max(data$obsescape))
  if (is.null(data$so_sd)) data$so_sd <- 0.5

  return(data)
}

get_report <- function(stanfit, fit, sims) {
  if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan package is needed.")

  pars <- rstan::extract(stanfit)
  if (missing(sims)) sims <- seq_len(length(pars[["lp__"]]))

  pars_samp <- lapply(pars[names(pars) != "lp__"], function(x) {
    if (is.matrix(x)) x[sims, , drop = FALSE] else x[sims]
  })

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
