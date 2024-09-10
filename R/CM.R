


#' Fit conditioning model to historical data
#'
#' Bayesian stock reconstruction model of natural and hatchery origin fish,
#' fitted to coded wire tag data, observed escapement, and hatchery releases.
#' Estimates time-varying maturity rate as well as time-varying ocean survival as a linear model of covariates (separate covariates
#' assumed between age 1 vs. ages 2+).
#' Currently models only a preterminal fishery.
#' `fit_CM()` generates the RTMB model from data which can then be passed to `sample_CM()` to run the MCMC in Stan.
#'
#' @param data  A list containing data inputs. See details.
#' @param start An optional list containing parameter starting values. See details.
#' @param lower_b1 Lower bound of coefficients for linear covariates that predict natural mortality for age 1. See details.
#' @param upper_b1 Upper bound of coefficients for linear covariates that predict natural mortality for age 1. See details.
#' @param lower_b Lower bound of coefficients for linear covariates that predict natural mortality for age 2+ See details.
#' @param upper_b Upper bound of coefficients for linear covariates that predict natural mortality for age 2+. See details.
#' @param do_fit Logical, whether to do the fit and estimate the Hessian.
#' @param silent Logical, whether to silence output from RTMB to the console.
#' @param ... For `fit_CM`, arguments to [`RTMB::MakeADFun()`]. For `sample_CM`, arguments to `rstan::sampling()`
#' @import RTMB
#' @importFrom stats nlminb
#' @returns
#' - `fit_CM()` returns a list containing the RTMB model, nlminb output, standard errors, and parameter bounds
#' - `sample_CM()` returns a `stanfit` object containing the MCMC chains
#'
#' @section Data:
#' Data should passed through a named list with the following:
#' - `Nages` Integer, number of age classes in the model
#' - `Ldyr` Integer, number of years in the model
#' - `lht` Integer, life history type. Should be 1 for now
#' - `hatchsurv` Numeric, survival of hatchery releases into the smolt life stage. Density-independent.
#' - `ssum` Numeric, proportion of spawners that is female
#' - `fhist` Numeric, historical fishing mortality for calculating the spawners at age in the first year of the model
#' - `bmatt` Vector length `Nages`. Proportion maturity at age, base values for calculating the unfished replacement line
#' - `fec` Vector length `Nages`. Fecundity, egg production at age
#' - `vul` Vector length `Nages`. Vulnerability at age to the preterminal fishery
#' - `mobase`. Vector length `Nages`. Natural mortality at age, base values for calculating the unfished replacement line
#' - `cwtrelease` Vector length `Ldyr`, coded wire tag (CWT) releases
#' - `cwtesc` Matrix `[Ldyr, Nages]`. CWT escapement **by brood year and age**. Poisson likelhood.
#' - `cwtcat` Matrix `[Ldyr, Nages]`. CWT catch, **by brood year and age**. Poisson likelhood.
#' - `RelRegF` Vector `Ldyr`. Trend in relative regional fishing mortality. Fishing mortality is estimated by estimating a scaling
#' coefficient and annual deviations from this vector.
#' - `obsescape` Vector length `Ldyr`, total observed escapement (all ages and both hatchery/natural fish). Lognormal likelhood.
#' - `propwildspawn` Vector length `Ldyr`, proportion of the escapement that spawn (accounts for en-route mortality and broodtake)
#' - `hatchrelease` Vector length `Ldyr+1`, number of hatchery juvenile fish released
#' - `cwtExp` *Optional*. Numeric, coefficient for scaling down the CWT predictions to match the observations. For example, `cwtExp = 10` means that
#' the observed values is ten times greater than the predicted value. Default is 1.
#' - `covariate1` *Optional*. Matrix `Ldyr, ncov1` of linear covariates that predict natural mortality for age 1.
#' - `covariate` *Optional*. Matrix `Ldyr, ncov` of linear covariates that predict natural mortality for ages 2+.
#' - `so_mu` Numeric, the prior mean for unfished spawners in logspace
#' - `so_sd` Numeric, the prior standard deviation for unfished spawners in logspace
#' - `so_min` *Optional*. Numeric, lower bound for the estimate of unfished spawners.
#' @section start:
#' Starting values for parameters can be provided through a named list:
#'
#' - `cr` Numeric, compensation ratio. Default is 3.
#' - `log_so` Numeric, unfished spawners in logspace. Default is 2.
#' - `moadd` Numeric, additive term to base natural mortality rate for age 1 juveniles. Default is zero.
#' - `wt` Vector `Ldyr`. Annual deviates in natural mortality during the freshwater life stage (affects survival to smolt life stage).
#' Estimated with normal prior with mean zero and standard deviation `p$wt_sd`. Default is zero.
#' - `wto` Vector `Ldyr`. Annual deviates in natural mortality for age 2+ juveniles (marine life stage).
#' Estimated with normal prior with mean zero and standard deviation `p$wto_sd`. Default is zero.
#' - `Fbase` Numeric, scaling coefficient to estimate fishing mortality from `data$RelRegF`. Default is 1.
#' - `fanomaly` Vector `Ldyr`. Annual deviates from `Fbase * data$RelRegF` to estimate fishing mortality.
#' Estimated with normal prior with mean zero and standard deviation `p$fanomaly_sd`. Default is zero.
#' - `lnS_sd` Numeric, lognormal standard deviation of the observed escapement. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `wt_sd` Numeric, lognormal standard deviation of the age 1 (freshwater) natural mortality deviates. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `wto_sd` Numeric, lognormal standard deviation of the age 2+ (marine) natural mortality deviates. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `fanomaly_sd` Numeric, lognormal standard deviation of `fanomaly`. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `logit_matt` Matrix `Ldyr, Nages-2` maturity at age in logit space. Age-1 maturity is fixed to zero and one at age 1 and the maximum age, respectively. Default is 0.1.
#' - `sd_matt` Vector `Nages-2`. Logit standard deviation of maturity (`logit_matt`) by age class. Default is 0.5.
#' - `b1` Vector `ncov1` of coefficients for linear covariates that predict natural mortality for age 1. Default is zero.
#' - `b` Vector `ncov` of coefficients for linear covariates that predict natural mortality for ages 2+. Default is zero.
#'
#' @section Covariates on natural mortality:
#'
#' Natural mortality is modeled as the sum of a base value \eqn{M^\textrm{base}}, additional scaling factor for age 1 \eqn{M^\textrm{add}},
#' a linear system of covariates \eqn{X} and coefficients \eqn{b}:
#'
#' \deqn{
#' M_{y,a} =
#' \begin{cases}
#' M^\textrm{base}_a + M^\textrm{add} + \sum_j b^1_j X^1_{y,j} & \quad a = 1
#' M^\textrm{base}_a + \sum_j b_j X_{y,j} & \quad a = 2, \ldots, A
#' \end{cases}
#' }
#' @author Q. Huynh with Stan code provided by J. Korman and C. Walters
#' @seealso [CM2SOM()]
#' @export
fit_CM <- function(data, start = list(), lower_b1, upper_b1, lower_b, upper_b, do_fit = TRUE, silent = TRUE, ...) {

  data <- check_data(data)

  p <- make_CMpars(start, data)

  f <- function(p) salmonMSE::CM_int(p, d = data)

  map <- list()
  if (!length(data$covariate1)) map$b1 <- factor(NA)
  if (!length(data$covariate)) map$b <- factor(NA)

  obj <- RTMB::MakeADFun(func = f, parameters = p, silent = silent, ...)

  lower <- structure(rep(-Inf, length(obj$par)), names = names(obj$par))
  upper <- structure(rep(Inf, length(obj$par)), names = names(obj$par))

  lower["cr"] <- 1
  if (is.null(data$maxcr)) {
    upper["cr"] <- 10
  } else {
    upper["cr"] <- data$maxcr
  }

  if (!is.null(data$so_min)) lower["log_so"] <- data$so_min

  lower["moadd"] <- -2/3 * data$mobase[1]
  lower["Fbase"] <- 1e-8

  if (!missing(lower_b1)) lower[names(lower) == "b1"] <- lower_b1
  if (!missing(upper_b1)) upper[names(upper) == "b1"] <- upper_b1

  if (!missing(lower_b)) lower[names(lower) == "b"] <- lower_b
  if (!missing(upper_b)) upper[names(upper) == "b"] <- upper_b

  lower["lnS_sd"] <- 1e-8
  lower["wt_sd"] <- 1e-8
  lower["wto_sd"] <- 1e-8
  lower["fanomaly_sd"] <- 1e-8

  lower[names(lower) == "sd_matt"] <- 1e-8

  lower[names(lower) == "wt"] <- lower[names(lower) == "wto"] <- lower[names(lower) == "fanomaly"] <- -3
  upper[names(upper) == "wt"] <- upper[names(upper) == "wto"] <- upper[names(upper) == "fanomaly"] <- 3

  if (do_fit) {
    opt <- nlminb(
      obj$par, obj$fn, obj$gr,
      lower = lower, upper = upper, control = list(eval.max = 1e5, iter.max = 1e5)
    )
    SD <- RTMB::sdreport(obj)
  } else {
    opt <- SD <- NULL
  }

  res <- list(obj = obj, opt = opt, SD = SD, lower = lower, upper = upper)
  return(res)
}

#' @rdname fit_CM
#' @param fit List of output from `fit_CM()`
#' @export
sample_CM <- function(fit, ...) {
  if (!requireNamespace("tmbstan", quietly = TRUE)) stop("tmbstan package is needed.")
  tmbstan::tmbstan(fit$obj, lower = fit$lower, upper = fit$upper, ...)
}


#' Convert conditioning model to operating model
#'
#' Creates an operating model from MCMC samples and data inputs of the conditioning model.
#' Management actions for habitat, hatchery production, and harvest still need to be specified in the operating model.
#'
#' @param stanfit Output from [`sample_CM()`]
#' @param fit Output from [`fit_CM()`]
#' @param sims Optional, a vector of integers indicating the MCMC iterations to convert to operating model simulations. Otherwise,
#' use argument `nsim` in order to sample a subset of the MCMC.
#' @param nsim Integer, total number of simulations in the operating model. Only used if `sims` is missing.
#' @param seed Integer, seed for sampling the MCMC output. Only used if `sims` is missing.
#' @param proyears Integer, the number of projection years in the operating model
#' @return \linkS4class{SOM} object.
#' @export
CM2SOM <- function(stanfit, fit, sims, nsim = 2, seed = 1, proyears = 40) {

  if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan package is needed.")
  pars <- rstan::extract(stanfit)

  if (missing(sims)) {
    nsim_stan <- length(pars[["lp__"]])

    set.seed(seed)
    sims <- sample(nsim_stan, nsim)
  }

  report <- get_report(stanfit, fit, sims)
  data <- get_CMdata(fit)
  nsim_om <- length(sims)
  nyears <- data$Ldyr

  SRbeta <- sapply(report, getElement, "beta")
  phi <- sapply(report, getElement, "epro")

  matt <- sapply(report, getElement, "matt", simplify = "array") %>%
    aperm(3:1)
  mo <- sapply(report, getElement, "mo", simplify = "array") %>%
    aperm(3:1)
  Njuv <- sapply(report, getElement, "N", simplify = "array") %>%
    aperm(c(4, 2, 1, 3))
  Spawner <- sapply(1:2, function(i) {
    sapply(1:nsim_om, function(x) {
      esc <- report[[x]]$N[1:nyears, , i] * report[[x]]$sfish * report[[x]]$matt
      spawn <- esc * data$propwildspawn
      return(spawn)
    }, simplify = "array")
  }, simplify = "array") %>%
    aperm(c(3:1, 4))
  FPT <- sapply(1:2, function(...) sapply(report, getElement, "ft"), simplify = "array") %>%
    aperm(c(2, 1, 3))
  FT <- array(0, c(nsim_om, nyears, 2))

  Bio <- new(
    "Bio",
    nsim = nsim_om,
    maxage = data$Nages,
    p_mature = expand_array(matt, proyears),
    SRrel = "Ricker",
    kappa = as.numeric(pars$cr[sims]),
    Smax = 1/SRbeta,
    phi = phi,
    Mocean_NOS = expand_array(mo, proyears),
    fec = data$fec,
    p_female = data$ssum
  )

  Hatchery <- new(
    "Hatchery",
    Mocean_HOS = expand_array(mo, proyears),
    gamma = 1 #data$gamma
  )

  Harvest <- new(
    "Harvest",
    vulPT = data$vul
  )

  Habitat <- new("Habitat")

  Historical <- new(
    "Historical",
    HistSpawner = Spawner,
    HistN = Njuv[, , 1:nyears, ],
    HistFPT = FPT,
    HistFT = FT
  )

  SOM <- new("SOM", Bio, Hatchery, Habitat, Harvest, Historical, nyears = nyears, proyears = proyears)

  return(SOM)
}
