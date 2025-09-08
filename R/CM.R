


#' Fit conditioning model to historical data
#'
#' @description{
#' Bayesian stock reconstruction model of natural and hatchery origin fish population.
#' Maturity and age-1 natural mortality are estimated from coded wire tag catch and escapement at age.
#' A separate series of observed escapement, and hatchery releases reconstructs the population of interest,
#' informed by natural mortality and maturity from CWT ([Korman and Walters 2024](https://publications.gc.ca/site/eng/9.940685/publication.html)).
#' The model estimates time-varying maturity rate as well as time-varying ocean survival as a linear model of covariates (separate covariates
#' for age 1 vs. ages 2+).
#' The model can include either a preterminal juvenile fishery, terminal return fishery, or both (see Data and start sections of the documentation).
#'
#' `fit_CM()` generates the RTMB model from data which can then be passed to `sample_CM()` to run the MCMC in Stan.
#'
#' More information is available on the [salmonMSE](https://docs.salmonmse.com/articles/conditioning.html) website
#' }
#'
#' @param data  A list containing data inputs. See details.
#' @param start An optional list containing parameter starting values. See details.
#' @param map An optional list that describes how parameters are fixed in the model. See [`TMB::MakeADFun()`].
#' @param lower Named list containing lower bounds for parameters. See details.
#' @param upper Named list containing upper bounds for parameters. See details.
#' @param do_fit Logical, whether to do the fit and estimate the Hessian.
#' @param silent Logical, whether to silence output from RTMB to the console.
#' @param control List, `control` argument to pass to [`stats::nlminb()`].
#' @param ... For `fit_CM`, arguments to [`RTMB::MakeADFun()`]. For `sample_CM`, arguments to `rstan::sampling()`
#' @import RTMB
#' @importFrom stats nlminb
#' @returns
#' - `fit_CM()` returns a named list containing the RTMB model (`obj`), nlminb output (`opt`), standard errors (`SD`),
#' and parameter bounds (`lower` and `upper`)
#' - `sample_CM()` returns a `stanfit` object containing the MCMC chains
#'
#' @section Data:
#' Data should passed through a named list with the following entries.
#'
#' - `Nages` Integer, number of age classes in the model
#' - `Ldyr` Integer, number of years in the model
#' - `lht` Integer, life history type. Should be 1 for now
#' - `n_r` Integer, number of release strategies for CWT, subset of a hatchery-origin brood year that differ in maturity rate. Default is 1.
#' - `cwtrelease` Matrix `[Ldyr, n_r]`, coded wire tag (CWT) releases by year and release strategy
#' - `cwtesc` Array `[Ldyr, Nages, n_r]`. CWT escapement **by brood year, age, and release strategy**. Poisson likelhood.
#' - `cwtcatPT` Array `[Ldyr, Nages, n_r]`. CWT preterminal catch (juvenile fish), **by brood year, age, and release strategy**. Poisson likelhood. Set all values to zero to turn off
#' parameters related to the preterminal fishery.
#' - `cwtcatT` Array `[Ldyr, Nages, n_r]`. CWT terminal catch (returning, mature fish), **by brood year, age, and release strategy**. Poisson likelhood. Set all values to zero to turn off
#' parameters related to the terminal fishery.
#'
#' - `bvulPT` Vector length `Nages`. Prior mean for the vulnerability at age to the preterminal fishery.
#' - `bvulT` Vector length `Nages`. Prior mean for the vulnerability at age to the terminal fishery.
#'
#' - `RelRegFPT` Vector `Ldyr`. Trend in relative regional preterminal fishing mortality. Fishing mortality is estimated by estimating a scaling
#' coefficient and annual deviations from this vector.
#' - `RelRegFT` Vector `Ldyr`. Trend in relative regional terminal fishing mortality.
#'
#' - `bmatt` Vector length `Nages`. Proportion maturity at age, base values for calculating the unfished replacement line. Also the prior means if year-specific
#' maturity rates are estimated.
#'
#' - `mobase`. Vector length `Nages`. Natural mortality at age, base values for calculating the unfished replacement line and the
#' the equilibrium spawners at age.
#'
#' - `covariate1` *Optional*. Matrix `Ldyr, ncov1` of linear covariates that predict natural mortality for age 1.
#' - `covariate` *Optional*. Matrix `Ldyr, ncov` of linear covariates that predict natural mortality for ages 2+.
#'
#' - `hatchsurv` Numeric, survival of hatchery releases into the smolt life stage. Density-independent.
#' - `gamma` *Optional*. Numeric, the relative spawning success of hatchery origin spawners. Default is 1.
#' - `ssum` Numeric, proportion of spawners that is female
#'
#' - `fec` Vector length `Nages`. Fecundity, egg production at age
#'
#' - `r_matt` Integer, the release strategy for which to use maturity parameter for the natural system. Default is 1.
#' - `obsescape` Vector length `Ldyr`, total observed escapement (all ages and both hatchery/natural fish). Lognormal likelhood.
#' - `propwildspawn` Vector length `Ldyr`, proportion of the escapement that spawn (accounts for en-route mortality and broodtake)
#' - `hatchrelease` Vector length `Ldyr+1`, number of hatchery juvenile fish released
#' - `s_enroute` Numeric, survival of escapement to spawning grounds. Default is 1.
#'
#' - `so_mu` Numeric, the prior mean for unfished spawners in logspace. Default is `log(3 * max(data$obsescape))`.
#' - `so_sd` Numeric, the prior standard deviation for unfished spawners in logspace. Default is 0.5.
#'
#' - `finitPT` Numeric, initial preterminal fishing mortality for calculating the equilibrium spawners at age in the first year of the model. Default is 0.
#' - `finitT` Numeric, initial terminal fishing mortality for calculating the equilibrium spawners at age in the first year of the model. Default is 0.
#'
#' - `cwtExp` Numeric, the CWT expansion factor, typically the reciprocal of the catch sampling rate (higher factors for lower sampling rate).
#' The model scales down the CWT predictions to match the observations. In other words,
#' the model assumes that the CWT catch and escapement are not expanded. For example, `cwtExp = 10` divides the CWT predictions by 10 for the likelihood. Default is 1.
#' The Poisson distribution is used for the likelihood of the CWT observations, and the expansion parameter can be used to downweight the CWT likelihood relative to the escapement time series.
#' However it requires adjustments of the CWT catches prior to fitting to ensure the proper population scale.
#' If the expanded catch is 100, then the input CWT catch should be 10 and 50 with `cwtExp` of 10 and 2, respectively, to maintain the same population scale.
#' The Poisson variance scales with the mean and is higher with `cwtExp = 2`.
#' - `fitness` Logical, whether to calculate fitness effects on survival. Default is `FALSE`.
#' - `theta` Vector length 2, the optimum phenotype value for the natural and hatchery environments. Default is 100 and 80, respectively. See
#' [online article](https://docs.salmonmse.com/articles/equations.html#fitness-effects-on-survival) for more information.
#' - `rel_loss` Vector length 3, the loss in fitness apportioned between the egg, fry (both prior to density-dependence), and smolt (after density-dependence) life stage. The three values should sum to 1.
#' - `zbar_start` Vector length 2, the mean phenotype of the spawners and broodtake in the natural and hatchery environment, respectively, at the start of the model. Default values of 100 and 100, implying maximum fitness at
#' for the natural environment at the start of the model.
#' - `fitness_variance` Numeric. The variance (omega-squared) of the fitness function. Assumed identical between the natural and hatchery environments. Default is 100.
#' - `phenotype_variance` Numeric. The variance (sigma-squared) of the phenotypic trait (with optimum theta). Assumed identical between the natural and hatchery environments. Default is 10.
#' - `heritability` Numeric. The heritability (h-squared) of the phenotypic trait. Between 0-1. Default is 0.5.
#' - `fitness_floor` Numeric. The minimum fitness value in the natural and hatchery environments. Default is 0.5.
#'
#' @section start:
#' Starting values for parameters can be provided through a named list:
#'
#' - `log_cr` Numeric, log of the compensation ratio (productivity). Default is 3.
#' - `log_so` Numeric, unfished spawners in logspace. Default is `log(3 * max(data$obsescape))`.
#' - `moadd` Numeric, additive term to base natural mortality rate for age 1 juveniles. Default is zero.
#' - `wt` Vector `Ldyr`. Annual deviates in natural mortality during the freshwater life stage (affects egg to smolt survival).
#' Estimated with normal prior with mean zero and standard deviation `p$wt_sd`. Default is zero.
#' - `wto` Vector `Ldyr`. Annual deviates in natural mortality for age 1 juveniles (marine life stage).
#' Estimated with normal prior with mean zero and standard deviation `p$wto_sd`. Default is zero.
#' - `log_FbasePT` Numeric, scaling coefficient to estimate preterminal fishing mortality from `data$RelRegFPT`. Default is `log(0.1)`.
#' - `log_FbaseT` Numeric, scaling coefficient to estimate preterminal fishing mortality from `data$RelRegFT`. Default is `log(0.1)`.
#' - `log_fanomalyPT` Vector `Ldyr`. Annual lognormal deviates from `exp(log_FbasePT) * data$RelRegFPT` to estimate preterminal fishing mortality.
#' Estimated with normal prior with mean zero and standard deviation `p$fanomaly_sd`. Default is zero.
#' - `log_fanomalyT` Vector `Ldyr`. Annual lognormal deviates from `exp(log_FbaseT) * data$RelRegFT` to estimate terminal fishing mortality.
#' Estimated with normal prior with mean zero and standard deviation `p$fanomalyPT_sd`. Default is zero.
#' - `lnE_sd` Numeric, lognormal standard deviation of the observed escapement. Estimated with hierarchical `gamma(2, 5)` prior. Default is 0.1.
#' - `wt_sd` Numeric, lognormal standard deviation of the egg to smolt (freshwater) natural mortality deviates. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `wto_sd` Numeric, lognormal standard deviation of the age 1 (marine) natural mortality deviates. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `fanomalyPT_sd` Numeric, lognormal standard deviation of `fanomalyPT`. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `fanomalyT_sd` Numeric, lognormal standard deviation of `fanomalyT`. Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.
#' - `logit_vulPT` Vector `Nages-2` of preterminal vulnerability at age in logit space. Fixed to zero and one at age 1 and the maximum age, respectively.
#' Default is `qlogis(data$bvul_PT[-c(1, data$Nages)])`.
#' - `logit_vulT` Vector `Nages-2` of terminal vulnerability at age in logit space. Fixed to zero and one at age 1 and the maximum age, respectively.
#' Default is `qlogis(data$bvul_T[-c(1, data$Nages)])`.
#' - `logit_matt` Matrix `Ldyr, Nages-2` maturity by year and age in logit space. Maturity is fixed to zero and one at age 1 and the maximum age, respectively.
#' Default is `matrix(qlogis(data$bmatt[-c(1, data$Nages)]), data$Ldyr, data$Nages-2, byrow = TRUE)`.
#' - `sd_matt` Vector `Nages-2`. Logit standard deviation of maturity (`logit_matt`) by age class. Default is 0.5.
#' - `b1` Vector `ncov1` of coefficients for linear covariates that predict natural mortality for age 1. Default is zero.
#' - `b` Vector `ncov` of coefficients for linear covariates that predict natural mortality for ages 2+. Default is zero.
#'
#' @section Bounds:
#'
#' By default, the standard deviation parameters and parameters in normal space (e.g., `FbasePT`, `Fbase_T`) have a lower bound of zero.
#' `moadd` has a lower bound of zero by default, but it is feasible that this parameter can be negative as well.
#' Deviation parameters centred around zero are bounded between -3 to 3.
#' The `log_cr` parameter has a lower bound of zero.
#'
#' All other parameters are unbounded.
#'
#' @section Covariates on natural mortality:
#'
#' Natural mortality is modeled as the sum of a base value \eqn{M^\textrm{base}}, additional scaling factor for age 1 \eqn{M^\textrm{add}},
#' a linear system of covariates \eqn{X} and coefficients \eqn{b}:
#'
#' \deqn{
#' M_{y,a} =
#' \begin{cases}
#' M^\textrm{base}_a + M^\textrm{add} + \sum_j b^1_j X^1_{y,j} & \quad a = 1\\
#' M^\textrm{base}_a + \sum_j b_j X_{y,j} & \quad a = 2, \ldots, A
#' \end{cases}
#' }
#' @author Q. Huynh with Stan code provided by J. Korman and C. Walters
#' @references
#' Korman, J. and Walters, C. 2024. A life cycle model for Chinook salmon population dynamics. Canadian Contractor Report of Hydrography and Ocean
#' Sciences 62: vi + 60 p.
#' @seealso [CM2SOM()]
#' @export
fit_CM <- function(data, start = list(), map = list(), lower = list(), upper = list(), do_fit = TRUE, silent = TRUE,
                   control = list(eval.max = 1e5, iter.max = 1e5), ...) {

  data <- check_data(data)
  p <- make_CMpars(start, data)
  map <- make_map(map, p, data)

  f <- function(p) salmonMSE::CM_int(p, d = data) # :: is needed for parallel MCMC sampling
  obj <- RTMB::MakeADFun(func = f, parameters = p, map = map, silent = silent, ...)

  bounds <- make_bounds(names(obj$par), data, lower, upper)

  if (do_fit) {
    opt <- nlminb(
      obj$par, obj$fn, obj$gr,
      lower = bounds$lower, upper = bounds$upper, control = control
    )
    SD <- RTMB::sdreport(obj)
  } else {
    opt <- SD <- NULL
  }

  res <- list(obj = obj, opt = opt, SD = SD, lower = bounds$lower, upper = bounds$upper)
  return(res)
}

#' @rdname fit_CM
#' @param fit List of output from `fit_CM()`
#' @export
sample_CM <- function(fit, ...) {
  if (!requireNamespace("tmbstan", quietly = TRUE)) stop("tmbstan package is needed.")
  samp <- tmbstan::tmbstan(fit$obj, lower = fit$lower, upper = fit$upper, ...)
  samp@.MISC$CMfit <- fit
  return(samp)
}


#' Convert conditioning model to operating model
#'
#' Creates an operating model from MCMC samples and data inputs of the conditioning model.
#' Management actions for habitat, hatchery production, and harvest still need to be specified in the operating model.
#'
#' @param stanfit Output from [`sample_CM()`]
#' @param sims Optional, a vector of integers indicating the MCMC iterations to convert to operating model simulations. Otherwise,
#' use argument `nsim` in order to sample a subset of the MCMC.
#' @param nsim Integer, total number of simulations in the operating model. Only used if `sims` is missing.
#' @param seed Integer, seed for sampling the MCMC output. Only used if `sims` is missing.
#' @param proyears Integer, the number of projection years in the operating model
#' @return \linkS4class{SOM} object.
#' @importFrom abind abind
#' @export
CM2SOM <- function(stanfit, sims, nsim = 2, seed = 1, proyears = 40) {

  if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan package is needed.")
  pars <- rstan::extract(stanfit)

  if (missing(sims)) {
    nsim_stan <- length(pars[["lp__"]])

    set.seed(seed)
    sims <- sample(nsim_stan, nsim)
  }

  fit <- stanfit@.MISC$CMfit
  if (is.null(fit)) stop("CM fitted object not found in stanfit@.MISC$CMfit")
  data <- get_CMdata(fit)

  report <- get_report(stanfit, sims)
  nsim_om <- length(sims)
  nyears <- data$Ldyr

  SRbeta <- sapply(report, getElement, "beta")
  phi <- sapply(report, getElement, "epro")

  matt <- sapply(report, getElement, "matt", simplify = "array") %>%
    aperm(3:1)

  mo <- sapply(report, getElement, "mo", simplify = "array") %>%
    aperm(3:1)
  mo_maxage <- mo[, dim(mo)[2], , drop = FALSE]
  mo <- abind::abind(mo, mo_maxage, along = 2)

  Njuv <- sapply(report, getElement, "N", simplify = "array") %>%
    aperm(c(4, 2, 1, 3))
  Spawner <- sapply(report, getElement, "syear", simplify = "array") %>%
    aperm(c(4, 2, 1, 3))
  FPT <- sapply(1:2, function(...) sapply(report, getElement, "FPT"), simplify = "array") %>%
    aperm(c(2, 1, 3))
  FT <- sapply(1:2, function(...) sapply(report, getElement, "FT"), simplify = "array") %>%
    aperm(c(2, 1, 3))

  Bio <- new(
    "Bio",
    maxage = data$Nages,
    n_g = 1,
    p_LHG = 1,
    p_mature = expand_array(matt, proyears),
    SRrel = "Ricker",
    kappa = as.numeric(exp(pars$log_cr[sims])),
    Smax = 1/SRbeta,
    phi = phi,
    Mjuv_NOS = expand_array(mo, proyears),
    fec = data$fec,
    p_female = data$ssum,
    s_enroute = data$s_enroute
  )

  Hatchery <- new(
    "Hatchery",
    Mjuv_HOS = expand_array(mo, proyears),
    gamma = data$gamma
  )

  Harvest <- new(
    "Harvest",
    vulPT = sapply(report, getElement, "vulPT") %>% t(),
    vulT =sapply(report, getElement, "vulT") %>% t()
  )

  Habitat <- new("Habitat")

  Historical <- new(
    "Historical",
    HistSpawner_NOS = Spawner[, , , 1],
    HistSpawner_HOS = Spawner[, , , 2],
    HistNjuv_NOS = Njuv[, , , 1, drop = FALSE],
    HistNjuv_HOS = Njuv[, , , 2, drop = FALSE],
    HistFPT = FPT,
    HistFT = FT
  )

  if (data$fitness) {
    Hatchery@fitness_type <- c("Ford", "none")
    Hatchery@theta <- data$theta
    Hatchery@rel_loss <- data$rel_loss

    zbar <- sapply(report, getElement, "zbar", simplify = "array") %>%
      aperm(c(3, 1, 2))
    Hatchery@zbar_start <- zbar[, nyears - seq(1, data$Nages) + 1, ]

    Hatchery@phenotype_variance <- data$phenotype_variance
    Hatchery@fitness_variance <- data$fitness_variance
    Hatchery@heritability <- data$heritability
    Hatchery@fitness_floor <- data$fitness_floor
  } else {
    Hatchery@fitness_type <- c("none", "none")
  }

  SOM <- new("SOM", Bio, Habitat, Hatchery, Harvest, Historical,
             nsim = nsim_om, nyears = nyears, proyears = proyears)

  return(SOM)
}
