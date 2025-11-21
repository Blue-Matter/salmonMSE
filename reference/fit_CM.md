# Fit conditioning model to historical data

Bayesian stock reconstruction model of natural and hatchery origin fish
population. Maturity and age-1 natural mortality are estimated from
coded wire tag catch and escapement at age. A separate series of
observed escapement, and hatchery releases reconstructs the population
of interest, informed by natural mortality and maturity from CWT
([Korman and Walters
2024](https://publications.gc.ca/site/eng/9.940685/publication.html)).
The model estimates time-varying maturity rate as well as time-varying
ocean survival as a linear model of covariates (separate covariates for
age 1 vs. ages 2+). The model can include either a preterminal juvenile
fishery, terminal return fishery, or both (see Data and start sections
of the documentation).`fit_CM()` generates the RTMB model from data
which can then be passed to `sample_CM()` to run the MCMC in Stan.
Generate a markdown report with
[`report_CM()`](https://docs.salmonmse.com/reference/report_CM.md).More
information is available on the
[salmonMSE](https://docs.salmonmse.com/articles/conditioning.html)
website

## Usage

``` r
fit_CM(
  data,
  start = list(),
  map = list(),
  lower = list(),
  upper = list(),
  do_fit = TRUE,
  silent = TRUE,
  control = list(eval.max = 1e+05, iter.max = 1e+05),
  ...
)

sample_CM(fit, ...)
```

## Arguments

- data:

  A list containing data inputs. See details.

- start:

  An optional list containing parameter starting values. See details.

- map:

  An optional list that describes how parameters are fixed in the model.
  See [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

- lower:

  Named list containing lower bounds for parameters. See details.

- upper:

  Named list containing upper bounds for parameters. See details.

- do_fit:

  Logical, whether to do the fit and estimate the Hessian.

- silent:

  Logical, whether to silence output from RTMB to the console.

- control:

  List, `control` argument to pass to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- ...:

  For `fit_CM`, arguments to
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).
  For `sample_CM`, arguments to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

- fit:

  List of output from `fit_CM()`

## Value

- `fit_CM()` returns a named list containing the RTMB model (`obj`),
  nlminb output (`opt`), standard errors (`SD`), and parameter bounds
  (`lower` and `upper`)

- `sample_CM()` returns a `stanfit` object containing the MCMC chains

## Data

Data should passed through a named list with the following entries.

- `Nages` Integer, number of age classes in the model

- `Ldyr` Integer, number of years in the model

- `lht` Integer, life history type. Should be 1 for now

- `n_r` Integer, number of release strategies for CWT, subset of a
  hatchery-origin brood year that differ in maturity rate. Default is 1.

- `cwtrelease` Matrix `[Ldyr, n_r]`, coded wire tag (CWT) releases by
  year and release strategy

- `cwtesc` Array `[Ldyr, Nages, n_r]`. CWT escapement **by brood year,
  age, and release strategy**. Poisson likelhood.

- `cwtcatPT` Array `[Ldyr, Nages, n_r]`. CWT preterminal catch (juvenile
  fish), **by brood year, age, and release strategy**. Poisson
  likelhood. Set all values to zero to turn off parameters related to
  the preterminal fishery.

- `cwtcatT` Array `[Ldyr, Nages, n_r]`. CWT terminal catch (returning,
  mature fish), **by brood year, age, and release strategy**. Poisson
  likelhood. Set all values to zero to turn off parameters related to
  the terminal fishery.

- `bvulPT` Vector length `Nages`. Prior mean for the vulnerability at
  age to the preterminal fishery.

- `bvulT` Vector length `Nages`. Prior mean for the vulnerability at age
  to the terminal fishery.

- `RelRegFPT` Vector `Ldyr`. Trend in relative regional preterminal
  fishing mortality. Fishing mortality is estimated by estimating a
  scaling coefficient and annual deviations from this vector. Default is
  `rep(1, d$Ldyr)` (no prior trend) if `cwtcatPT` is provided, otherwise
  zero.

- `RelRegFT` Vector `Ldyr`. Trend in relative regional terminal fishing
  mortality. Default is `rep(1, d$Ldyr)` (no prior trend) if `cwtcatT`
  is provided, otherwise zero.

- `bmatt` Vector length `Nages`. Proportion maturity at age, base values
  for calculating the unfished replacement line. Also the prior means if
  year-specific maturity rates are estimated.

- `mobase`. Vector length `Nages`. Natural mortality at age, base values
  for calculating the unfished replacement line and the the equilibrium
  spawners at age.

- `covariate1` *Optional*. Matrix `Ldyr, ncov1` of linear covariates
  that predict natural mortality for age 1.

- `covariate` *Optional*. Matrix `Ldyr, ncov` of linear covariates that
  predict natural mortality for ages 2+.

- `hatchsurv` Numeric, survival of hatchery releases into the smolt life
  stage. Density-independent. Default is 1. If less than 1, then
  hatchery origin fish have lower survival to age 2 (after first year of
  marine life stage) compared to natural origin fish.

- `gamma` *Optional*. Numeric, the relative spawning success of hatchery
  origin spawners. Default is 1.

- `ssum` Numeric, proportion of spawners that is female. Can also be a
  vector `Nages`

- `fec` Vector length `Nages`. Fecundity, egg production at age

- `r_matt` Integer, the release strategy for which to use maturity
  parameter for the natural system. Default is 1.

- `obsescape` Vector length `Ldyr`, total observed escapement from
  fisheries, i.e., return to river (all ages and both hatchery/natural
  fish). Lognormal likelhood.

- `propwildspawn` Vector length `Ldyr`, proportion of the escapement
  that spawn (accounts for en-route mortality and broodtake)

- `obs_pHOS` Vector length `Ldyr`, proportion of hatchery origin
  spawners (census) (between 0-1). Logistic-normal likelihood.

- `pHOS_sd` Numeric, logistic-normal standard deviation. Default is 0.1.

- `hatch_init` Numeric, equilibrium hatchery releases used to initialize
  the model. Helpful if hatchery production starts prior to CWT time
  series. Default is zero.

- `s_enroute` Numeric, survival of escapement to spawning grounds.
  Default is 1.

- `so_mu` Numeric, the prior mean for spawners at unfished replacement
  in logspace. Default is `log(3 * max(data$obsescape))`.

- `so_sd` Numeric, the prior standard deviation for spawners at unfished
  replacement in logspace. Default is 0.5.

- `finitPT` Numeric, initial preterminal fishing mortality for
  calculating the equilibrium spawners at age in the first year of the
  model. Default is 0. Set to `"estimate"` to allow the model to
  estimate the equilibrium condition.

- `finitT` Numeric, initial terminal fishing mortality for calculating
  the equilibrium spawners at age in the first year of the model.
  Default is 0. Set to `"estimate"` to allow the model to estimate the
  equilibrium condition.

- `cwtExp` Numeric, the CWT expansion factor, typically the reciprocal
  of the catch sampling rate (higher factors for lower sampling rate).
  The model scales down the CWT predictions to match the observations.
  In other words, the model assumes that the CWT catch and escapement
  are not expanded. For example, `cwtExp = 10` divides the CWT
  predictions by 10 for the likelihood. Default is 1. The Poisson
  distribution is used for the likelihood of the CWT observations, and
  the expansion parameter can be used to downweight the CWT likelihood
  relative to the escapement time series. However it requires
  adjustments of the CWT catches prior to fitting to ensure the proper
  population scale. If the expanded catch is 100, then the input CWT
  catch should be 10 and 50 with `cwtExp` of 10 and 2, respectively, to
  maintain the same population scale. The Poisson variance scales with
  the mean and is higher with `cwtExp = 2`.

- `fitness` Logical, whether to calculate fitness effects on survival.
  Default is `FALSE`.

- `theta` Vector length 2, the optimum phenotype value for the natural
  and hatchery environments. Default is 100 and 80, respectively. See
  [online
  article](https://docs.salmonmse.com/articles/equations.html#fitness-effects-on-survival)
  for more information.

- `rel_loss` Vector length 3, the loss in fitness apportioned between
  the egg, fry (both prior to density-dependence), and smolt (after
  density-dependence) life stage. The three values should sum to 1.

- `zbar_start` Vector length 2, the mean phenotype of the spawners and
  broodtake in the natural and hatchery environment, respectively, at
  the start of the model. Default values of 100 and 100, implying
  maximum fitness at for the natural environment at the start of the
  model.

- `fitness_variance` Numeric. The variance (omega-squared) of the
  fitness function. Assumed identical between the natural and hatchery
  environments. Default is 100.

- `phenotype_variance` Numeric. The variance (sigma-squared) of the
  phenotypic trait (with optimum theta). Assumed identical between the
  natural and hatchery environments. Default is 10.

- `heritability` Numeric. The heritability (h-squared) of the phenotypic
  trait. Between 0-1. Default is 0.5.

- `fitness_floor` Numeric. The minimum fitness value in the natural and
  hatchery environments. Default is 0.5.

## start

Starting values for parameters can be provided through a named list:

- `log_cr` Numeric, log of the compensation ratio (productivity).
  Default is 3.

- `log_so` Numeric, unfished spawners in logspace. Default is
  `log(3 * max(data$obsescape))`.

- `moadd` Numeric, additive term to base natural mortality rate for age
  1 juveniles. Default is zero.

- `wt` Vector `Ldyr`. Annual deviates in natural mortality during the
  freshwater life stage (affects egg to smolt survival). Estimated with
  normal prior with mean zero and standard deviation `p$wt_sd`. Default
  is zero.

- `wto` Vector `Ldyr`. Annual deviates in natural mortality for age 1
  juveniles (marine life stage). Estimated with normal prior with mean
  zero and standard deviation `p$wto_sd`. Default is zero.

- `log_FbasePT` Numeric, scaling coefficient to estimate preterminal
  fishing mortality from `data$RelRegFPT`. Default is `log(0.1)`.

- `log_FbaseT` Numeric, scaling coefficient to estimate preterminal
  fishing mortality from `data$RelRegFT`. Default is `log(0.1)`.

- `log_fanomalyPT` Vector `Ldyr`. Annual lognormal deviates from
  `exp(log_FbasePT) * data$RelRegFPT` to estimate preterminal fishing
  mortality. Estimated with normal prior with mean zero and standard
  deviation `p$fanomaly_sd`. Default is zero.

- `log_fanomalyT` Vector `Ldyr`. Annual lognormal deviates from
  `exp(log_FbaseT) * data$RelRegFT` to estimate terminal fishing
  mortality. Estimated with normal prior with mean zero and standard
  deviation `p$fanomalyPT_sd`. Default is zero.

- `lnE_sd` Numeric, lognormal standard deviation of the observed
  escapement. Estimated with hierarchical `gamma(2, 5)` prior. Default
  is 0.1.

- `wt_sd` Numeric, lognormal standard deviation of the egg to smolt
  (freshwater) natural mortality deviates. Estimated with hierarchical
  `gamma(2, 5)` prior. Default is 1.

- `wto_sd` Numeric, lognormal standard deviation of the age 1 (marine)
  natural mortality deviates. Estimated with hierarchical `gamma(2, 5)`
  prior. Default is 1.

- `fanomalyPT_sd` Numeric, lognormal standard deviation of `fanomalyPT`.
  Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.

- `fanomalyT_sd` Numeric, lognormal standard deviation of `fanomalyT`.
  Estimated with hierarchical `gamma(2, 5)` prior. Default is 1.

- `logit_vulPT` Vector `Nages-2` of preterminal vulnerability at age in
  logit space. Fixed to zero and one at age 1 and the maximum age,
  respectively. Default is `qlogis(data$bvul_PT[-c(1, data$Nages)])`.

- `logit_vulT` Vector `Nages-2` of terminal vulnerability at age in
  logit space. Fixed to zero and one at age 1 and the maximum age,
  respectively. Default is `qlogis(data$bvul_T[-c(1, data$Nages)])`.

- `logit_matt` Matrix `Ldyr, Nages-2` maturity by year and age in logit
  space. Maturity is fixed to zero and one at age 1 and the maximum age,
  respectively. Default is
  `matrix(qlogis(data$bmatt[-c(1, data$Nages)]), data$Ldyr, data$Nages-2, byrow = TRUE)`.

- `sd_matt` Vector `Nages-2`. Logit standard deviation of maturity
  (`logit_matt`) by age class. Default is 0.5.

- `b1` Vector `ncov1` of coefficients for linear covariates that predict
  natural mortality for age 1. Default is zero.

- `b` Vector `ncov` of coefficients for linear covariates that predict
  natural mortality for ages 2+. Default is zero.

## Bounds

By default, the standard deviation parameters and parameters in normal
space (e.g., `FbasePT`, `Fbase_T`) have a lower bound of zero. `moadd`
has a lower bound of zero by default, but it is feasible that this
parameter can be negative as well. Deviation parameters centred around
zero are bounded between -3 to 3. The `log_cr` parameter has a lower
bound of zero.

All other parameters are unbounded.

## Covariates on natural mortality

Natural mortality is modeled as the sum of a base value
\\M^\textrm{base}\\, additional scaling factor for age 1
\\M^\textrm{add}\\, a linear system of covariates \\X\\ and coefficients
\\b\\:

\$\$ M\_{y,a} = \begin{cases} M^\textrm{base}\_a + M^\textrm{add} +
\sum_j b^1_j X^1\_{y,j} & \quad a = 1\\ M^\textrm{base}\_a + \sum_j b_j
X\_{y,j} & \quad a = 2, \ldots, A \end{cases} \$\$

## References

Korman, J. and Walters, C. 2024. A life cycle model for Chinook salmon
population dynamics. Canadian Contractor Report of Hydrography and Ocean
Sciences 62: vi + 60 p.

## See also

[`report_CM()`](https://docs.salmonmse.com/reference/report_CM.md)

[`CM2SOM()`](https://docs.salmonmse.com/reference/CM2SOM.md)

## Author

Q. Huynh from Stan code provided by J. Korman and C. Walters
