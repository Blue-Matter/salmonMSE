# Figures for conditioning model results

Functions used by the markdown report to generate summary figures from
the age-structured conditoning model

## Usage

``` r
get_report(stanfit, sims, inc_warmup = FALSE)

get_CMdata(fit)

CM_trace(stanfit, vars, inc_warmup = FALSE)

CM_pairs(stanfit, vars, inc_warmup = FALSE)

CM_fit_esc(report, d, year)

CM_fit_CWTesc(report, d, year1 = 1, rs_names)

CM_fit_CWTcatch(report, d, PT = TRUE, year1 = 1, rs_names)

CM_maturity(
  report,
  d,
  year1 = 1,
  r = 1,
  brood = TRUE,
  annual = FALSE,
  rs_names
)

CM_vul(report, type = c("vulPT", "vulT"))

CM_SRR(report, year1 = 1)

CM_prod(report, d, year1 = 1)

CM_Srep(report, d, year1 = 1, type = c("spawner", "egg"))

CM_M(report, year1 = 1, ci = TRUE)

CM_Megg(report, year1 = 1, ci = TRUE, surv = FALSE)

CM_Njuv(report, year1 = 1, ci = TRUE)

CM_recr(report, year1 = 1, ci = TRUE)

CM_esc(report, year1 = 1, ci = TRUE)

CM_F(report, PT = TRUE, year1 = 1, ci = TRUE)

CM_surv(report, year1 = 1, ci = TRUE)

CM_wt(stanfit, year1 = 1, ci = TRUE)

CM_surv2(report, year1 = 1, ci = TRUE, ylab = "Survival to age 2")

CM_wto(stanfit, year1 = 1, ci = TRUE)

CM_ER(
  report,
  brood = TRUE,
  type = c("PT", "T", "all"),
  year1 = 1,
  ci = TRUE,
  at_age = TRUE,
  r = 1
)

CM_CWT_ER(
  report,
  brood = TRUE,
  type = c("PT", "T", "all"),
  year1 = 1,
  ci = TRUE,
  rs_names
)

CM_covariate(x, names, year1 = 1, b, ylab = "Covariate")
```

## Arguments

- stanfit:

  Output from
  [`sample_CM()`](https://docs.salmonmse.com/reference/fit_CM.md)

- sims:

  Optional integer vector for subset of MCMC iterations

- inc_warmup:

  Logical, whether to include warmup MCMC samples

- fit:

  Output from `[fit_CM()]`

- vars:

  Character vector for variable names (see
  `names(stanfit@sim$samples[[1]])`). Regex and partial matching
  supported because it is passed to the `pattern` argument of
  [`grepl()`](https://rdrr.io/r/base/grep.html)

- report:

  List, output of state variables from individual MCMC samples, obtained
  with `get_report()`

- d:

  List of data variables, obtained with `get_CMdata()`

- year:

  Vector of years

- year1:

  Numeric, first year of model

- rs_names:

  Character vector of hatchery release strategies

- PT:

  Logical, whether to plot preterminal catch, otherwise (plot terminal
  catch)

- r:

  Integer, the release strategy for the figure (only if
  `annual = FALSE`)

- brood:

  Logical, whether to show results by brood year or return year (FALSE)

- annual:

  Logical, whether to show panel figure by individual year (TRUE) or a
  single time series figure

- type:

  Character, indicates type of variable to plot

- ci:

  Logical whether to show posterior intervals in addition to posterior
  median

- surv:

  Logical, whether to plot survival (values between 0 - 1) or
  instantaneous mortality rates

- ylab:

  Character y axis label

- at_age:

  Logical, whether to make figure by individual age

- x:

  Matrix of covariates by year x covariate

- names:

  Character of covariate names

- b:

  Matrix of fixed effect coefficients by simulation x covariate. If
  missing only the covariates (`x`) are plotted, otherwise, the dot
  product `sum(x * b)` is calculated by individual simulation and
  quantiles are plotted

## Value

- `get_report()` returns the list of state variables by individual MCMC
  samples

&nbsp;

- `get_CMdata()` returns the list of data variables used in the
  conditioning model

&nbsp;

- `CM_trace()` returns a ggplot showing the MCMC trace plot (aka
  wormplot)

&nbsp;

- `CM_pairs()` returns output from
  [`graphics::pairs()`](https://rdrr.io/r/graphics/pairs.html), a matrix
  of scatterplots of MCMC posterior samples

&nbsp;

- `CM_fit_esc()` returns base graphics with fit to total escapement time
  series

&nbsp;

- `CM_fit_CWTesc()` returns ggplot of fit to CWT escapement at age

&nbsp;

- `CM_fit_CWTcatch()` returns ggplot of fit to CWT catch at age

&nbsp;

- `CM_maturity()` returns ggplot of estimated maturity at age

&nbsp;

- `CM_vul()` returns ggplot of estimated fishery vulnerability at age

&nbsp;

- `CM_SRR()` returns ggplot of estimated stock-recruit relationship
  (density-dependent juvenile production from egg production) with
  average relationship and realized annual values. Years correspond to
  years of egg production (assumes juvenile production for the following
  brood year).

&nbsp;

- `CM_prod()` returns ggplot of realized productivity in the absence of
  fishery harvest, annual values are based on natural mortality and
  maturity at age

&nbsp;

- `CM_Srep()` returns ggplot of realized spawner or egg production at
  replacement, annual values are based on natural mortality and maturity
  at age

&nbsp;

- `CM_M()` returns ggplot of estimated natural mortality time series by
  age (marine stage)

&nbsp;

- `CM_Megg()` returns ggplot of egg-juvenile mortality time series

&nbsp;

- `CM_Njuv()` returns ggplot of juvenile abundance

&nbsp;

- `CM_recr()` returns ggplot of recruitment (mature return)

&nbsp;

- `CM_esc()` returns ggplot of escapement (after terminal harvest)

&nbsp;

- `CM_F()` returns ggplot of instantaneous fishing mortality

&nbsp;

- `CM_surv()` returns ggplot of natural survival (converting from
  instantaneous units of natural mortality)

&nbsp;

- `CM_wt()` returns ggplot of annual deviations in egg-juvenile
  mortality from the Ricker function

&nbsp;

- `CM_surv2()` returns ggplot of annual survival to age 2, which
  includes age-1 mortality (marine life stage) for both natural and
  hatchery origin fish. Hatchery fish experience additional mortality
  specified by release mortality.

&nbsp;

- `CM_wt()` returns ggplot of annual deviations in age 1 natural
  mortality (first year in marine life stage, deviations from time
  series average)

&nbsp;

- `CM_ER()` returns ggplot of exploitation rate either by individual age
  or aggregate values using adult equivalents

&nbsp;

- `CM_CWT_ER()` returns ggplot of CWT exploitation rate (by release
  strategy)

&nbsp;

- `CM_covariate()` returns ggplot of mortality covariates
