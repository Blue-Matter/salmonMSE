# Changelog

## Version 3.0.0

CRAN release: 2026-08-20

- salmonMSE no longer uses openMSE as a dependency.
  [`ProjectSOM()`](https://docs.salmonmse.com/reference/salmonMSE.md) is
  now the primary internal function that organizes the projection.
  Various internal functions have been replaced and updated.
- Parallel processing by individual operating model is now supported
  with the `ncores` argument to
  [`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md).
- Check for NAs in start values for `fit_CM`.
- Update
  [`.CM_MSY()`](https://docs.salmonmse.com/reference/dot-CM_prod.md) and
  `.CM_Sgen()` functions for parallel computation.
- MSY calculations and catch calculations use adult equivalents (AEQs)
  for preterminal fisheries.
- Add new slot `Hatchery@p_female_brood` for sex-selective brood.
- `SMSE` output object reports brood numbers and catch with age
  dimension. Harvest rates and exploitation rates are aggregated over
  age classes.
- Add custom catch function (for harvest control rules, etc.)
- Add `simple_SOM` class and
  [`simple_salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md)
  function for modeling recruits-spawners

## Version 2.1.0

CRAN release: 2026-04-25

- Update terminal year escapement, brood, and egg production
  calculations for multiple release strategies and life cycle groups
- Fix typos in `smolt_func()` and
  [`compare_statevar_ts()`](https://docs.salmonmse.com/reference/compare_statevar_ts.md)
- Broodtake optimization uses
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) with tolerance to
  machine precision (need this when population is very large and
  broodtake is small)
- `plot_tradeoff` can plot `x1` as a continuous variable (numerics were
  previously coerced to factors)

## Version 2.0.0

CRAN release: 2026-04-16

- Add example objects
- New arguments in Hatchery object to determine if hatchery releases
  compete with natural-origin juveniles
- [`plot_tradeoff()`](https://docs.salmonmse.com/reference/plot_tradeoff.md)
  can create panels with scenario rows and columns
- [`plot_decision_table2()`](https://docs.salmonmse.com/reference/plot_decision_table.md)
  allows more control of colour scheme. Add examples as well in
  documentation
- Update units of Ricker Smax to be units of spawners, `tau` argument
  convert to Emax (corresponding egg production)
- Projection now reports terminal year escapement, brood, and egg
  production (long-standing issue)

## Version 1.0.0

CRAN release: 2026-03-20

### salmonMSE projections

- `SOM@Bio@p_female` parameter can now be age-specific. Intended for
  cases where older spawners are predominantly female.
- Fix markdown reporting for projections
- If release target is exceeded in custom brood rule, projection will
  proportionally reduce broodtake to return to target
- Fix calculation in freshwater survival when habitat features are used
- Hard code in openMSE historical period to 2 real years. Historical
  reconstruction no longer supported.
- M vector now reduce to length `maxage-1`
- Allow NO removals at the spawning grounds (facilitates in-river
  fishery operations)
- Add `compare` function to generate markdown report to compare several
  model runs

### Conditioning model

- Allow NAs in target population escapement time series in likelihood
- Conditioning model reports realized productivity (adults/spawner)
  annually with time-varying maturity and M
- Add initial hatchery production assumption
- Add likelihood for observed pHOS (by calendar year or brood year)
- Export CM plotting functions
- Allow estimation of finitPT and finitT (initial fishing mortality),
  but likely slows down MCMC convergence
- Allow initial population to be specified by `spawn_init` (total
  spawners) and `pHOS_init` (fraction hatchery-origin)

### Other

- Decision tables and tradeoff figures can now be a grid
- Remove [`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html) import
  for compatibility with RTMB 1.9

## Version 0.1.0

CRAN release: 2025-09-23

- Initial CRAN release
