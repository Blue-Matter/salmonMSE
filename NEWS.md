The current version of `salmonMSE` package is available on [CRAN](https://cran.r-project.org/package=salmonMSE).

# 2.1.0

- Update terminal year escapement, brood, and egg production calculations for multiple release strategies and life cycle groups
- Fix typos in `smolt_func()` and `compare_statevar_ts()`
- Broodtake optimization uses `uniroot()` with tolerance to machine precision (need this when population is very large and broodtake is small)
- `plot_tradeoff` can plot `x1` as a continuous variable (numerics were previously coerced to factors)

# 2.0.0

- Add example objects
- New arguments in Hatchery object to determine if hatchery releases compete with natural-origin juveniles
- `plot_tradeoff()` can create panels with scenario rows and columns
- `plot_decision_table2()` allows more control of colour scheme. Add examples as well in documentation
- Update units of Ricker Smax to be units of spawners, `tau` argument convert to Emax (corresponding egg production)
- Projection now reports terminal year escapement, brood, and egg production (long-standing issue)

# 1.0.0

## salmonMSE projections

- `SOM@Bio@p_female` parameter can now be age-specific. Intended for cases where older spawners are predominantly female. 
- Fix markdown reporting for projections
- If release target is exceeded in custom brood rule, projection will proportionally reduce broodtake to return to target
- Fix calculation in freshwater survival when habitat features are used
- Hard code in openMSE historical period to 2 real years. Historical reconstruction no longer supported.
- M vector now reduce to length `maxage-1`
- Allow NO removals at the spawning grounds (facilitates in-river fishery operations)
- Add `compare` function to generate markdown report to compare several model runs

## Conditioning model

- Allow NAs in target population escapement time series in likelihood
- Conditioning model reports realized productivity (adults/spawner) annually with time-varying maturity and M
- Add initial hatchery production assumption
- Add likelihood for observed pHOS (by calendar year or brood year)
- Export CM plotting functions
- Allow estimation of finitPT and finitT (initial fishing mortality), but likely slows down MCMC convergence
- Allow initial population to be specified by `spawn_init` (total spawners) and `pHOS_init` (fraction hatchery-origin)

## Other

- Decision tables and tradeoff figures can now be a grid
- Remove `stats::uniroot` import for compatibility with RTMB 1.9

# 0.1.0

- Initial CRAN release
