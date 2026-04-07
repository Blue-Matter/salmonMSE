The current version of `salmonMSE` package is available on [CRAN](https://cran.r-project.org/package=salmonMSE).

# 1.1.0

- Example objects
- New arguments in Hatchery object to determine if hatchery releases compete with natural-origin juveniles

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
