The current version of `salmonMSE` package is available on [CRAN](https://cran.r-project.org/package=salmonMSE).

# 0.1.1

## salmonMSE projections

- `SOM@Bio@p_female` parameter can now be age-specific. Intended for cases where older spawners are predominantly female. 
- Fix markdown reporting for projections
- If release target is exceeded in custom brood rule, projection will proportionally reduce broodtake to return to target
- Fix calculation in freshwater survival when habitat features are used

## Conditioning model

- Allow NAs in target population escapement time series in likelihood
- Conditioning model reports realized productivity annually with time-varying maturity and M
- Add initial hatchery production assumption
- Add likelihood for observed pHOS
- Export CM plotting functions

# 0.1.0

- Initial CRAN release
