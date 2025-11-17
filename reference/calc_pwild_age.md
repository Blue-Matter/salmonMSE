# Proportion wild spawners

Calculate the proportion of wild spawners from a time series of spawners

- `calc_pwild()` is the simple calculation based on the proportion of
  hatchery spawners

- `calc_pwild_age()` performs the calculation weighted by age class
  fecundity

## Usage

``` r
calc_pwild(pHOS_cur, pHOS_prev, gamma)

calc_pwild_age(NOS_a, HOS_a, fec, gamma)
```

## Arguments

- pHOS_cur:

  Numeric, proportion of hatchery spawners in current generation

- pHOS_prev:

  Numeric, proportion of hatchery spawners in previous generation

- gamma:

  Numeric, reduced reproductive success of hatchery spawners

- NOS_a:

  Array `[nsim, maxage, years]` for natural spawners

- HOS_a:

  Array `[nsim, maxage, years]` for hatchery spawners

- fec:

  Array `[nsim, maxage, years]` for age class fecundity

## Value

`calc_pwild_age()` a matrix of pWILD by simulation and year.
`calc_pwild()` returns a numeric
