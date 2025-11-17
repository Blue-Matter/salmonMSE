# Reference points

Calculate MSY and Sgen reference points for the operating model. Uses
the biological parameters (maturity, natural mortality) in the last year
of the projection.

- `calc_MSY()` calculates the MSY reference points from a set of
  biological and fishery parameters

- `calc_Sgen()` calculates the Sgen, the spawner abundance that would
  reach the spawner abundance at MSY after one generation without
  fishing

- `calc_ref()` is a wrapper function that calculates MSY and Sgen for an
  operating model

## Usage

``` r
calc_ref(SOM, rel_F, check = TRUE, maximize = c("MSY", "MER"))

calc_MSY(
  Mjuv,
  fec,
  p_female,
  rel_F,
  vulPT,
  vulT,
  p_mature,
  s_enroute,
  n_g = 1,
  p_LHG = 1,
  SRRpars,
  maximize = c("MSY", "MER"),
  F_search = c(1e-08, 5)
)

calc_Sgen(
  Mjuv,
  fec,
  p_female,
  rel_F,
  vulPT,
  vulT,
  p_mature,
  s_enroute,
  n_g = 1,
  p_LHG = 1,
  SRRpars,
  SMSY,
  F_search = c(1e-08, 100),
  nyears
)
```

## Arguments

- SOM:

  An object of class
  [SOM](https://docs.salmonmse.com/reference/SOM-class.md)

- rel_F:

  Numeric length 2, indicates the relative effort in the preterminal and
  terminal fisheries, with a maximum value of 1. The default is
  `c(0, 1)` which indicates a yield calculation with only the terminal
  fishery.

- check:

  Logical, whether to check the SOM object using
  [`check_SOM()`](https://docs.salmonmse.com/reference/check_SOM.md)

- maximize:

  Character, whether the MSY calculation is the optimum that maximizes
  catch (`"MSY"`) or excess recruitment (`"MER"`). The two methods
  should be equivalent when `rel_F = c(0, 1)`.

- Mjuv:

  Numeric `maxage` for juvenile natural mortality. Can be a matrix
  `[maxage, n_g]`.

- fec:

  Numeric `maxage` for fecundity. Can be a matrix `[maxage, n_g]`.

- p_female:

  Numeric for proportion female spawners

- vulPT:

  Numeric `maxage` for preterminal vulnerability at age

- vulT:

  Numeric `maxage` for terminal vulnerability at age

- p_mature:

  Numeric `maxage` for maturity proportions at age. Can be a matrix
  `[maxage, n_g]`.

- s_enroute:

  Numeric for en-route survival of escapement to spawning grounds

- n_g:

  Integer, number of life history groups within a cohort

- p_LHG:

  Numeric `n_g` for proportion of the total egg production assigned to
  each life history group within a cohort

- SRRpars:

  Data frame, one row, that contains the stock recruit parameters that
  predicts density-dependent survival at the egg-smolt life stage

- F_search:

  Numeric, length 2 for the range of F values to search for the
  instantaneous fishing mortality that produces MSY

- SMSY:

  Numeric, spawning abundance at MSY

- nyears:

  Integer, number of years to project the population with no fishing to
  reach `SMSY`. Default is the minimum age of maturity.

## Value

- `calc_MSY` returns a vector of various state variables (catch,
  exploitation rate, egg production, spawners) at MSY

- `calc_Sgen` returns a numeric

- `calc_ref` returns a list by stock, each containing a matrix of MSY
  state variables and Sgen by simulation

## See also

[`calc_Smsy_Ricker()`](https://docs.salmonmse.com/reference/calc_Smsy_Ricker.md)
