# Update natural mortality of juveniles

Internal function that updates juvenile natural mortality in the marine
environement due to fitness

- `SAR_fitness()` calculates the new natural mortality value

- `makeRel_SAR` generates a list for openMSE to use in the simulations

## Usage

``` r
SAR_fitness(
  x = -1,
  y = 1,
  envir = c("natural", "hatchery"),
  rel_loss = 1,
  s = 1,
  nyears,
  Mbase
)

makeRel_SAR(
  p_smolt = 1,
  s = 1,
  envir = c("natural", "hatchery"),
  rel_loss,
  nyears,
  Mbase
)
```

## Arguments

- x:

  Integer, simulation number from openMSE

- y:

  Integer, simulation year (including historical years)

- envir:

  Character, whether to obtain the fitness value for the natural or
  hatchery environment.

- rel_loss:

  Numeric, the loss exponent for the juveniles

- s:

  Integer, the salmonMSE population index. Used to search for the
  fitness value

- nyears:

  Integer, the number of historical years in the operating model

- Mbase:

  Array `[nsim, n_age, proyears]` the base natural mortality value in
  the openMSE operating model.

- p_smolt:

  Integer, the population index for the juvenile population in the
  openMSE model

## Value

- [`smolt_func()`](https://docs.salmonmse.com/reference/smolt_func.md)
  returns a numeric for the ratio of the realized smolt production vs.
  the hypothetical value if there were no hatchery, en route mortality,
  or habitat improvement

- `makeRel_smolt` returns a list that is passed to openMSE as a
  inter-population relationship
