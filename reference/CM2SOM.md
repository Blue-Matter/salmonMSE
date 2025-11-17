# Convert conditioning model to operating model

Creates an operating model from MCMC samples and data inputs of the
conditioning model. Management actions for habitat, hatchery production,
and harvest still need to be specified in the operating model.

**Note: the function assumes the maturity values in the last
conditioning year for the projection, which are likely not well informed
by CWT data. Consider updating the maturity using some historical
average (e.g., across most recent completed brood years).**

## Usage

``` r
CM2SOM(stanfit, sims, nsim = 2, seed = 1, proyears = 40)
```

## Arguments

- stanfit:

  Output from
  [`sample_CM()`](https://docs.salmonmse.com/reference/fit_CM.md)

- sims:

  Optional, a vector of integers indicating the MCMC iterations to
  convert to operating model simulations. Otherwise, use argument `nsim`
  in order to sample a subset of the MCMC.

- nsim:

  Integer, total number of simulations in the operating model. Only used
  if `sims` is missing.

- seed:

  Integer, seed for sampling the MCMC output. Only used if `sims` is
  missing.

- proyears:

  Integer, the number of projection years in the operating model

## Value

[SOM](https://docs.salmonmse.com/reference/SOM-class.md) object.
