# Plot reference points from conditioning model

Plots time series of MSY reference points. Not used in automated
reporting as they can be computationally expensive and individual case
studies may require specific assumptions.

## Usage

``` r
CM_MSY(
  report,
  d,
  year1 = 1,
  simple = FALSE,
  index = NULL,
  mean_bio = FALSE,
  type = c("spawner", "egg", "u"),
  maximize = c("MSY", "MER"),
  AEQ = TRUE,
  ncores = 1,
  na.rm = FALSE
)

CM_Sgen(
  report,
  d,
  year1 = 1,
  simple = FALSE,
  index = NULL,
  mean_bio = FALSE,
  ncores = 1,
  na.rm = FALSE
)
```

## Arguments

- report:

  List, output of state variables from individual MCMC samples, obtained
  with
  [`get_report()`](https://docs.salmonmse.com/reference/CMfigures.md)

- d:

  List of data variables, obtained with
  [`get_CMdata()`](https://docs.salmonmse.com/reference/CMfigures.md)

- year1:

  Numeric, first year of model

- simple:

  Logical, whether to use Ricker lambert equations for MSY reference
  points (TRUE) or age-structured optimization (FALSE)

- index:

  Integer vector to subset years with which to calculate reference
  points. Can be used to reduce computation or average biological
  parameters from a subset of years, see `mean_bio` argument. If `NULL`,
  uses all years of model.

- mean_bio:

  Logical, whether to average the natural mortality and maturity
  parameters across years indicated in `index`

- type:

  Character, the type of reference point to calculate

- maximize:

  Character, whether the numerical optimization maximizes catch (MSY) or
  excess recruitment (MER). For testing only, should not impact results.
  Only used if `simple = FALSE`

- AEQ:

  Logical, whether to use adult equivalents when calculating MSY for
  preterminal fisheries#' Only used if `simple = FALSE`. Should always
  be `TRUE`

- ncores:

  Numeric, number of processors for parallel computation. Useful if
  calculating numerically for many MCMC samples

- na.rm:

  Logical, whether to exclude negative values from the median in figures

## Value

ggplot object

## See also

[`CM_prod()`](https://docs.salmonmse.com/reference/CMfigures.md)
[`CM_Srep()`](https://docs.salmonmse.com/reference/CMfigures.md)
[`.CM_MSY()`](https://docs.salmonmse.com/reference/dot-CM_prod.md)
