# Reference points for conditioning model

Internal functions that calculate productivity (`.CM_prod()`) and
reference points (`.CM_MSY()`) from the conditioning model. These
functions can be used to calculate values for a subset of years
(productivity can vary in time with natural mortality and maturity), or
from average biological parameters in a subset of years.

## Usage

``` r
.CM_prod(report, d, index = NULL, mean_bio = FALSE)

.CM_MSY(
  report,
  d,
  simple = TRUE,
  index = NULL,
  mean_bio = FALSE,
  type = c("spawner", "egg", "u", "Sgen"),
  AEQ = TRUE,
  maximize = c("MSY", "MER"),
  ncores = 1
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

- index:

  Integer vector to subset years with which to calculate reference
  points. Can be used to reduce computation or average biological
  parameters from a subset of years, see `mean_bio` argument. If `NULL`,
  uses all years of model.

- mean_bio:

  Logical, whether to average the natural mortality and maturity
  parameters across years indicated in `index`

- simple:

  Logical, whether to use Ricker lambert equations for MSY reference
  points (TRUE) or age-structured optimization (FALSE)

- type:

  Character, the type of reference point to calculate

- AEQ:

  Logical, whether to use adult equivalents when calculating MSY for
  preterminal fisheries#' Only used if `simple = FALSE`. Should always
  be `TRUE`

- maximize:

  Character, whether the numerical optimization maximizes catch (MSY) or
  excess recruitment (MER). For testing only, should not impact results.
  Only used if `simple = FALSE`

- ncores:

  Numeric, number of processors for parallel computation. Useful if
  calculating numerically for many MCMC samples

## Value

Matrix, dimension `[length(index), length(report)]`. If
`mean_bio = TRUE`, matrix has 1 row.

## See also

[`CM_MSY()`](https://docs.salmonmse.com/reference/CM_MSY.md)
