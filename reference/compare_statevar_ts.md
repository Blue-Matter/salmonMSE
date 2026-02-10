# Compare state variables from simulation runs

Compare outputs from multiple simulations to evaluate performance across
states of nature and/or management levers (identified by colour):

- `compare_statevar_ts()` produces a time series for all simulations, or
  with medians and 95th percentile intervals

- `compare_statevar_hist()` produces a histogram or density plot across
  all simulations for a particular year

## Usage

``` r
compare_statevar_ts(
  SMSE_list,
  var = "PNI",
  s = 1,
  figure = TRUE,
  xlab = "Projection Year",
  quant = FALSE,
  ylab = var,
  ylim,
  agg.fun = sum,
  names,
  col_vec,
  ...
)

compare_statevar_hist(
  SMSE_list,
  var = "PNI",
  s = 1,
  y,
  figure = TRUE,
  xlab = var,
  names,
  col_vec,
  type = c("density", "hist"),
  ...
)
```

## Arguments

- SMSE_list:

  List of SMSE objects for multiple model runs returned by
  [`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md)

- var:

  Character. Slot for the state variable in `SMSE` object. See
  `slotNames(SMSE)` for options. Additional supported options are:
  `"ESS"` (egg-smolt survival), `"pbrood"` (broodtake to escapement
  ratio), `"pNOSesc"` (NOS/natural escapement), `"pHOSesc"`
  (HOS/hatchery escapement), `Total Spawners` (NOS + HOS), `NOS/SMSY`,
  `S/SMSY`, and `NOS/Sgen`.

- s:

  Integer. Population index for multi-population model (e.g., `s = 1` is
  the first population in the model)

- figure:

  Logical, whether to generate a figure (set to FALSE if only using the
  function to return the data matrix)

- xlab:

  Character. Name of time variable for the figure

- quant:

  Logical, whether to plot individual simulations (FALSE) or the median
  with 95 percent confidence intervals (TRUE)

- ylab:

  Character. Name of the state variable for the figure

- ylim:

  Vector. Y-axis limits

- agg.fun:

  Function. Defines how to aggregate state variables that are reported
  by age. Typically, `sum` is used but `max` is also possible for
  reporting apical exploitation rates.

- names:

  Character vector `length(SMSE_list)` to label individual model runs

- col_vec:

  Character vector `length(SMSE_list)` for custom colour schemes for
  comparing across model scenarios in figures

- ...:

  Additional arguments to base plot function

- y:

  Integer. Projection year for the state variable to plot the histogram.
  If missing, the last projection year is used.

- type:

  Character, whether to generate a density figure or histogram

## Value

An array invisibly. Also generates base graphics if `figure = TRUE`

## See also

[`plot_statevar_ts()`](https://docs.salmonmse.com/reference/plot_statevar_ts.md)
[`compare()`](https://docs.salmonmse.com/reference/compare.md)
