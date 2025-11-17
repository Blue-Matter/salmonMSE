# Plot core output from salmonMSE

Various functions that plot the state variables from salmonMSE
projections:

- `plot_statevar_ts()` produces a time series for all simulations, or
  with confidence intervals

- `plot_statevar_hist()` produces a histogram across all simulations for
  a particular year

- `plot_spawners()` produces a summary barplot of spawners, including
  NOS, HOS, and wild spawners

- `plot_escapement()` produces a summary figure of the proportion of
  spawners and broodtake to escapement

- `plot_fitness()` produces a summary figure of metrics (fitness, PNI,
  pHOS, and pWILD) related to hatchery production

- `plot_fishery()` produces a summary figure of metrics related to the
  fishery, e.g., median catch, exploitation rate or harvest rate

## Usage

``` r
plot_statevar_ts(
  SMSE,
  var = "PNI",
  s = 1,
  figure = TRUE,
  xlab = "Projection Year",
  quant = FALSE,
  ylab = var,
  ylim,
  agg.fun = sum,
  ...
)

plot_statevar_hist(SMSE, var = "PNI", s = 1, y, figure = TRUE, xlab = var, ...)

plot_spawners(SMSE, s = 1, prop = TRUE, FUN = median, figure = TRUE, ylim)

plot_fitness(SMSE, s = 1, FUN = median, figure = TRUE, ylim)

plot_escapement(SMSE, s = 1, FUN = median, figure = TRUE, ylim)

plot_fishery(
  SMSE,
  s = 1,
  type = c("catch", "exploit", "harvest"),
  FUN = median,
  figure = TRUE,
  ylim,
  ylab,
  ...
)

plot_Kobe(
  SMSE,
  s = 1,
  FUN = median,
  figure = TRUE,
  xlim,
  ylim,
  xlab = expression(NOS/S[MSY]),
  ylab = expression(U/U[MSY]),
  type = c("T", "PT")
)
```

## Arguments

- SMSE:

  Class [SMSE](https://docs.salmonmse.com/reference/SMSE-class.md)
  object returned by
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

- ...:

  Additional arguments to base plot function

- y:

  Integer. Projection year for the state variable to plot the histogram.
  If missing, the last projection year is used.

- prop:

  Logical, whether to plot proportions or absolute numbers

- FUN:

  Summarizing function across simulations, typically
  [`median()`](https://rdrr.io/r/stats/median.html) or
  [`mean()`](https://rdrr.io/r/base/mean.html)

- type:

  For `plot_Kobe`, the fishery state variable to plot. Whether to plot
  the exploitation rate for the terminal (T) or pre-terminal fishery
  (PT).

- xlim:

  Vector. X-axis limits

## Value

Functions return the matrix of plotted values invisibly. Figure plotted
from base graphics

## See also

[`plot_decision_table()`](https://docs.salmonmse.com/reference/plot_decision_table.md)
[`plot_LHG()`](https://docs.salmonmse.com/reference/plot_LHG.md)
