# Compare simulation runs

Create figures that compare results across two dimensions

## Usage

``` r
compare_spawners(SMSE_list, Design, prop = FALSE, FUN = median)

compare_fitness(SMSE_list, Design, FUN = median)

compare_escapement(SMSE_list, Design, FUN = median)
```

## Arguments

- SMSE_list:

  A list of SMSE objects returned by
  [`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md)

- Design:

  A data frame with two columns that describes the factorial design of
  the simulations. Used to label the figure. Rows correspond to each
  object in `SMSE_list`. There two columns are variables against which
  to plot the result. See example in
  <https://docs.salmonmse.com/articles/decision-table.html>.

- prop:

  Logical, whether to plot absolute numbers over proportions

- FUN:

  Summarizing function across simulations, typically
  [`stats::median()`](https://rdrr.io/r/stats/median.html) or
  [`base::mean()`](https://rdrr.io/r/base/mean.html)

## Value

A ggplot object

## Details

- `compare_spawners()` generates a time series of the composition of
  spawners

- `compare_fitness()` generates a time series of metrics (fitness, PNI,
  pHOS, and pWILD) related to hatchery production

- `compare_escapement()` generates a time series of the proportion of
  spawners and broodtake to escapement
