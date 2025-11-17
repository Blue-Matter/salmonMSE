# Plot life history groups and release strategies

Plot the annual proportions of life history groups (natural origin fish)
or release strategies (hatchery origin) at various life stages

## Usage

``` r
plot_LHG(
  SMSE,
  var = "NOS",
  type = c("prop", "abs"),
  s = 1,
  FUN = median,
  figure = TRUE,
  xlab = "Projection Year",
  ylab,
  name,
  ylim
)

plot_RS(
  SMSE,
  var = "HOS",
  type = c("prop", "abs"),
  s = 1,
  FUN = median,
  figure = TRUE,
  xlab = "Projection Year",
  ylab,
  name,
  ylim
)
```

## Arguments

- SMSE:

  Class [SMSE](https://docs.salmonmse.com/reference/SMSE-class.md)
  object returned by
  [`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md)

- var:

  Character. Slot for the state variables in `SMSE@Misc$LHG[[1]]` or
  `SMSE@Misc$RS[[1]]`.

- type:

  Character to indicate whether to plot proportion or absolute numbers

- s:

  Integer. Population index for multi-population model (e.g., `s = 1` is
  the first population in the model)

- FUN:

  Summarizing function across simulations, typically
  [`median()`](https://rdrr.io/r/stats/median.html) or
  [`mean()`](https://rdrr.io/r/base/mean.html)

- figure:

  Logical, whether to generate a figure (set to FALSE if only using the
  function to return the data matrix)

- xlab:

  Character. Name of time variable for the figure

- ylab:

  Character. Name of the state variable for the figure

- name:

  Character. Vector of names for the life history groups or release
  strategies

- ylim:

  Vector length 2, y-axis limits

## Value

Base graphics figure, barplot of distribution or total numbers by LHG or
RS. Returns invisibly the matrix of plotted values

## See also

[`plot_statevar_ts()`](https://docs.salmonmse.com/reference/plot_statevar_ts.md)
