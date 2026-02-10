# Decision table of performance metrics

Generates a coloured table of a performance metric across two axes,
which may be a population dynamics variable (e.g., productivity) or a
management action (e.g., hatchery production levels or harvest
strategy). See example at
<https://docs.salmonmse.com/articles/decision-table.html>

## Usage

``` r
plot_decision_table(
  x,
  y,
  z,
  title,
  xlab,
  ylab,
  scenario,
  ncol = NULL,
  dir = "v"
)
```

## Arguments

- x:

  Atomic, vector of values for the x axis (same length as z). Will be
  converted to factors

- y:

  Atomic, vector of values for the y axis (same length as z). Will be
  converted to factors

- z:

  Numeric, vector of values for the performance metric

- title:

  Character, optional title of figure

- xlab:

  Character, optional x-axis label

- ylab:

  Character, optional y-axis label

- scenario:

  Atomic, vector of faceting variables (same length as z) used to
  generate a grid of decision tables

- ncol:

  Integer, number of columns for decision table grid, only used if
  `scenario` is provided

- dir:

  Character, either "h" or "v" to describe how the grid of tables should
  be organized (horizontally or vertically)

## Value

ggplot object

## See also

[`plot_statevar_ts()`](https://docs.salmonmse.com/reference/plot_statevar_ts.md)
[`plot_tradeoff()`](https://docs.salmonmse.com/reference/plot_tradeoff.md)
