# Tradeoff figure

Generates a tradeoff figure, a comparison between two performance
metrics, across two variables which may represent a population dynamics
variable (e.g., productivity) or a management action (e.g., hatchery
production levels or harvest strategy). See example at
<https://docs.salmonmse.com/articles/decision-table.html>

## Usage

``` r
plot_tradeoff(
  pm1,
  pm2,
  x1,
  x2,
  xlab,
  ylab,
  x1lab,
  x2lab,
  scenario,
  ncol = NULL
)
```

## Arguments

- pm1:

  Numeric or matrix. A vector of values for the first performance metric
  on the x-axis. Alternatively, provide a three column matrix
  corresponding to the lower bound, central tendency, and upper bound.

- pm2:

  Numeric or matrix. A vector of values for the second performance
  metric on the y-axis (same length as pm1). Alternatively, provide a
  three column matrix corresponding to the lower bound, central
  tendency, and upper bound.

- x1:

  Atomic, vector of values for the first grouping variable. Various
  levels are represented by colours. Same length as pm1.

- x2:

  Numeric, vector of values for the second grouping variable. Various
  levels are represented by shapes. Same length as pm1.

- xlab:

  Character, optional x-axis label

- ylab:

  Character, optional y-axis label

- x1lab:

  Character, optional label for the first grouping variable

- x2lab:

  Character, optional label for the second grouping variable

- scenario:

  Atomic, vector of faceting variables (same length as `pm1`, `pm2`)
  used to generate a grid of decision tables

- ncol:

  Integer, number of columns for decision table grid, only used if
  `scenario is provided`

## Value

ggplot object

## See also

[`plot_statevar_ts()`](https://docs.salmonmse.com/reference/plot_statevar_ts.md)
[`plot_decision_table()`](https://docs.salmonmse.com/reference/plot_decision_table.md)
