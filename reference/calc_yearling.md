# Calculate hatchery releases

From egg target, calculate production of yearlings and subyearlings
conditional on survival and proportion of releases.

`.yearling_func()` is the optimization function that determines
proportions.

## Usage

``` r
calc_yearling(egg_target, s_yearling, s_subyearling, p_yearling, p_subyearling)

.yearling_func(
  p_egg_yearling,
  egg_target,
  s_yearling,
  s_subyearling,
  p_yearling,
  opt = TRUE
)
```

## Arguments

- egg_target:

  Numeric

- s_yearling:

  Numeric, survival from egg to yearling

- s_subyearling:

  Numeric, survival from egg to subyearling

- p_yearling:

  Vector `[n_r]`, proportion of releases at subyearling stage, where
  `sum(p_yearling, p_subyearling) = 1`

- p_egg_yearling:

  Numeric, proportion to eggs to become yearlings

- opt:

  Logical, whether the function is used for optimization (TRUE) or
  reporting (FALSE)

## Value

`calc_yearling()` returns a named list returned by `.yearling_func()`

`.yearling_func()` returns a numeric if `opt = TRUE`:
`yearling/(yearling + subyearling) - p_yearling`. Otherwise,
log(p_unmarked) - log(ptarget_NOB)\`, otherwise a named list:

- `yearling` Vector `[n_r]` of yearling releases

- `subyearling` Vector `[n_r]` of subyearling releases
