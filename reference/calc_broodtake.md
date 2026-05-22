# Solve for broodtake numbers

Internal functions to calculate broodtake with various constraints,
called by
[`brood_func()`](https://docs.salmonmse.com/reference/brood_func.md).
`calc_broodtake()` is a wrapper function.

`calc_broodtake_custom()` uses a user-provided function to generate
brood numbers and adjusts downwards if egg production exceeds the
target.

`.broodtake_func()` is the optimization function used to calculate
broodtake

## Usage

``` r
calc_broodtake(
  NO,
  HO,
  stray,
  brood_import,
  ptarget_NOB,
  pmax_NOB,
  phatchery,
  egg_target,
  p_female,
  fec,
  s_prespawn,
  m
)

calc_broodtake_custom(
  f_brood,
  NO,
  HO,
  stray,
  p_female,
  fec,
  s_prespawn,
  m,
  egg_target
)

.broodtake_func(
  ptake_unmarked,
  NO,
  HO,
  stray,
  brood_import,
  phatchery,
  p_female,
  fec,
  egg_target,
  s_prespawn,
  ptarget_NOB,
  m = 1,
  opt = TRUE
)
```

## Arguments

- NO:

  Matrix `[nage, n_g]`, natural-origin fish available for broodtake

- HO:

  Matrix `[nage, n_r]`, hatchery-origin fish, potentially available for
  broodtake

- stray:

  Matrix `[nage, n_r]`, hatchery-origin strays available for broodtake

- brood_import:

  Vector `[nage]` of imported brood

- ptarget_NOB:

  Numeric, target proportion of `NOB/(NOB + HOB)`. If `m = 1`, then the
  realized pNOB should be ptarget_NOB. If `m < 1`, then the system
  achieves `unmarked brood/(NOB + HOB) = ptarget_NOB`.

- pmax_NOB:

  Numeric, maximum proportion of `NOB/NO`

- phatchery:

  Numeric, proportion of `HO` that return to hatchery instead of the
  spawning ground. HOB is taken from this subset. Set to `NA` to obtain
  HOB on the way to spawning ground.

- egg_target:

  Numeric, target egg production from which to back-calculate brood
  numbers

- p_female:

  Vector `[nage]` of proportion female to calculate egg production

- fec:

  Vector `[nage]`, egg production per female

- s_prespawn:

  Numeric, survival of brood prior to egg production

- m:

  Numeric, mark rate of `HO` used. Discounts `ptarget_NOB` and
  `pmax_NOB` based on marked fish

- f_brood:

  Function that calculates the brood numbers

- ptake_unmarked:

  Numeric, proportion of unmarked fish used for brood

- opt:

  Logical, whether the function is used for optimization (TRUE) or
  reporting (FALSE)

## Value

`calc_broodtake()` returns a list from `.broodtake_func()`

`calc_broodtake_custom()` returns a named list, same format as
`calc_broodtake()`

`.broodtake_func()` returns a numeric if `opt = TRUE`:
`log(p_unmarked) - log(ptarget_NOB)`, otherwise a named list of brood
and egg production.

- `egg_NOB` Matrix `[nage, n_g]`

- `egg_HOB_unmarked` Matrix `[nage, n_r]` (including strays)

- `egg_HOB_marked` Matrix `[nage, n_r]`

- `egg_HOB_import` Vector `[nage]`

- `ptake_unmarked` Numeric

- `ptake_marked` Numeric

- `pNOB` Numeric

- `NOB` Numeric

- `HOB_unmarked` Matrix `[nage, n_r]`

- `HOB_marked` Matrix `[nage, n_r]`

- `HOB_import` Matrix `[nage, n_r]`

- `HOB_stray` Vector `[nage]`
