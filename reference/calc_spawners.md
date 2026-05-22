# Calculate spawners after broodtake

Internal function that calculates remaining spawners after broodtake and
additional removals.

## Usage

``` r
calc_spawners(broodtake, NO, HO, stray, phatchery, premove_HOS, premove_NOS, m)
```

## Arguments

- broodtake:

  List, returned by
  [`calc_broodtake()`](https://docs.salmonmse.com/reference/calc_broodtake.md)

- NO:

  Matrix `[nage, n_g]`, natural-origin in-river return

- HO:

  Matrix `[nage, n_r]`, hatchery-origin in-river return

- stray:

  Matrix `[nage, n_r]`, hatchery-origin strays in-river return

- phatchery:

  Numeric, proportion of `HO` that return to hatchery instead of the
  spawning ground. HOB is taken from this subset. Set to `NA` to obtain
  HOB on the way to spawning ground.

- premove_HOS:

  Numeric, proportion of HO fish to remove (after brood). Can also be a
  function

- premove_NOS:

  Numeric, proportion of NO fish to remove (after brood). Can also be a
  function

- m:

  Numeric, mark rate

## Value

Named list:

- NOS `[nage, n_g]` natural-origin spawners

- HOS `[nage, n_r]` hatchery-origin spawners

- HOS_stray `[nage, n_r]` spawners that are strays

- NO_remove `[nage, n_g]` natural-origin fish removed before spawning

- HO_remove `[nage, n_r]` hatchery-origin fish removed before spawning
