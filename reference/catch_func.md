# Catch function

Internal function that does catch calculations for the preterminal and
terminal marine fisheries.

## Usage

``` r
catch_func(
  NO,
  HO,
  type = c("u", "catch"),
  U,
  K,
  V,
  MSF = FALSE,
  m = 1,
  release_mort = 0
)
```

## Arguments

- NO:

  Array `[ns, nage, n_g]` of natural-origin fish

- HO:

  Array `[ns, nage, n_r]` of hatchery-origin fish

- type:

  Character. Whether to calculate removals from harvest rate `"u"`or
  catch `"catch"`

- U:

  Numeric or function. Harvest rate of fishery

- K:

  Numeric or function. Total catch of the fishery

- V:

  Matrix `[ns, nage]` Relative vulnerability by age class to fishery

- MSF:

  Logical, whether fishing is mark-selective

- m:

  Numeric vector length `[ns]`, mark rate of hatchery-origin fish. Only
  used if `MSF = TRUE`.

- release_mort:

  Numeric vector length `[ns]`, proportion of released unmarked fish
  that die. Only used if `MSF = TRUE`.

## Value

Named list:

- `K_NO` kept natural-origin catch, same dimension as `NO`

- `K_HO` kept hatchery-origin catch, same dimension as `HO`

- `D_NO` discarded (live + dead) natural-origin catch, same dimension as
  `NO`

- `D_HO` discarded (live + dead) hatchery-origin catch, same dimension
  as `HO`

- `DD_NO` dead discarded natural-origin catch, same dimension as `NO`

- `DD_HO` dead discarded hatchery-origin catch, same dimension as `HO`

- `U_NO` natural-origin harvest rate (ratio of kept catch and
  abundance), same dimension as `NO`

- `U_HO` hatchery-origin harvest rate (ratio of dead catch and
  abundance), same dimension as `HO`

- `Ex_NO` natural-origin exploitation rate (ratio of dead catch and
  abundance), same dimension as `NO`

- `Ex_HO` hatchery-origin exploitation rate (ratio of dead catch and
  abundance), same dimension as `HO`
