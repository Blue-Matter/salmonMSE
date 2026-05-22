# Brood function

Calculates broodtake and in-river removals from escapement of marine
fisheries. This function also applies en-route mortality.

## Usage

``` r
brood_func(
  Nage_NOS,
  Nage_HOS,
  stray_external,
  m,
  m_stray,
  s_enroute,
  hatchery_args
)
```

## Arguments

- Nage_NOS:

  Array `[nage, n_g]` of natural-origin fish

- Nage_HOS:

  Array `[nage, n_r]` of hatchery-origin fish

- stray_external:

  Array `[nage, n_r]` of hatchery-origin strays

- m:

  Numeric, mark rate of `Nage_HOS`

- m_stray:

  Numeric, mark rate of `stray_external`

- s_enroute:

  Numeric, en-route survival

- hatchery_args:

  List of various hatchery arguments created by
  [`define_hatchery_args()`](https://docs.salmonmse.com/reference/define_hatchery_args.md)
  and adjusted by
  [`ProjectSOM()`](https://docs.salmonmse.com/reference/salmonMSE.md).

## Value

Named list:

- `broodtake` list returned by
  [`calc_broodtake()`](https://docs.salmonmse.com/reference/calc_broodtake.md)

- `hatchery_production` list returned by
  [`calc_yearling()`](https://docs.salmonmse.com/reference/calc_yearling.md)

- `spawners` list returned by
  [`calc_spawners()`](https://docs.salmonmse.com/reference/calc_spawners.md)

## See also

[`calc_broodtake()`](https://docs.salmonmse.com/reference/calc_broodtake.md)
