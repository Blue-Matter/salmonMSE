# Spawning and early life stage function

Calculates egg production from spawners arriving at spawning grounds,
broodtake and in-river removals from escapement of marine fisheries.
This function also applies en-route mortality.

## Usage

``` r
nextgen_SRR_func(
  Brood_Calcs,
  fec,
  p_female,
  hatchery_args = list(),
  fitness_args = list(),
  zbar_brood,
  SRRpars,
  p_LHG
)

nextgen_habitat_func(
  Brood_Calcs,
  Habitat,
  fec,
  p_female,
  hatchery_args = list(),
  fitness_args = list(),
  zbar_brood,
  p_LHG
)
```

## Arguments

- Brood_Calcs:

  List, returned by
  [`brood_func()`](https://docs.salmonmse.com/reference/brood_func.md)

- fec:

  Vector `[nage]`, fecundity at age of female fish

- p_female:

  Vector `[nage]`, proportion female by age class

- hatchery_args:

  List of arguments created by
  [`define_hatchery_args()`](https://docs.salmonmse.com/reference/define_hatchery_args.md)

- fitness_args:

  List of arguments created by
  [`define_fitness_args()`](https://docs.salmonmse.com/reference/define_hatchery_args.md)

- Habitat:

  [Habitat](https://docs.salmonmse.com/reference/Habitat-class.md)
  object, modified by
  [`ProjectSOM()`](https://docs.salmonmse.com/reference/salmonMSE.md)

## Value

Named list:

- `NOS` Matrix `[nage, n_g]`

- `HOS` Matrix `[nage, n_r]`

- `HOS_effective` Matrix `[nage, n_r]`

- `HOS_stray` Matrix `[nage, n_r]`

- `pHOSeff` Numeric

- `pHOScensus` Numeric

- `Egg_NOS` `[nage, n_g]`

- `Egg_HOS` `[nage, n_r]`

- `Smolt_RelOut` Vector `[n_r]`

- `Fry_NOS` Vector `[n_g]`

- `Fry_HOS` Vector `[n_g]`

- `Smolt_NOS` Vector `[n_g]`

- `Smolt_HOS` Vector `[n_g]`
