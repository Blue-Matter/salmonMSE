# Class `"Habitat"`

The component of the operating model that controls survival in the
freshwater environment. Includes changes in survival from either
environmental/climate effects or habitat mitigation.

## Slots

- `Name`:

  Character. Identifying name

- `use_habitat`:

  Logical. If `TRUE`, utilize stage-specific density-dependent functions
  from egg production from incubation mortality, egg-to-fry production,
  and fry-to-smolt production with annual deviations. Otherwise, the
  density-dependence is modeled for egg-to-smolt survival. See `Bio`
  object.

- `prespawn_rel`:

  Character, density-dependent function for pre-spawn mortality, e.g.,
  for spawners to reach spawning sites. Choices are "BH" (Beverton-Holt)
  or "HS" (hockey stick). Default is "BH".

- `prespawn_prod`:

  Numeric, productivity for pre-spawn mortality. Default is 1. Default
  if `Inf`.

- `prespawn_capacity`:

  Numeric, capacity for pre-spawn mortality. Default is `Inf`, i.e.,
  density-independence. Default is `Inf`.

- `egg_rel`:

  Character, density-dependent function for egg production from total
  spawning output. Choices are "BH" (Beverton-Holt) or "HS" (hockey
  stick). Default is "BH".

- `egg_prod`:

  Numeric, productivity for egg production from total spawning output
  (incubation). Default is 1. Default if `Inf`.

- `egg_capacity`:

  Numeric, capacity for egg production from total spawning output
  (incubation). Default is `Inf`, i.e., density-independence. Default is
  `Inf`.

- `fry_rel`:

  Character, density-dependent function for egg-to-fry production.
  Choices are "BH" (Beverton-Holt) or "HS" (hockey stick). Default is
  "BH".

- `fry_prod`:

  Numeric between 0-1, productivity for egg production from total
  spawning output, i.e., maximum survival as egg production approaches
  zero. Default is 0.4.

- `fry_capacity`:

  Numeric, capacity for fry production from egg production. Default is
  `Inf`, i.e., for density-independence. Default is `Inf`.

- `fry_sdev`:

  Matrix `[nsim, proyears]`, deviations from the density-dependent
  egg-fry survival. Can be utilized to incorporate time-varying
  environmental, climate, or habitat mitigation effects. Default is
  `matrix(1, nsim, proyears)`.

- `smolt_rel`:

  Character, density-dependent function for fry-to-smolt production.
  Choices are "BH" (Beverton-Holt) or "HS" (hockey stick). Default is
  "BH".

- `smolt_prod`:

  Numeric between 0-1, productivity for smolt production from fry, i.e.,
  maximum survival as fry production approaches zero. Default is 1.

- `smolt_capacity`:

  Numeric, capacity for smolt production from fry production. Set to
  `Inf` for density-independence. Default is `Inf`.

- `smolt_sdev`:

  Matrix `[nsim, proyears]`, deviations from the density-dependent
  fry-smolt survival. Can be utilized to incorporate time-varying
  environmental, climate, or habitat mitigation effects. Default is
  `matrix(1, nsim, proyears)`.

## Creating Object

Objects can be created by calls of the form `new("Habitat")`

## Examples

``` r
showClass("Habitat")
#> Class "Habitat" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                               
#> Name:               Name       use_habitat      prespawn_rel     prespawn_prod
#> Class:         character           logical         character           numeric
#>                                                                               
#> Name:  prespawn_capacity           egg_rel          egg_prod      egg_capacity
#> Class:           numeric         character           numeric           numeric
#>                                                                               
#> Name:            fry_rel          fry_prod      fry_capacity          fry_sdev
#> Class:         character           numeric           numeric            matrix
#>                                                                               
#> Name:          smolt_rel        smolt_prod    smolt_capacity        smolt_sdev
#> Class:         character           numeric           numeric            matrix
#> 
#> Extends: "Habitat.list"
```
