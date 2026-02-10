# Class `"Historical"`

Optional component of the operating model that specifies the historical
dynamics. The default option starts the projection with 1000
natural-origin juveniles and 1000 hatchery-origin juveniles (if there is
hatchery production) in the oldest age class, representing single-brood
year returns since there is only one age class in the population.
Specify the abundance in all age class to simulate multiple brood-year
returns.

## Slots

- `Name`:

  Character. Identifying name

- `InitNjuv_NOS`:

  Array by `[nsim, maxage, n_g]`. The abundance of immature natural
  origin fish at the beginning of the projection. Default assumes 1000
  in the oldest age class, which creates a population with single brood
  year returns.

- `InitNjuv_HOS`:

  Array by `[nsim, maxage, n_r]`. The abundance of immature hatchery
  origin fish at the beginning of the projection. Default assumes 1000
  in the oldest age class, which creates a population with single brood
  year returns.

## Creating Object

Objects can be created by calls of the form `new("Historical")`

## Examples

``` r
showClass("Historical")
#> Class "Historical" [package "salmonMSE"]
#> 
#> Slots:
#>                                              
#> Name:          Name InitNjuv_NOS InitNjuv_HOS
#> Class:    character    num.array    num.array
#> 
#> Extends: "Historical.list"
```
