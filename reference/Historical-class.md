# Class `"Historical"`

Optional component of the operating model that specifies the historical
dynamics.

## Details

Several approaches are possible:

- No set up. Default option sets 1000 natural-origin juveniles (age 1),
  and 1000 hatchery-origin juveniles (age 1) if there is hatchery
  production (otherwise, zero).

- *Recommended option*: specify the initial spawning abundance in the
  terminal age class.

- Detailed setup that reconstructs a historical population by specifying
  the juvenile abundance (at the beginning of the year), annual fishing
  mortality rates, and spawner abundance. Typically used if there an
  estimation/conditioning model is used to inform parameters of the
  operating model.

## Slots

- `Name`:

  Character. Identifying name

- `HistSpawner_NOS`:

  Natural origin spawners at age. Either a numeric to specify the total
  natural spawners (in the oldest age class) at the beginning of the
  projection, otherwise, an array by `[nsim, maxage, nyears, n_g]`.
  Default is 1,000 spawners.

- `HistSpawner_HOS`:

  Hatchery origin spawners at age. Either a numeric to specify the total
  hatchery spawners (in the oldest age class) at the beginning of the
  projection, otherwise, an array by `[nsim, maxage, nyears, n_r]`.
  Default is 1,000 spawners if there is hatchery production or zero
  otherwise.

- `HistNjuv_NOS`:

  Array by `[nsim, maxage, nyears+1, n_g]`. The abundance of immature
  natural origin fish at the beginning of the annual time step. Default
  assumes 1000 smolts (age-1) fish annually.

- `HistNjuv_HOS`:

  Array by `[nsim, maxage, nyears+1, n_r]`. The abundance of immature
  hatchery origin fish at the beginning of the annual time step. Default
  assumes 1000 smolts (age-1) fish annually.

- `HistFPT`:

  Vector by historical years (`nyears`) or an array by dimension
  `[nsim, nyears, 2]`. The instantaneous fishing mortality in the
  preterminal fishery. The first array slice corresponds to F for
  natural origin fish and the second array slice corresponds to hatchery
  origin fish. Default is zero.

- `HistFT`:

  Vector by historical years (`nyears`) or an array by dimension
  `[nsim, nyears, 2]`. The instantaneous fishing mortality in the
  terminal fishery. The first array slice corresponds to F for natural
  origin fish and the second array slice corresponds to hatchery origin
  fish. Default is zero.

## Creating Object

Objects can be created by calls of the form `new("Historical")`

## Examples

``` r
showClass("Historical")
#> Class "Historical" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                       
#> Name:             Name HistSpawner_NOS HistSpawner_HOS    HistNjuv_NOS
#> Class:       character       num.array       num.array           array
#>                                                       
#> Name:     HistNjuv_HOS         HistFPT          HistFT
#> Class:           array       num.array       num.array
#> 
#> Extends: "Historical.list"
```
