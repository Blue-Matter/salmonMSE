# Class `"simpleSOM"`

A simple operating model for modeling recruit-spawners and a terminal
marine fishery, no hatchery production.

## Details

Various parameters can be stochastic (length `nsim`) or input as a
single numeric (value identical across all simulations).

## Slots

- `Name`:

  Character. Identifying name

- `nsim`:

  Integer. Number of simulations.

- `maxage`:

  Integer. The maximum age of the population age structure.

- `ngen`:

  Integer. The number of generations (life cycles) to run the
  projection.

- `kappa`:

  Vector length `nsim`. The adult productivity ratio for the
  stock-recruit function. **Units of recruits per spawner.** Natural
  per-capita production of recruits as the population approaches zero
  (density-independent component).

- `Smax`:

  Vector length `nsim`. The spawner abundance that maximizes smolt
  production in the Ricker stock-recruit function. **Units of
  spawners.**

- `sigmaR`:

  Vector length `nsim`. The lognormal standard deviation in the Ricker
  stock-recruit relationship. Recruitment anomalies will be simulated
  with `[stats::rlnorm()]`.

- `Recdev`:

  Optional matrix `[nsim, ngen-1]` of recruitment anomalies in the
  Ricker stock-recruit relationship. Overrides the `sigmaR` argument and
  provides more user flexibility, for example, simulating with
  auto-correlation.

- `type_T`:

  Character. Whether to manage terminal fishery catch from exploitation
  rate ("u") or catch target ("catch"). Default is "u".

- `u_terminal`:

  Numeric, matrix `[nsim, proyears]`, or function. If `type_T = "u"`,
  the harvest rate (ratio of kept catch to of the terminal marine
  fishery. Function should be of the form
  `function(NO, HO, m) return(u)`.

- `K_T`:

  Numeric or function. If `type_T = "catch"`, the catch target of the
  return in the terminal fishery. Function should be of the form
  `function(NO, HO, m) return(K)`.

- `InitReturn`:

  Single numeric or vector `[nsim]`. The return at the beginning of the
  projection. Default assumes 1000.

## Creating Object

Objects can be created by calls of the form `new("simpleSOM")`

## Examples

``` r
showClass("simpleSOM")
#> Class "simpleSOM" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                   
#> Name:                 Name                nsim              maxage
#> Class:           character             numeric             numeric
#>                                                                   
#> Name:                 ngen               kappa                Smax
#> Class:             numeric             numeric             numeric
#>                                                                   
#> Name:               sigmaR              Recdev              type_T
#> Class:             numeric              matrix           character
#>                                                                   
#> Name:           u_terminal                 K_T          InitReturn
#> Class: num.matrix.function        num.function             numeric
```
