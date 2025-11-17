# Class `"Harvest"`

The component of the operating model that controls harvest.

## Slots

- `Name`:

  Character. Identifying name

- `type_PT`:

  Character. Whether to manage preterminal fishery catch from
  exploitation rate ("u") or catch target ("catch"). Default is "u",

- `type_T`:

  Character. Whether to manage terminal fishery catch from exploitation
  rate ("u") or catch target ("catch"). Default is "u",

- `u_preterminal`:

  Numeric. If `type_PT = "u"`, the exploitation rate of the immature
  stock in the pre-terminal fishery. This will be converted to an
  instantaneous fishing mortality rate, i.e.,
  `F_preterminal = -log(1 - u_preterminal)`.

- `u_terminal`:

  Numeric. If `type_T = "u"`, The exploitation rate of the return in the
  terminal fishery. This will be converted to an instantaneous fishing
  mortality rate, i.e., `F_terminal = -log(1 - u_terminal)`.

- `K_PT`:

  Numeric. If `type_PT = "catch"`, the catch target of the immature
  stock in the pre-terminal fishery.

- `K_T`:

  Numeric. If `type_T = "catch"`, the catch target of the return in the
  terminal fishery.

- `MSF_PT`:

  Logical. Whether to implement mark-selective fishing in the
  preterminal fishery, with no retention on unmarked fish.

- `MSF_T`:

  Logical. Whether to implement mark-selective fishing in the terminal
  fishery, with no retention on unmarked fish.

- `release_mort`:

  Vector length 2. The proportion of released fish that die after
  release, in the pre-terminal and terminal fishery. Implemented to
  model mark-selective fishing. Not used if either `MSF_PT` or `MSF_T`
  is ` FALSE`.

- `vulPT`:

  Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability
  schedule (between 0-1) in the preterminal fishery. Values indicate the
  proportion of fishing mortality experienced by each age class, where
  `F_preterminal = -log(1 - u_preterminal)`.

- `vulT`:

  Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability
  schedule (between 0-1) in the terminal fishery. Values indicate the
  proportion of fishing mortality experienced by each age class, where
  `F_terminal = -log(1 - u_terminal)`.

## Creating Object

Objects can be created by calls of the form `new("Harvest")`

## Examples

``` r
showClass("Harvest")
#> Class "Harvest" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                             
#> Name:           Name       type_PT        type_T u_preterminal    u_terminal
#> Class:     character     character     character       numeric       numeric
#>                                                                             
#> Name:           K_PT           K_T        MSF_PT         MSF_T  release_mort
#> Class:       numeric       numeric       logical       logical       numeric
#>                                   
#> Name:          vulPT          vulT
#> Class:    num.matrix    num.matrix
#> 
#> Extends: "Harvest.list"
```
