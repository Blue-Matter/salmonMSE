# Class `"Harvest"`

The component of the operating model that controls marine harvest.

## Slots

- `Name`:

  Character. Identifying name

- `type_PT`:

  Character. Whether to manage preterminal fishery catch from
  exploitation rate ("u") or catch target ("catch"). Default is "u".

- `type_T`:

  Character. Whether to manage terminal fishery catch from exploitation
  rate ("u") or catch target ("catch"). Default is "u".

- `u_preterminal`:

  Numeric, matrix `[nsim, proyears]`, or function. If `type_PT = "u"`,
  the harvest rate of the immature component of the population in the
  pre-terminal fishery. The harvest rate is the ratio to kept AEQ catch
  to (kept AEQ catch + return), where AEQ are adult equivalents.
  Function should be of the form `function(NO, HO, m) return(u)`. *Not
  used if `SOM@UPT_complex` is specified.*

- `u_terminal`:

  Numeric, matrix `[nsim, proyears]`, or function. If `type_T = "u"`,
  the harvest rate (ratio of kept catch to of the terminal marine
  fishery. Function should be of the form
  `function(NO, HO, m) return(u)`. *Not used if `SOM@UT_complex` is
  specified.*

- `K_PT`:

  Numeric or function. If `type_PT = "catch"`, the catch target of the
  immature component of the population in the pre-terminal fishery.
  Function should be of the form `function(NO, HO, m) return(K)`. *Not
  used if `SOM@UPT_complex` is specified.*

- `K_T`:

  Numeric or function. If `type_T = "catch"`, the catch target of the
  return in the terminal fishery. Function should be of the form
  `function(NO, HO, m) return(K)`. *Not used if `SOM@UT_complex` is
  specified.*

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
  proportion of fishing intensity experienced by each age class.

- `vulT`:

  Vector length `maxage` or matrix `[nsim, maxage]`. Vulnerability
  schedule (between 0-1) in the terminal fishery. Values indicate the
  proportion of fishing intensity experienced by each age class.

- `ForeErr_PT`:

  Numeric or matrix `[nsim, proyears]`. Multiplicative forecast error of
  the juvenile abundance in the preterminal fishery, i.e., observation
  error. Only used if `u_preterminal` or `K_PT` is a function where
  fishing intensity is determined by abundance. Default is 1 (no error).

- `ForeErr_T`:

  Numeric or matrix `[nsim, proyears]`. Multiplicative forecast error of
  the return size for the terminal fishery, i.e., observation error.
  Only used if `u_terminal` or `K_T` is a function where fishing
  intensity is determined by return size. Default is 1 (no error).

## Creating Object

Objects can be created by calls of the form `new("Harvest")`

## Examples

``` r
showClass("Harvest")
#> Class "Harvest" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                   
#> Name:                 Name             type_PT              type_T
#> Class:           character           character           character
#>                                                                   
#> Name:        u_preterminal          u_terminal                K_PT
#> Class: num.matrix.function num.matrix.function        num.function
#>                                                                   
#> Name:                  K_T              MSF_PT               MSF_T
#> Class:        num.function             logical             logical
#>                                                                   
#> Name:         release_mort               vulPT                vulT
#> Class:             numeric          num.matrix          num.matrix
#>                                               
#> Name:           ForeErr_PT           ForeErr_T
#> Class:          num.matrix          num.matrix
#> 
#> Extends: "Harvest.list"
```
