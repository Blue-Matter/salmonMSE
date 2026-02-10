# Class `"SMSE"`

Stores the outputs from the simulation of salmon operating models.

## Details

In generation \\t\\, proportionate natural influence (PNI) is defined
as:

\$\$\textrm{PNI}\_t = \dfrac{p^\textrm{NOB}\_t}{p^\textrm{NOB}\_t +
p^\textrm{HOSeff}\_t}\$\$

with \\p^\textrm{HOSeff} = \textrm{HOSeff}/(\textrm{NOS} +
\textrm{HOSeff})\\.

The proportion of wild salmon is defined as:

\$\$p^{\textrm{WILD}}\_t = q^\textrm{HOScen}\_t
\dfrac{(q^\textrm{HOScen}\_{t-1})^2} {(q^\textrm{HOScen}\_{t-1})^2 +
2\gamma \times p^\textrm{HOScen}\_{t-1} q^\textrm{HOScen}\_{t-1} +
\gamma^2 (p^\textrm{HOScen}\_{t-1})^2}\$\$

where \\q = 1-p\\ and \\p^\textrm{HOScen} = \textrm{HOS}/(\textrm{NOS} +
\textrm{HOS})\\.

## Slots

- `Name`:

  Character. Identifying name

- `proyears`:

  Integer. The number of projected years

- `nsim`:

  Integer. The number of simulations

- `nstocks`:

  Integer. The number of stocks

- `Snames`:

  Character. Stock names

- `Egg_NOS`:

  Array `[nsim, nstocks, proyears]`. Spawning output, i.e., egg
  production, of natural origin spawners.

- `Egg_HOS`:

  Array `[nsim, nstocks, proyears]`. Spawning output of hatchery origin
  spawners.

- `Fry_NOS`:

  Array `[nsim, nstocks, proyears]`. Fry that are offspring of natural
  origin spawners.

- `Fry_HOS`:

  Array `[nsim, nstocks, proyears]`. Fry that are offspring of hatchery
  origin spawners.

- `Smolt_NOS`:

  Array `[nsim, nstocks, proyears]`. Smolts that are offspring of
  natural origin spawners.

- `Smolt_HOS`:

  Array `[nsim, nstocks, proyears]`. Smolts that are offspring of
  hatchery origin spawners.

- `Smolt_Rel`:

  Array `[nsim, nstocks, proyears]`. Smolts that are offspring of
  broodtake, i.e., hatchery releases.

- `Njuv_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Abundance of juvenile natural
  origin fish at the beginning of the year.

- `Njuv_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Abundance of juvenile
  hatchery origin fish at the beginning of the year.

- `Return_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Mature fish that will be
  natural origin spawners.

- `Return_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Mature fish that will be
  hatchery origin spawners.

- `Escapement_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. The escapement of mature fish
  that will be natural origin spawners.

- `Escapement_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. The escapement of mature fish
  that will be hatchery origin spawners.

- `NOB`:

  Array `[nsim, nstocks, proyears]`. Natural origin broodtake.

- `HOB`:

  Array `[nsim, nstocks, proyears]`. Hatchery origin broodtake (local +
  strays).

- `HOB_stray`:

  Array `[nsim, nstocks, proyears]`. Hatchery origin broodtake (strays
  only).

- `HOB_import`:

  Array `[nsim, nstocks, proyears]`. Imported hatchery origin broodtake
  used for hatchery production.

- `NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Natural origin spawners.

- `HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Hatchery origin spawners
  (local + strays).

- `HOS_stray`:

  Array `[nsim, nstocks, nage, proyears]`. Hatchery origin spawners
  (strays only).

- `HOS_effective`:

  Array `[nsim, nstocks, nage, proyears]`. Hatchery origin spawners
  (local + strays) discounted by `gamma`.

- `KPT_NOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery kept catch of
  natural origin spawners.

- `KT_NOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery kept catch of
  natural origin spawners.

- `KPT_HOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery kept catch of
  hatchery origin spawners.

- `KT_HOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery kept catch of
  hatchery origin spawners.

- `DPT_NOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery released catch
  (live and dead) of natural origin spawners.

- `DT_NOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery released catch
  (live and dead) of natural origin spawners.

- `DPT_HOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery released catch
  (live and dead) of hatchery origin spawners.

- `DT_HOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery released catch
  (live and dead) hatchery origin spawners.

- `UPT_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery harvest
  rate (from kept catch) of natural origin spawners.

- `UT_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Terminal fishery harvest rate
  of natural origin spawners.

- `UPT_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery harvest
  rate of hatchery origin spawners.

- `UT_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Terminal fishery harvest rate
  of hatchery origin spawners.

- `ExPT_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery
  exploitation rate (from kept catch and dead releases) of natural
  origin spawners.

- `ExT_NOS`:

  Array `[nsim, nstocks, nage, proyears]`. Terminal fishery exploitation
  rate of natural origin spawners.

- `ExPT_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery
  exploitation rate of hatchery origin spawners.

- `ExT_HOS`:

  Array `[nsim, nstocks, nage, proyears]`. Terminal fishery exploitation
  rate of hatchery origin spawners.

- `fitness`:

  Array `[nsim, nstocks, 2, proyears]`. Fitness of the population in the
  natural (1) and hatchery (2) environments.

- `pNOB`:

  Array `[nsim, nstocks, proyears]`. Proportion of natural fish in the
  brood.

- `pHOS_census`:

  Array `[nsim, nstocks, proyears]`. Proportion of spawners of hatchery
  origin, weighted by age class fecundity.

- `pHOS_effective`:

  Array `[nsim, nstocks, proyears]`. Proportion of spawners of hatchery
  origin, discounted by `gamma`, weighted by age class fecundity.

- `PNI`:

  Array `[nsim, nstocks, proyears]`. Proportionate natural influence,
  index of gene flow from hatchery to the natural environment.

- `p_wild`:

  Array `[nsim, nstocks, proyears]`. Proportion of wild spawners,
  natural spawners whose parents were also produced in the natural
  environment assuming non-assortative mating, defined under Canada's
  Wild Salmon Policy.

- `Mjuv_loss`:

  Array `[nsim, nstocks, nage, proyears]`. Realized juvenile natural
  mortality, which may differ from inputs due to fitness loss.

- `Misc`:

  List. Miscellaneous output:

  - `Ref` for reference points

  - `SHist` for the
    [SHist](https://docs.salmonmse.com/reference/SHist-class.md) object

  - `SOM` for the
    [SOM](https://docs.salmonmse.com/reference/SOM-class.md) object.

  - `LHG` list `nstocks` long containing state variables by life history
    group

## Creating Object

Objects can be created by calls of the form `new("SMSE")`

## References

Withler et al. 2018. Genetically Based Targets for Enhanced
Contributions to Canadian Pacific Chinook Salmon Populations. DFO Can.
Sci. Advis. Sec. Res. Doc. 2018/019. xii + 88 p.

## Examples

``` r
showClass("SMSE")
#> Class "SMSE" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                   
#> Name:            Name       proyears           nsim        nstocks
#> Class:      character        numeric        numeric        numeric
#>                                                                   
#> Name:          Snames        Egg_NOS        Egg_HOS        Fry_NOS
#> Class:      character          array          array          array
#>                                                                   
#> Name:         Fry_HOS      Smolt_NOS      Smolt_HOS      Smolt_Rel
#> Class:          array          array          array          array
#>                                                                   
#> Name:        Njuv_NOS       Njuv_HOS     Return_NOS     Return_HOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:  Escapement_NOS Escapement_HOS            NOB            HOB
#> Class:          array          array          array          array
#>                                                                   
#> Name:       HOB_stray     HOB_import            NOS            HOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:       HOS_stray  HOS_effective        KPT_NOS         KT_NOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:         KPT_HOS         KT_HOS        DPT_NOS         DT_NOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:         DPT_HOS         DT_HOS        UPT_NOS         UT_NOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:         UPT_HOS         UT_HOS       ExPT_NOS        ExT_NOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:        ExPT_HOS        ExT_HOS        fitness           pNOB
#> Class:          array          array          array          array
#>                                                                   
#> Name:     pHOS_census pHOS_effective      Mjuv_loss            PNI
#> Class:          array          array          array          array
#>                                     
#> Name:          p_wild           Misc
#> Class:          array           list
```
