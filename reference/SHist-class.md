# Class `"SHist"`

Stores the outputs from the historical reconstruction of salmon
operating models.

## Slots

- `Name`:

  Character. Identifying name

- `nyears`:

  Integer. The number of historical years

- `nsim`:

  Integer. The number of simulations

- `nstocks`:

  Integer. The number of stocks

- `Snames`:

  Character. Stock names

- `Egg_NOS`:

  Array `[nsim, nstocks, nyears]`. Spawning output, i.e., egg
  production, of natural origin spawners.

- `Egg_HOS`:

  Array `[nsim, nstocks, nyears]`. Spawning output of hatchery origin
  spawners.

- `Smolt`:

  Array `[nsim, nstocks, nyears]`. Natural smolt production (sum of
  offspring of natural and hatchery spawners).

- `Smolt_Rel`:

  Array `[nsim, nstocks, proyears]`. Smolts that are offspring of
  broodtake, i.e., hatchery releases.

- `Njuv_NOS`:

  Array `[nsim, nstocks, nage, nyears]`. Abundance of juvenile natural
  origin fish at the beginning of the year.

- `Njuv_HOS`:

  Array `[nsim, nstocks, nage, nyears]`. Abundance of juvenile hatchery
  origin fish at the beginning of the year.

- `Return_NOS`:

  Array `[nsim, nstocks, nage, nyears]`. Mature fish that will be
  natural origin spawners.

- `Return_HOS`:

  Array `[nsim, nstocks, nage, nyears]`. Mature fish that will be
  hatchery origin spawners.

- `Escapement_NOS`:

  Array `[nsim, nstocks, nage, nyears]`. The escapement of mature fish
  that will be natural origin spawners.

- `Escapement_HOS`:

  Array `[nsim, nstocks, nage, nyears]`. The escapement of mature fish
  that will be hatchery origin spawners.

- `NOS`:

  Array `[nsim, nstocks, proyears]`. Natural origin spawners.

- `HOS`:

  Array `[nsim, nstocks, proyears]`. Hatchery origin spawners.

- `HOS_effective`:

  Array `[nsim, nstocks, proyears]`. Hatchery origin spawners discounted
  by `gamma`.

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

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery harvest rate
  (from kept catch) of natural origin spawners.

- `UT_NOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery harvest rate of
  natural origin spawners.

- `UPT_HOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery harvest rate
  of hatchery origin spawners.

- `UT_HOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery harvest rate of
  hatchery origin spawners.

- `ExPT_NOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery exploitation
  rate (from kept catch and dead releases) of natural origin spawners.

- `ExT_NOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery exploitation rate
  of natural origin spawners.

- `ExPT_HOS`:

  Array `[nsim, nstocks, proyears]`. Pre-terminal fishery exploitation
  rate of hatchery origin spawners.

- `ExT_HOS`:

  Array `[nsim, nstocks, proyears]`. Terminal fishery exploitation rate
  of hatchery origin spawners.

- `Misc`:

  List. Miscellaneous output

## Examples

``` r
showClass("SHist")
#> Class "SHist" [package "salmonMSE"]
#> 
#> Slots:
#>                                                                   
#> Name:            Name         nyears           nsim        nstocks
#> Class:      character        numeric        numeric        numeric
#>                                                                   
#> Name:          Snames        Egg_NOS        Egg_HOS          Smolt
#> Class:      character          array          array          array
#>                                                                   
#> Name:       Smolt_Rel       Njuv_NOS       Njuv_HOS     Return_NOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:      Return_HOS Escapement_NOS Escapement_HOS            NOS
#> Class:          array          array          array          array
#>                                                                   
#> Name:             HOS  HOS_effective        KPT_NOS         KT_NOS
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
#> Name:        ExPT_HOS        ExT_HOS           Misc
#> Class:          array          array           list
```
