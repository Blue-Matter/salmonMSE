# Ricker reference points

Compute reference points (Umsy, Smsy, and Sgen) from Ricker
stock-recruit function based on Scheuerell (2016).

## Usage

``` r
calc_Smsy_Ricker(loga, b)

calc_Umsy_Ricker(loga)

calc_Sgen_Ricker(loga, b)
```

## Arguments

- loga:

  Numeric, alpha parameter (returns per spawner) in the Ricker function:
  \\R=S\exp(\log(a)-bS)\\ where `S` is the number of spawners and `R` is
  the return

- b:

  Numeric, beta parameter

## Value

All three functions return a numeric

## References

Scheuerell, M.D. 2016. An explicit solution for calculating optimum
spawning stock size from Ricker’s stock recruitment model. PeerJ
4:e1623. [doi:10.7717/peerj.1623](https://doi.org/10.7717/peerj.1623)

## See also

[`calc_ref()`](https://docs.salmonmse.com/reference/calc_ref.md)
