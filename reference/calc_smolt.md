# Smolt production

Calculate smolt production from base stock-recruit parameters and
fitness loss

## Usage

``` r
calc_smolt(
  N1,
  N2 = N1,
  kappa,
  capacity,
  Smax,
  phi = 1,
  fitness_loss = 1,
  SRrel = c("BH", "Ricker"),
  per_recruit = FALSE
)
```

## Arguments

- N1:

  Egg production for the density-independent component of the
  stock-recruit relationship. Can be the number of spawners if `phi = 1`
  and `Smax` is in units of spawners.

- N2:

  Egg production for the density-dependent component of the
  stock-recruit relationship (only used if `per_recruit = FALSE`)

- kappa:

  Base productivity parameter

- capacity:

  Base capacity parameter if `SRrel = "BH"`

- Smax:

  Base Smax parameter if `SRrel = "Ricker"`

- phi:

  Unfished egg per smolt (`1/phi` is the replacement line)

- fitness_loss:

  Survival term to reduce smolt production due to fitness, between 0-1

- SRrel:

  Character for the stock-recruit function

- per_recruit:

  Logical, whether N1 is a per recruit quantity (TRUE) or in absolute
  numbers (FALSE)

## Value

Numeric
