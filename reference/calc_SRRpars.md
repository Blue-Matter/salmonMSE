# Convert density-dependent survival parameters

Converts from capacity/productivity parameters to alpha/beta
stock-recruit parameters where productivity is in terms of smolts per
spawner and alpha is terms of smolts per egg.

## Usage

``` r
calc_SRRpars(p, capacity, f = 1, p_female = 1, type = c("BH", "Ricker", "HS"))
```

## Arguments

- p:

  Numeric, the productivity parameter that sets the maximum survival as
  the initial abundance approaches zero

- capacity:

  Numeric, the capacity parameter that set the maximum survivors

- f:

  Fecundity, the spawning output per mature female

- p_female:

  The proportion of females per spawner

- type:

  Character, the functional form of the stock-recruit relationship

## Value

Numeric vector length 2 for alpha and beta value, respectively

## Details

\$\$\alpha = \dfrac{P}{f \times p\_{female}}\$\$

For the Beverton-Holt stock recruit relationship: \$\$\beta =
\dfrac{\alpha}{C}\$\$

For the Ricker stock recruit relationship: \$\$\beta =
\dfrac{\alpha}{Ce}\$\$, \\e\\ is Euler's number.

## See also

[`calc_SRR()`](https://docs.salmonmse.com/reference/calc_SRR.md)
