# Calculate abundance from density-dependent mortality

Calculates the abundance of survivors after applying either a
Beverton-Holt or Ricker stock-recruit relationship.

## Usage

``` r
calc_SRR(N1, N2 = N1, p, capacity, type = c("BH", "Ricker", "HS"))
```

## Arguments

- N1:

  Numeric, the initial abundance that scales the density-independent
  survival term

- N2:

  Numeric, the initial abundance that scales the density-dependent
  survival term

- p:

  Numeric, the productivity parameter that sets the maximum survival as
  the initial abundance approaches zero

- capacity:

  Numeric, the capacity parameter that set the maximum survivors

- type:

  Character, the functional form of the stock-recruit relationship

## Value

Numeric, the abundance of survivors

## Details

The Beverton-Holt stock recruit relationship is of the following form:
\$\$\textrm{Smolt} = \dfrac{\alpha N_1}{1 + \beta N_2}\$\$ where
\\\alpha = P\\, \\\beta = P/C\\.

The Ricker stock recruit relationship is of the following form:
\$\$\textrm{Smolt} = \alpha N_1 \exp(-\beta N_2)\$\$ where \\\alpha =
P\\, \\\beta = P/(Ce)\\, \\e\\ is Euler's number.

Productivity \\P\\ is in terms of abundance per unit of \\N_1\\ and
\\N_2\\.

The hockey stick is of the following form:

\$\$ \textrm{Smolt} = \begin{cases} p N_1 &, N_1 \le \frac{N_1}{N_2}
\times C\\ \frac{N_1}{N_2} \times C &, \textrm{otherwise} \end{cases}
\$\$

## See also

[`calc_SRRpars()`](https://docs.salmonmse.com/reference/calc_SRRpars.md)
