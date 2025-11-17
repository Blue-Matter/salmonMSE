# Estimation function for conditioning model

Population dynamics model of an age structured salmon population. Used
with RTMB to estimate historical reconstruction from data.

## Usage

``` r
CM_int(p, d)
```

## Arguments

- p:

  List of parameter variables. See
  [`fit_CM()`](https://docs.salmonmse.com/reference/fit_CM.md).

- d:

  List of data variables. See
  [`fit_CM()`](https://docs.salmonmse.com/reference/fit_CM.md).

## Value

Numeric, objective function value (log-posterior)

## Author

Q. Huynh with Stan code provided by J. Korman and C. Walters
