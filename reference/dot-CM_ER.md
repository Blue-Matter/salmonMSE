# Calculate exploitation rate in conditioning model

Internal functions for calculating aggregate exploitation rate
(`.CM_ER()`) or by age class (`.CM_ER()`).

## Usage

``` r
.CM_ER(
  report,
  type = c("PT", "T", "all"),
  r = 1,
  index_AEQ = NULL,
  brood = FALSE,
  simplify = TRUE
)

.CM_ERage(report, type = c("PT", "T"), brood = FALSE, simplify = TRUE)
```

## Arguments

- report:

  List, output of state variables from individual MCMC samples, obtained
  with
  [`get_report()`](https://docs.salmonmse.com/reference/CMfigures.md)

- type:

  Character, indicates type of variable to plot

- r:

  Integer, the release strategy for the figure (only if
  `annual = FALSE`)

- index_AEQ:

  Optional integer vector to identify years from which to borrow natural
  mortality and maturity to calculate adult equivalents for incomplete
  brood years. Only used if `at_age = FALSE`.

- brood:

  Logical, whether to show results by brood year or return year (FALSE)

- simplify:

  Logical, will return a matrix or array if TRUE, otherwise a list of
  output by MCMC simulation

## Value

An array if `simplify = TRUE`. Otherwise, a list.
