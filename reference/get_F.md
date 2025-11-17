# Calculate F from harvest rate

Solves for apical instantaneous fishing mortality rate (F), proportional
to fishing effort, from harvest rate (total retained catch over total
abundance). The apical F can be greater than the realized F, if
retention \< 1.

## Usage

``` r
get_F(
  u = 0,
  K = 0,
  type = c("u", "catch"),
  M,
  N = 1,
  vul = 1,
  ret = 1,
  release_mort = 0,
  Fmax = 20
)
```

## Arguments

- u:

  Harvest rate, between 0-1

- K:

  Catch, between 0-Inf

- type:

  Character, either `"catch"`, or `"u"`, whether to solve for catch or
  harvest rate, respectively

- M:

  Instantaneous natural mortality rate

- N:

  Abundance

- vul:

  Vulnerability

- ret:

  Retention rate

- release_mort:

  Release mortality as a proportion, between 0-1. Only relevant if
  `ret < 1`.

- Fmax:

  Maximum allowable value of F

## Value

Numeric for the apical F
