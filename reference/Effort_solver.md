# Solve for fishing effort

Internal solver used by
[`catch_func()`](https://docs.salmonmse.com/reference/catch_func.md) to
calculate the fishing effort needed to achieve target harvest rate or
catch rate, subject to partial retention due to mark-selective fishing.
Harvest rate is discounted by adult equivalents for preterminal
fisheries.

## Usage

``` r
Effort_solver(
  Eff,
  N,
  vul,
  ret,
  release_mort,
  type = c("u", "catch"),
  u = 0,
  K = 0,
  AEQ = array(1, dim(N)),
  p_mature = array(1, dim(N))
)
```

## Arguments

- Eff:

  Numeric, fishing effort

- N:

  Array `[ns, nage, n_r]`, total abundance (juvenile abundance for
  preterminal, return for terminal)

- vul:

  Array `[ns, nage, n_r]`, fishery vulnerability

- ret:

  Vector `[ns]`, retention rate

- release_mort:

  Vector `[ns]`, release mortality as a proportion, between 0-1. Only
  relevant if `ret < 1`.

- type:

  Character, either `"catch"`, or `"u"`, whether to solve for kept catch
  or harvest rate, respectively

- u:

  Numeric, harvest rate target

- K:

  Numeric, kept catch target

- AEQ:

  Array `[ns, nage, n_r]`, adult equivalents of catch

- p_mature:

  Array `[ns, nage, n_r]`, proportion mature by age class (used to
  calculate adult equivalent escapement)

## Value

Numeric. Returns the difference between the realized harvest rate (for
the given value of `Eff`) and the target
