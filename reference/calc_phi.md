# Calculate equilibrium quantities with life history groups

Calculate eggs/smolt or spawners/smolt based on life history parameters
(survival, maturity, fecundity)

## Usage

``` r
calc_phi(
  Mjuv,
  p_mature,
  p_female,
  fec,
  s_enroute = 1,
  n_g = 1,
  p_LHG,
  output = c("egg", "spawner")
)
```

## Arguments

- Mjuv:

  Matrix `[maxage, n_g]`, but can be a vector if `n_g = 1`. Juvenile
  natural mortality

- p_mature:

  Matrix `[maxage, n_g]`, but can be a vector if `n_g = 1`. Maturity at
  age

- p_female:

  Numeric. Proportion female

- fec:

  Matrix `[maxage, n_g]`, but can be a vector if `n_g = 1`. Fecundity at
  age. Only used if `output = "egg"`

- s_enroute:

  Numeric, en-route survival of escapement to spawning grounds

- n_g:

  Integer. Number of life history groups

- p_LHG:

  Vector length `n_g` of proportion of life history groups per recruit.
  Default is `rep(1/n_g, n_g)`

- output:

  Character to indicate the output units, e.g., "egg" returns eggs per
  smolt, and "spawner" returns spawners per smolt

## Value

Numeric, units depend on `"output"` argument
