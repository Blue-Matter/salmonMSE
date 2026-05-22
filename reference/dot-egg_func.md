# Egg production function

Simple wrapper function to calculate hatchery egg production

## Usage

``` r
.egg_func(
  ptake = 1,
  N,
  gamma = 1,
  fec,
  p_female,
  s_prespawn,
  val = 0,
  opt = TRUE
)
```

## Arguments

- ptake:

  Numeric, proportion of spawners that spawn

- N:

  Numeric, spawners

- gamma:

  Numeric, relative reproductive success of spawners

- p_female:

  Numeric, proportion female

- s_prespawn:

  Numeric, survival of spawners prior to egg production

- val:

  Numeric, target egg production. Used to optimize for `ptake` if
  `opt = TRUE`

- opt:

  Logical, whether the function is used to optimize for `ptake`

## Value

Numeric
