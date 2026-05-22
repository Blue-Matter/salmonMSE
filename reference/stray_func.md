# Internal stray function

Calculates the number of fish that move from donor to recipient
populations from a movement matrix

## Usage

``` r
stray_func(N, stray_matrix, m)
```

## Arguments

- N:

  Array `[ns, nage, n_r]` of hatchery-origin fish

- stray_matrix:

  Matrix `[ns, ns]` that specifies proportion that move from population
  in the i-th row to the j-th column. The diagonal informs non-straying
  proportion

- m:

  Vector `[ns]` mark rater of hatchery-origin fish by donor population

## Value

Named list:

- `N_remain` Abundance of hatchery-origin fish that do not stray, same
  dimension as `N`

- `N_stray` Abundance of strays by recipient population, same dimension
  as `N`

- `m_stray` Mark rate of strays by recipient population, vector `[ns]`
