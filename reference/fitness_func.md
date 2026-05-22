# Calculate mean phenotype and fitness

Quantitative genetics model of mean phenotypic trait and fitness of next
brood year.

`calc_zbar()` is the internal function that calculates mean phenotype.

## Usage

``` r
fitness_func(
  Egg_NOS,
  Egg_HOS,
  Egg_NOB,
  Egg_HOB,
  zbar_brood,
  fitness_args = list()
)

calc_zbar(
  Egg_NOS,
  Egg_HOS,
  Egg_NOB,
  Egg_HOB,
  zbar_prev,
  fitness_variance,
  theta,
  phenotype_variance,
  heritability
)
```

## Arguments

- Egg_NOS:

  Numeric, egg production by natural-origin spawners

- Egg_HOS:

  Numeric, egg production by hatchery-origin brood

- Egg_NOB:

  Numeric, egg production by natural-origin brood

- fitness_args:

  List, containing `fitness_variance`, `theta`, `phenotype_variance`,
  `heritability`, `rel_loss`, and `fitness_floor`

- zbar_prev:

  Numeric, mean phenotype of parents

- fitness_variance:

  Numeric, variance (omega-squared) of the fitness function

- theta:

  Numeric length 2, optimum phenotype value for the natural and hatchery
  environments, respectively.

- phenotype_variance:

  Numeric, variance (sigma-squared) of the phenotypic trait (theta)

- heritability:

  Numeric, heritability (h-squared) of the phenotypic trait

## Value

`fitness_func` returns a list:

- `zbar` Numeric length 2, mean phenotype of next generation to natural
  `zbar[1]` and hatchery `zbar[2]` environments

- `fitness` Numeric length 2, fitness of next generation to natural and
  hatchery environments

- `fitness_loss` Matrix `[2, 3]` penalty in survival of next generation
  due to fitness effects

`calc_zbar` returns numeric length 2, mean phenotype of next generation
