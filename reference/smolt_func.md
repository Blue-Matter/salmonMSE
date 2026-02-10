# Predict smolt production

Internal function that predicts the natural origin or hatchery origin
smolt production from the escapement, and saves state variables to
[salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md).

- `smolt_func()` is the population dynamics function

- `makeRel_smolt()` generates a list for openMSE to use in the
  simulations

## Usage

``` r
smolt_func(
  Nage_NOS,
  Nage_HOS,
  Nage_stray,
  x = -1,
  y,
  output = c("natural", "hatchery"),
  s_enroute,
  p_female,
  fec,
  SRRpars,
  hatchery_args,
  fitness_args,
  habitat_args,
  stray_args,
  s,
  g,
  prop_LHG,
  r
)

makeRel_smolt(
  p_smolt = 1,
  s = 1,
  p_natural,
  p_hatchery = NULL,
  p_stray = NULL,
  maxage,
  output = c("natural", "hatchery"),
  s_enroute,
  p_female,
  fec,
  SRRpars,
  hatchery_args,
  fitness_args,
  habitat_args,
  stray_args,
  g,
  prop_LHG,
  r
)
```

## Arguments

- Nage_NOS:

  Matrix `[n_age, n_g]` of natural escapement from openMSE. `n_g` is the
  number of life history groups (sub groups within cohorts) that
  contribute to spawning.

- Nage_HOS:

  Matrix `[n_age, n_r]` of hatchery escapement from openMSE. `n_r` is
  the number of release strategies (sub groups within cohorts) that
  contribute to spawning.

- Nage_stray:

  Matrix `[nage, n_stray]` of hatchery escapement from donor populations
  that will stray to the recipient population. `stray_args$stray_matrix`
  determines the proportions of the donor escapement that will move to
  the recipient population.

- x:

  Integer, simulation number from openMSE

- y:

  Integer, simulation year (including historical years)

- output:

  Character, whether to predict the natural origin or hatchery origin
  smolt production

- s_enroute:

  Numeric, en route survival of the escapement to the spawning grounds

- p_female:

  Vector `[maxage]`, proportion female for calculating the egg
  production.

- fec:

  Array `[nsim, maxage, nyears+proyears]`. The fecundity at age schedule
  of spawners

- SRRpars:

  Data frame containing stock recruit parameters for natural smolt
  production from egg production. Column names include: SRrel, kappa,
  capacity, Smax, phi

- hatchery_args:

  Named list containing various arguments controlling broodtake and
  hatchery production. See details below.

- fitness_args:

  Named list containing various arguments controlling population fitness
  from hatchery production. Names include: fitness_type, omega2, theta,
  fitness_variance, heritability, zbar_start, fitness_floor, rel_loss

- habitat_args:

  Named list with arguments sontrolling freshwater life stage survival.
  Names include: use_habitat (logical), Habitat
  ([Habitat](https://docs.salmonmse.com/reference/Habitat-class.md)
  object)

- stray_args:

  Named list with arguments controlling strays in and out of the local
  population Names include:

- s:

  Integer, salmonMSE population index. Used to report variables to
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md).

- g:

  Integer for the life history group of natural origin fish to pass the
  parameter back to openMSE (if `output = "natural"`)

- prop_LHG:

  Numeric vector, proportion of the egg production assign to life
  history groups (sums up to one) (only used if `output = "natural"`)

- r:

  Integer for the release strategy of hatchery origin fish to pass the
  parameter back to openMSE (if `output = "hatchery"`)

- p_smolt:

  Integer, the population index for the smolt production in the openMSE
  model, corresponding to `output`

- p_natural:

  Integer vector, the population index for the natural origin escapement
  in the openMSE model. Can be more than one if spawning is from
  multiple life history groups

- p_hatchery:

  Integer vector, the population index for the hatchery origin
  escapement in the openMSE model. Can be more than one if there are
  multiple release strategies. Set to `NULL` for no hatchery production

- p_stray:

  Integer vector, population index for the hatchery strays in
  multi-system models

## Value

- `smolt_func()` returns a numeric for the ratio of the realized smolt
  production vs. the hypothetical value if there were no hatchery, en
  route mortality, or habitat improvement

- `makeRel_smolt()` returns a list that is passed to openMSE as a
  inter-population relationship

## hatchery_args

Hatchery control parameters are included in a named list with the
following arguments:

- `egg_target` Numeric, target egg production for hatchery. Set to zero
  for no hatchery production.

- `ptarget_NOB` Numeric, the target proportion of the natural origin
  broodtake relative to the overall broodtake

- `pmax_NOB` Numeric, the maximum proportion of the natural origin
  escapement to be used as broodtake

- `fec_brood` Array `[nsim, maxage, nyears+proyears]`. the fecundity at
  age schedule of broodtake to calculate the total hatchery egg
  production

- `s_yearling` Numeric, the survival of eggs to the smolt life stage
  (for yearling release)

- `s_subyearling` Numeric. the survival of eggs to subyearling life
  stage (for subyearling release)

- `p_yearling` Numeric, the proportion of annual releases at the
  yearling life stage (vs. subyearling)

- `phatchery` Numeric, the proportion of the hatchery origin escapement
  that return to the hatchery, for example, by removal from spawning
  grounds or swim-in facilities. These fish are available for broodtake.

- `premove_HOS` Numeric or function, the proportion of the hatchery
  origin fish to be removed from the spawning grounds (in order to
  ensure a high proportion of NOS or for an in-river fishery). These
  fish are not available for broodtake. For example, a value less than
  one can represent imperfect implementation of weir removal.

- `premove_NOS` Numeric or function, the proportion of the natural
  origin fish to be removed from the spawning grounds, for example, an
  in-river fishery.

- `s_prespawn` Numeric, the survival of broodtake prior to egg
  production. `1 - s_prespawn` is the proportion of fish not used for
  hatchery purposes, e.g., mortality or other resesarch purposes. Used
  to back-calculate the broodtake.

- `gamma` Numeric, the relative reproductive success of hatchery origin
  spawners (relative to natural origin spawners)

- `m` Numeric, mark rate for selective broodtake
