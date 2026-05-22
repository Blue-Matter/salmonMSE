# Create list of parameters

Internal functions that convert the operating model inputs into a list
or data frame of parameters to pass on to other internal functions.

## Usage

``` r
define_hatchery_args(SOM)

define_habitat_args(SOM)

define_fitness_args(SOM)

define_SRRpars(SOM)
```

## Arguments

- SOM:

  [SOM](https://docs.salmonmse.com/reference/SOM-class.md) operating
  model object

## Value

`define_hatchery_args()` returns a length `nstocks` of hatchery
parameters, in-river removals, and behavior of hatchery-origin fish in
the freshwater environment (fecundity, reproductive success, etc.)

`define_habitat_args()` returns a list of
[Habitat](https://docs.salmonmse.com/reference/Habitat-class.md)
objects.

`define_fitness_args()` returns a length `nstocks` of fitness
parameters.

`define_SRRpars_args()` returns a length `nstocks`, each of which is a
data frame of stock-recruit parameters.
