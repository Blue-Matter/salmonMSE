# All-H Analyzer

Wrapper function for an implementation of All-H Analyzer
([AHA](https://www.streamnet.org/home/data-maps/hatchery-reform/hsrg-tools/))
in R. Can be used to compare outputs between AHA and salmonMSE.

## Usage

``` r
AHA(SOM, ngen = 100, silent = FALSE)
```

## Arguments

- SOM:

  An object of class
  [SOM](https://docs.salmonmse.com/reference/SOM-class.md)

- ngen:

  Integer, the number of generations for which to run the simulation

- silent:

  Logical, indicates whether to silence messages to the R console

## Value

A named list containing vectors of state variables (by simulation,
population, and generation). See
[SMSE](https://docs.salmonmse.com/reference/SMSE-class.md) object
description.

## References

Hatchery Scientific Review Group. 2020. All-H Analyzer Tool Guide and
Documentation. May 2020.
