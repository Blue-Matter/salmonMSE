# Introduction

## Background

[**salmonMSE**](https://salmonmse.com/) is a quantitative and stochastic
decision-support tool for Pacific salmon focusing on strategic
trade-offs among harvest, hatchery and habitat management levers.
salmonMSE can be used for risk-based analyses to evaluate the
performance of and prioritize management actions and identify trade-offs
towards achieving biological and harvest objectives

Initial development of salmonMSE is based on the [All-H
Analyzer](https://www.streamnet.org/wp-content/uploads/2021/01/3.-All-H-Analyzer-Guide-and-Documentation-051120.pdf)
(AHA) spreadsheet which models salmon dynamics of individual stocks by
life stages over successive generations to obtain the long-term
equilibrium properties of the state dynamics.

Currently, salmonMSE is intended to be used similarly, but expands upon
the population dynamics modeling in several ways. To accommodate more
complex life histories, salmonMSE is age-structured which is useful for
species that have multiple cohorts overlapping in the return.
Multi-stock models can be developed to evaluate outcomes at a stock
management unit. Finally, stochasticity can be incorporated in the model
to incorporate uncertainty in our understanding of the stock
productivity.

[openMSE](https://openmse.com/) is the population dynamics model
underlying salmonMSE. While openMSE is intended for the evaluation of
management strategies in marine fisheries, the functionality of the
software can be re-purposed to model salmon life stages. Most state
variables are therefore directly modeled within openMSE, while
additional state variables that are needed for salmon life history are
calculated in salmonMSE. For more information, see the article on the
[population dynamics
model](https://docs.salmonmse.com/articles/model.md).

Use of openMSE allows for rigorous review of the source code by a wider
user base. However, typical users of salmonMSE will not need to see the
internal conversion between salmonMSE and openMSE operating models.
Inputs and outputs of salmonMSE use salmon-specific terminology.

## Getting started

A salmon operating model contains the parameters for the population
dynamics and the management levers to be implemented. In salmonMSE, an
object of class `SOM` can be created from constituent objects of class
`Bio`, `Habitat`, `Hatchery`, and `Harvest` as follows:

``` r
library(salmonMSE)

Bio <- new("Bio", ...)
Habitat <- new("Habitat", ...)
Hatchery <- new("Hatchery", ...)
Harvest <- new("Harvest", ...)

SOM <- new("SOM", Bio, Hatchery, Habitat, Harvest)
```

- The `Bio` class specifies the natural production, for example,
  maturity, fecundity, stock-recruit relationship, and marine survival.
- The `Habitat` class (optional) specifies freshwater survival from egg,
  fry, and smolt life stages as a series of density-dependent functions,
  with options for time-varying survival, for example, as a function of
  environmental or habitat mitigation/restoration actions.
- The `Hatchery` class specifies the parameters surrounding hatchery
  production, such as the number of target releases, removal of hatchery
  spawners to maintain high proportions of natural spawners, and the
  population fitness parameters arising from interbreeding of hatchery
  and natural spawners.
- The `Harvest` class specifies the exploitation rate and harvest
  control rules for the fishery.

Additional slots in the `SOM` class control the projections, for
example, the number of years and simulation replicates.

Details on the slots of the various S4 classes can be obtained by typing
[`class?SOM`](https://docs.salmonmse.com/reference/SOM-class.md) in the
R console.

The simulation can then run with the
[`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md)
function:

``` r
MSE <- salmonMSE(SOM)
```

The output is a class `SMSE` object containing the state variables and
some performance metrics pertaining to hatchery dynamics (fitness, PNI,
etc.) as arrays typically indexed by simulation, stock, age, and year.
For example, `SMSE@NOS` reports the natural origin spawners.

For convenience and comparison purposes, salmonMSE distributes an
implementation of AHA in R as well:

``` r
SAHA <- AHA(SOM)
```

The resulting output is a named list following the format of the `SMSE`
object, but indexed by generation instead of year.

More details are also provided in the
[example](https://docs.salmonmse.com/articles/example.md) article.

## Future development

Development of salmonMSE is currently in progress. Stochasticity is
incorporated in the stock-recruit relationship, which predicts
density-dependent smolt production from eggs, and in the SAR parameter
(survival to adult return). Multi-stock models are not yet supported.
