# Run salmonMSE

`salmonMSE()` runs a salmon management strategy evaluation from an
operating model object
([SOM](https://docs.salmonmse.com/reference/SOM-class.md)), by checking
the operating model object with
[`check_SOM()`](https://docs.salmonmse.com/reference/check_SOM.md),
running the projection in `ProjectSOM()` (parallel if called upon, then
stitches together the output in a single object), and calculates
reference points with
[`calc_ref()`](https://docs.salmonmse.com/reference/calc_ref.md).

`ProjectSOM()` is the internal projection function.

## Usage

``` r
salmonMSE(SOM, ncores = 1, silent = FALSE)

ProjectSOM(SOM, sims, check = FALSE)
```

## Arguments

- SOM:

  An object of class
  [SOM](https://docs.salmonmse.com/reference/SOM-class.md)

- ncores:

  Integer, maximum number of processors to run projection with parallel
  processing

- silent:

  Logical, whether to report progress in console

- sims:

  Optional integer vector to run projection for a subset of simulations.
  Intended for parallel processing.

- check:

  Logical, whether to check the structure of the input object with
  [`check_SOM()`](https://docs.salmonmse.com/reference/check_SOM.md)

## Value

[SMSE](https://docs.salmonmse.com/reference/SMSE-class.md) object

## Examples

``` r
if (FALSE) { # \dontrun{
SMSE <- salmonMSE(simple_SOM)
} # }
```
