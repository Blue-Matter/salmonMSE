# Run salmonMSE

`salmonMSE()` runs a salmon management strategy evaluation through the
following steps:

- Converts a salmon operating model
  ([SOM](https://docs.salmonmse.com/reference/SOM-class.md)) to a
  multi-stock operating model
  ([MSEtool::MOM](https://msetool.openmse.com/reference/MOM-class.html))
  via
  [`SOM2MOM()`](https://docs.salmonmse.com/reference/salmonMSE-int.md)

- Creates a harvest management procedure specifying the harvest control
  rule

- Generates the historical reconstruction of the state variables

- Runs projection (if `Hist = FALSE`)

- Converts the openMSE output, along with additional state variables
  recorded in
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md),
  into a salmon MSE object (SMSE) via
  [`MMSE2SMSE()`](https://docs.salmonmse.com/reference/salmonMSE-int.md)

## Usage

``` r
salmonMSE(SOM, Hist = FALSE, silent = FALSE, trace = FALSE, convert = TRUE)
```

## Arguments

- SOM:

  An object of class
  [SOM](https://docs.salmonmse.com/reference/SOM-class.md)

- Hist:

  Logical, whether to stop the function stop after historical
  simulations?

- silent:

  Logical, whether to report progress in console

- trace:

  Logical, whether to report additional messages from openMSE

- convert:

  Logical, whether to convert the output into a salmon MSE (SHist or
  SMSE, depending on `Hist`) object

## Value

If `Hist = TRUE`: if `convert = TRUE`, a
[SHist](https://docs.salmonmse.com/reference/SHist-class.md) object or
if `convert = FALSE`, a multiHist object (list).

If `Hist = FALSE`: if `convert = TRUE`, a
[SMSE](https://docs.salmonmse.com/reference/SMSE-class.md) object or if
`convert = FALSE`, a
[MSEtool::MMSE](https://msetool.openmse.com/reference/MMSE-class.html)
object.
