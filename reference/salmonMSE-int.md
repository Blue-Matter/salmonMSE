# Internal salmonMSE functions for converting operating model inputs and outputs

- `SOM2MOM()` converts a salmon operating model
  ([SOM](https://docs.salmonmse.com/reference/SOM-class.md)) to a
  multi-stock operating model
  ([MSEtool::MOM](https://msetool.openmse.com/reference/MOM-class.html))

- `make_Stock()` creates the
  [MSEtool::Stock](https://msetool.openmse.com/reference/Stock-class.html)
  object (openMSE) corresponding to salmon life stage

- `make_Fleet()` creates the
  [MSEtool::Fleet](https://msetool.openmse.com/reference/Fleet-class.html)
  object (openMSE) corresponding to the fishery that interacts with the
  various salmon life stages

- `multiHist2SHist()` converts the openMSE historical reconstruction
  into a salmon Hist object
  ([SHist](https://docs.salmonmse.com/reference/SHist-class.md))

- `MMSE2SMSE()` converts the openMSE projection output, along with
  additional state variables recorded in
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md),
  into a salmon MSE object
  ([SMSE](https://docs.salmonmse.com/reference/SMSE-class.md))

- `make_Harvest_MMP()` creates a multi-stock management procedure for
  the harvest component of the operating model by specifying
  exploitation rates through updating the formal arguments for
  [`Harvest_MMP()`](https://docs.salmonmse.com/reference/Harvest_MMP.md)

[`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md) is
the wrapper function that coordinates the simulation and the output.

## Usage

``` r
make_Harvest_MMP(SOM, check = TRUE)

MMSE2SMSE(MMSE, SOM, Harvest_MMP, N, stateN, Ford, H, stateH)

SOM2MOM(SOM, check = TRUE)

make_Stock(
  SOM,
  s = 1,
  g = 1,
  r = 1,
  NOS = TRUE,
  stage = c("immature", "return", "escapement")
)

make_Fleet(SOM, s, NOS = TRUE, stage = c("immature", "return", "escapement"))

multiHist2SHist(multiHist, SOM, check = TRUE)
```

## Arguments

- SOM:

  An object of class
  [SOM](https://docs.salmonmse.com/reference/SOM-class.md)

- check:

  Logical, whether to check the `SOM` object using
  [`check_SOM()`](https://docs.salmonmse.com/reference/check_SOM.md)

- MMSE:

  Object of class
  [MSEtool::MMSE](https://msetool.openmse.com/reference/MMSE-class.html)
  returned from MSEtool

- Harvest_MMP:

  Optional harvest function created by `make_Harvest_MMP()`

- N:

  Data frame of natural origin abundance at age saved in the
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md)
  environment during the simulation

- stateN:

  Data frame of natural origin state variables saved in the
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md)
  environment during the simulation

- Ford:

  Data frame of phenotypic trait values saved in the
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md)
  environment during the simulation

- H:

  Data frame of hatchery origin abundance at age saved in the
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md)
  environment during the simulation

- stateH:

  Data frame of hatchery origin state variables saved in the
  [salmonMSE_env](https://docs.salmonmse.com/reference/salmonMSE_env.md)
  environment during the simulation

- s:

  Integer, the population integer for which to create the Stock or Fleet
  object

- g:

  Integer, the life history group for which to create the Stock object.
  Not relevant if `NOS = FALSE`

- r:

  Integer, the hatchery release group for which to create the Stock
  object. Not relevant if `NOS = TRUE`

- NOS:

  Logical, whether the Stock or Fleet object corresponds to natural
  origin or hatchery origin fish

- stage:

  Character indicating the corresponding salmon life stage of the Stock
  or Fleet object

- multiHist:

  Class multiHist object returned from MSEtool

## Value

`make_Harvest_MMP`: Function of class "MMP" by updating the formal
arguments for
[`Harvest_MMP()`](https://docs.salmonmse.com/reference/Harvest_MMP.md)

`MMSE2SMSE`: [SMSE](https://docs.salmonmse.com/reference/SMSE-class.md)
object

`SOM2MOM`:
[MSEtool::MOM](https://msetool.openmse.com/reference/MOM-class.html)
object

`make_Stock`: List containing a
[MSEtool::Stock](https://msetool.openmse.com/reference/Stock-class.html)
object and accompanying custom parameters list

`make_Stock`: List containing a
[MSEtool::Fleet](https://msetool.openmse.com/reference/Fleet-class.html)
object and accompanying custom parameters list

`multiHist2SHist`:
[SHist](https://docs.salmonmse.com/reference/SHist-class.md) object
