# Class `"SOM"`

An object containing all the parameters for a salmon operating model
(SOM).

## Slots

- `Name`:

  Character. Identifying name

- `nsim`:

  Integer. Number of simulations

- `proyears`:

  Integer. The number of projected years

- `seed`:

  Integer. A random seed to ensure users can reproduce results exactly.
  Not currently used.

- `Bio`:

  [Bio](https://docs.salmonmse.com/reference/Bio-class.md) object
  informing biological parameters and natural production. Provide a list
  of Bio objects for multi-population models.

- `Habitat`:

  [Habitat](https://docs.salmonmse.com/reference/Habitat-class.md)
  object containing management levers for controlling survival in the
  freshwater environment. Provide a list of Habitat objects for
  multi-population models.

- `Hatchery`:

  [Hatchery](https://docs.salmonmse.com/reference/Hatchery-class.md)
  object containing management levers for hatchery production and
  in-river removals. Provide a list of Hatchery objects for
  multi-population models.

- `Harvest`:

  [Harvest](https://docs.salmonmse.com/reference/Harvest-class.md)
  object containing management levers for marine harvest. Provide a list
  of Harvest objects for multi-population models.

- `Historical`:

  [Historical](https://docs.salmonmse.com/reference/Historical-class.md)
  object to inform historical reconstruction and informing starting
  abundance for the projection. Provide a list of Historical objects for
  multi-population models.

- `stray`:

  For multi-population models: matrix `[np, np]` where
  `np = length(Bio)` and row `p` indicates the re-assignment of hatchery
  fish to each population when they mature (at the recruitment life
  stage). For example,
  `SOM@stray <- matrix(c(0.75, 0.25, 0.25, 0.75), 2, 2)` indicates that
  75 percent of mature fish return to their natal river and 25 percent
  stray in both populations. By default, an identity matrix is used (no
  straying).

- `UPT_complex`:

  *Optional* For multi-population models: aggregate preterminal harvest
  rate on stock complex. Can be numeric, matrix `[nsim, proyears]`, or
  function of the form `function(NO, HO, m) return(u)`. Leave empty to
  `numeric(0)` to operate on individual population basis.

- `UT_complex`:

  *Optional* For multi-population models: aggregate terminal harvest
  rate on stock complex. Can be numeric, matrix `[nsim, proyears]`, or
  function of the form `function(NO, HO, m) return(u)` Leave empty to
  `numeric(0)` to operate on individual population basis.

- `MSF_PT_complex`:

  *Optional* For multi-population models: whether the preterminal
  fishery on stock complex is mark-selective. Default is `FALSE`

- `MSF_T_complex`:

  *Optional* For multi-population models: whether the terminal fishery
  on stock complex is mark-selective. Default is `FALSE`

- `ForeErr_PT_complex`:

  *Optional* For multi-population models: numeric or matrix
  `[nsim, proyears]`. Multiplicative forecast error for preterminal
  fishery (observation error). Only used if `UPT_complex` is a function
  where fishing intensity is determined by abundance. Default is 1 (no
  error).

- `ForeErr_T_complex`:

  *Optional* For multi-population models: numeric or matrix
  `[nsim, proyears]`. Multiplicative forecast error for terminal fishery
  (observation error). Only used if `UT_complex` is a function where
  fishing intensity is determined by return size. Default is 1 (no
  error).

## Objects from the Class

Objects can be created by calls of the form
`new("SOM", Bio, Habitat, Hatchery, Harvest, Historical)`.
