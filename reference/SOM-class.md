# Class `"SOM"`

An object containing all the parameters for a salmon operating model
(SOM).

## Slots

- `Name`:

  Character. Identifying name

- `nsim`:

  Integer. Number of simulations

- `nyears`:

  Integer. The number of historical years

- `proyears`:

  Integer. The number of projected years

- `seed`:

  Integer. A random seed to ensure users can reproduce results exactly

- `Bio`:

  [Bio](https://docs.salmonmse.com/reference/Bio-class.md) object
  informing biological parameters, natural production, and habitat
  effects. Provide a list of Bio objects for multi-population models.

- `Habitat`:

  [Habitat](https://docs.salmonmse.com/reference/Habitat-class.md)
  object containing management levers for habitat mitigation. Provide a
  list of Habitat objects for multi-population models.

- `Hatchery`:

  [Hatchery](https://docs.salmonmse.com/reference/Hatchery-class.md)
  object containing management levers for hatchery production. Provide a
  list of Hatchery objects for multi-population models.

- `Harvest`:

  [Harvest](https://docs.salmonmse.com/reference/Harvest-class.md)
  object containing management levers for harvest. Provide a list of
  Harvest objects for multi-population models.

- `Historical`:

  [Historical](https://docs.salmonmse.com/reference/Historical-class.md)
  object to inform historical reconstruction and informing starting
  abundance for the projection. Provide a list of Historical objects for
  multi-population models.

- `stray`:

  Matrix `[np, np]` where `np = length(Bio)` and row `p` indicates the
  re-assignment of hatchery fish to each population when they mature (at
  the recruitment life stage). For example,
  `SOM@stray <- matrix(c(0.75, 0.25, 0.25, 0.75), 2, 2)` indicates that
  75 percent of mature fish return to their natal river and 25 percent
  stray in both populations. By default, an identity matrix is used (no
  straying).

## Objects from the Class

Objects can be created by calls of the form
`new("SOM", Bio, Habitat, Hatchery, Harvest, Historical)`.
