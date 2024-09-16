
# Compare salmonMSE and AHA outputs

# Install packages
remotes::install_github("Blue-Matter/MSEtool", ref = "salmonMSE")
remotes::install_github("Blue-Matter/salmonMSE")


library(salmonMSE)

class?SOM # Definition of inputs

Bio <- new(
  "Bio",
  nsim = 3,
  maxage = 3,
  p_mature = c(0, 0, 1),
  capacity_smolt = 17250,   # Beverton-Holt asymptote. Not unfished capacity!!
  prod_smolt = 300,         # Productivity adult/SAR. At unfished equilibrium the realized smolt productivity = 1/SAR
  SAR = 0.01,               # Future feature to allow for time-varying (PDO forcing)
  fec = 5040,               # Spawning fecundity of NOS and HOS
  p_female = 0.49,
  s_enroute = 1
  #strays = 0
)

Hatchery <- new(
  "Hatchery",
  n_yearling = 0,                   # Management lever. No hatchery if both this line and next line are zero
  n_subyearling = 2e6,              # Management lever. No hatchery if both this line and previous line are zero
  #n_subyearling = 0,
  s_prespawn = 1,                   # Survival prior to spawning
  s_egg_smolt = 1e-6,               # Survival of eggs in hatchery
  s_egg_subyearling = 0.92,
  gamma = 0.8,
  pmax_NOB = 0.7,
  ptarget_NOB = 0.51,
  phatchery = 0.8,
  premove_HOS = 0,
  theta = c(100, 80),
  rel_loss = c(0.5, 0.4, 0.1),
  fec_brood = 5040,
  fitness_type = "Ford",
  zbar_start = c(93.1, 92),
  fitness_variance = 10,
  selection_strength = 3,
  heritability = 0.5,
  fitness_floor = 0.5
)

Habitat <- new(
  "Habitat",
  capacity_smolt_improve = 1,    # Keep carrying capacity (SR alpha/beta) constant
  prod_smolt_improve = 1         # Keep productivity (SR alpha) constant
)

Harvest <- new(
  "Harvest",
  u = sum(c(0.038, 0.025, 0, 0.14)),        # Specify fixed harvest rate of mature fish
  m = 0                                     # Mark rate of hatchery releases (future feature)
)

# Stitched salmon operating model
SOM <- new("SOM", Bio, Habitat, Hatchery, Harvest)

# Base OM. Small natural population. No hatchery.
SOM1 <- SOM
SAHA1 <- AHA(SOM1)

# Alternate OM with larger capacity. Test hatchery dynamics.
SOM2 <- SOM
SOM2@capacity_smolt <- 17250 * 1e3



