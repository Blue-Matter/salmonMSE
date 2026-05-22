
library(salmonMSE)

class?SOM # Definition of inputs

SAR <- c(0.01, 0.02)
ns <- 2

maxage <- 3
nsim <- 3
proyears <- 50

Bio_list <- lapply(1:ns, function(s) {
  Mjuv_NOS <- array(0, c(nsim, maxage-1, proyears, ifelse(s == 1, 2, 1)))
  if (s == 1) {
    Mjuv_NOS[, 2, , 1] <- -log(0.01)
    Mjuv_NOS[, 2, , 2] <- -log(0.02)
  } else {
    Mjuv_NOS[, 2, , ] <- -log(0.01)
  }

  new(
    "Bio",
    maxage = maxage,
    n_g = ifelse(s == 1, 2, 1),
    p_LHG = if (s == 1) c(0.9, 0.1) else 1,
    p_mature = c(0, 0, 1),
    SRrel = "BH",
    capacity = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
    kappa = 3,                  # Productivity in recruits per spawner
    Mjuv_NOS = Mjuv_NOS,
    fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
    p_female = 0.49,
    s_enroute = 1
    #strays = 0
  )
})

Hatchery <- new(
  "Hatchery",
  n_yearling = 0,             # Management lever. No hatchery if both this line and next line are zero
  n_subyearling = 0,              # Management lever. No hatchery if both this line and previous line are zero
  s_prespawn = 1,                 # Survival prior to spawning
  s_egg_smolt = 0.92,             # Survival of eggs in hatchery
  s_egg_subyearling = 1,
  Mjuv_HOS = Bio_list[[1]]@Mjuv_NOS[, , , 1],
  gamma = 0.8,
  m = 1,                          # Mark rate of hatchery releases
  pmax_esc = 1,                   # Maximum proportion of escapement (after en route mortality) that could be used as broodtake
  pmax_NOB = 0.7,
  ptarget_NOB = 0.51,
  phatchery = 0.8,
  premove_HOS = 0,
  theta = c(100, 80),
  rel_loss = c(0.5, 0.4, 0.1),
  fec_brood = c(0, 0, 5040),
  fitness_type = c("Ford", "none"),
  zbar_start = c(93.1, 92),
  fitness_variance = 100,
  phenotype_variance = 10,
  heritability = 0.5,
  fitness_floor = 0.5
)


Habitat <- new(
  "Habitat",
  use_habitat = FALSE
)


Harvest <- new(
  "Harvest",
  u_preterminal = 0,             # No pre-terminal fishery
  u_terminal = 0.203,            # Specify fixed harvest rate of mature fish
  MSF_PT = FALSE,
  MSF_T = FALSE,
  release_mort = c(0.1, 0.1),
  vulPT = c(0, 0, 0),
  vulT = c(1, 1, 1)
)

# Return of 1000 natural and hatchery fish each for the first generation
Historical <- new(
  "Historical",
  InitNjuv_NOS = 1000,
  InitNjuv_HOS = 1000
)

# Stitched salmon operating model
SOM_2LHG <- new("SOM",
                Bio_list[[1]], Habitat, Hatchery, Harvest, Historical,
                nsim = nsim,
                proyears = proyears)

SOM_1LHG <- new("SOM",
                Bio_list[[2]], Habitat, Hatchery, Harvest, Historical,
                nsim = nsim,
                proyears = proyears)


SMSE_2LHG <- salmonMSE(SOM_2LHG)
SMSE_1LHG <- salmonMSE(SOM_1LHG)

compare_statevar_ts(list(SMSE_2LHG, SMSE_1LHG), var = "NOS", names = c("2 LHG", "1 LHG"))

plot_LHG(SMSE_2LHG)
plot_LHG(SMSE_2LHG, var = "Smolt")
