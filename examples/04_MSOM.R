
library(salmonMSE)

class?SOM # Definition of inputs

SAR <- 0.01
ns <- 2
nsim <- 3
Bio <- lapply(1:ns, function(s) {
  new(
    "Bio",
    maxage = 3,
    p_mature = c(0, 0, 1),
    SRrel = "BH",
    capacity_smolt = ifelse(s == 1, 17250, 1000),     # Beverton-Holt asymptote. Not unfished capacity!!
    kappa = ifelse(s == 1, 3, 2),                     # Productivity in recruits per spawner
    Mjuv_NOS = c(0, -log(SAR), 0),
    fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
    p_female = 0.49,
    s_enroute = 1
    #strays = 0
  )
})

Hatchery <- lapply(1:ns, function(s) {
  new(
    "Hatchery",
    n_yearling = ifelse(s == 1, 10000, 0),             # Management lever. No hatchery if both this line and next line are zero
    n_subyearling = 0,              # Management lever. No hatchery if both this line and previous line are zero
    s_prespawn = 1,                 # Survival prior to spawning
    s_egg_smolt = 0.92,             # Survival of eggs in hatchery
    s_egg_subyearling = 1,
    Mjuv_HOS = Bio[[s]]@Mjuv_NOS,
    gamma = 0.8,
    m = 1,                          # Mark rate of hatchery releases
    pmax_esc = 1,                   # Maximum proportion of escapement (after en route mortality) that could be used as broodtake
    pmax_NOB = 0.7,
    ptarget_NOB = 0.51,
    phatchery = ifelse(s == 1, 0.8, 0)
    premove_HOS = 0,
    theta = c(100, 80),
    rel_loss = c(0.5, 0.4, 0.1),
    fec_brood = c(0, 0, 5040),
    fitness_type = c("Ford", "none"),
    zbar_start = if (s == 1) c(93.1, 92) else c(100, 92),
    fitness_variance = 10,
    selection_strength = 3,
    heritability = 0.5,
    fitness_floor = 0.5
  )
})


Habitat <- lapply(1:ns, function(s) {
  new(
    "Habitat",
    capacity_smolt_improve = 1,    # Keep carrying capacity (SR alpha/beta) constant
    kappa_improve = 1              # Keep productivity (SR alpha) constant
  )
})


Harvest <- lapply(1:ns, function(s) {
  new(
    "Harvest",
    u_preterminal = 0,             # No pre-terminal fishery
    u_terminal = 0.203,            # Specify fixed harvest rate of mature fish
    MSF = FALSE,
    release_mort = c(0.1, 0.1),
    vulPT = c(0, 0, 0),
    vulT = c(1, 1, 1)
  )
})

# Return of 1000 natural and hatchery fish each for the first generation
nyears <- 2

Historical <- lapply(1:ns, function(s) {
  HistN <- array(0, c(nsim, Bio[[s]]@maxage, nyears+1, 2))
  if (s == 1) {
    HistN[, 1, 1, ] <- HistN[, 2, 2, ] <- 1000/SAR
  } else {
    HistN[, 1, 1, 1] <- HistN[, 2, 2, 1] <- 100/SAR
  }

  new(
    "Historical",
    HistN = HistN
  )
})

# Stitched salmon operating model
SOM <- new("SOM",
           Bio, Habitat, Hatchery, Harvest, Historical,
           nsim = nsim,
           nyears = 2,
           proyears = 50)

# Stray
SOM@stray <- matrix(c(0.75, 0.25, 0, 1), 2, 2, byrow = TRUE)

MOM <- SOM2MOM(SOM)

out <- salmonMSE(SOM, convert = FALSE)
saveRDS(out, file = 'examples/MSOM.rds')



