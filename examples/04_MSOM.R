
library(salmonMSE)

class?SOM # Definition of inputs

ns <- 2
nsim <- 3
SAR <- 0.01

Bio <- lapply(1:ns, function(s) {
  new(
    "Bio",
    maxage = 3,
    p_mature = c(0, 0, 1),
    SRrel = "BH",
    capacity = ifelse(s == 1, 17250, 1000),     # Beverton-Holt asymptote. Not unfished capacity!!
    kappa = ifelse(s == 1, 3, 2),                     # Productivity in recruits per spawner
    Mjuv_NOS = c(0, -log(SAR), 0),
    fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
    p_female = 0.49,
    s_enroute = 1
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
    phatchery = ifelse(s == 1, 0.8, 0),
    premove_HOS = 0,
    theta = c(100, 80),
    rel_loss = c(0.5, 0.4, 0.1),
    fec_brood = c(0, 0, 5040),
    fitness_type = c("Ford", "none"),
    zbar_start = if (s == 1) c(93.1, 92) else c(100, 92),
    fitness_variance = 100,
    phenotype_variance = 10,
    heritability = 0.5,
    fitness_floor = 0.5
  )
})


Habitat <- lapply(1:ns, function(s) {
  new(
    "Habitat",
    use_habitat = FALSE
  )
})


Harvest <- lapply(1:ns, function(s) {
  new(
    "Harvest",
    u_preterminal = 0,             # No pre-terminal fishery
    u_terminal = 0.203,            # Specify fixed harvest rate of mature fish
    MSF_PT = FALSE,
    MSF_T = FALSE,
    release_mort = c(0, 0),
    vulPT = c(0, 0, 0),
    vulT = c(1, 1, 1)
  )
})

Historical <- lapply(1:ns, function(s) {
  new(
    "Historical",
    HistSpawner_NOS = 100,
    HistSpawner_HOS = ifelse(s == 1, 100, 0)
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

out <- salmonMSE(SOM)
saveRDS(out, file = "examples/SMSE/SMSE_MSOM.rds")

report(out, dir = "examples/reports", filename = "MSOM")

if (FALSE) {

  png("man/figures/multipop_spawners.png", height = 3.5, width = 4.5, units = "in", res = 400)
  par(mar = c(5, 4, 1, 1))
  plot_spawners(out, s = 2, prop = TRUE) # Spawner composition
  dev.off()

  # As of June 26, 2025
  max(out@p_wild[, 2, ], na.rm = TRUE) %>% round(3) # 0.444
  max(out@p_wild[, 1, ], na.rm = TRUE) %>% round(3) # 0.626
}
