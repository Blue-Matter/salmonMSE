
# A simple example, same as example 1, with mark-selective fishing (see Harvest object)
library(salmonMSE)

class?SOM # Definition of inputs

SAR <- 0.01
nsim <- 3
Bio <- new(
  "Bio",
  maxage = 3,
  p_mature = c(0, 0, 1),
  SRrel = "BH",
  capacity = 17250,           # Beverton-Holt asymptote. Not unfished capacity!!
  kappa = 3,                  # Productivity in recruits per spawner
  Mjuv_NOS = c(0, -log(SAR), 0),
  fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
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
  Mjuv_HOS = Bio@Mjuv_NOS,
  gamma = 0.8,
  m = 0.5,                          # Mark rate of hatchery releases
  pmax_esc = 1,                     # Maximum proportion of escapement (after en route mortality) that could be used as broodtake
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
  MSF_T = TRUE,
  release_mort = c(0.1, 0.1),
  vulPT = c(0, 0, 0),
  vulT = c(1, 1, 1)
)

# 1000 natural and hatchery spawners each for the first generation
Historical <- new(
  "Historical",
  HistSpawner_NOS = 1000,
  HistSpawner_HOS = 1000
)

# Stitched salmon operating model
SOM <- new("SOM",
           Bio, Habitat, Hatchery, Harvest, Historical,
           nsim = nsim, nyears = 2, proyears = 50)

# run salmonMSE
MOM <- SOM2MOM(SOM)
H <- MSEtool::SimulateMOM(MOM, parallel = FALSE)
MMSE <- salmonMSE(SOM, convert = FALSE)

N <- apply(MMSE@N, c(1, 2, 3, 5), sum)
N[1, 1, , ]

MMSE@Misc$MICE$M_ageArray[1, 1, , 1, ]

SMSE <- salmonMSE(SOM)
class?SMSE # Definitions of arrays

# run AHA - list of vectors by generation
SAHA <- AHA(SOM, ngen = 20)

# Compare AHA (x) and salmonMSE (y) output
compare <- function(x, y, ylab = "y", ylim, yline) {
  y[y < 1e-8] <- NA
  if (missing(ylim)) ylim <- range(x, y, na.rm = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
  matplot(t(x), xlab = "Generation", ylab = paste("AHA:", ylab), ylim = ylim, typ = 'o', pch = 1, lwd = 1, col = 1)
  if (!missing(yline)) abline(h = yline, lty = 2)
  matplot(t(y), xlab = "Projection year", ylab = paste("salmonMSE:", ylab),
          ylim = ylim, typ = 'o', pch = 1, lwd = 1, col = 1)
  if (!missing(yline)) abline(h = yline, lty = 2)
}

compare(SAHA$Egg_NOS, SMSE@Egg_NOS[, 1, ], "Egg_NOS", ylim = c(0, 150000))
NOS_AHA <- SAHA$NOS[, 1, ]
NOS_SMSE <- apply(SMSE@NOS[, 1, , ], c(1, 3), sum)
compare(NOS_AHA, NOS_SMSE, "NOS", ylim = c(0, 60), yline = 28.9731)


