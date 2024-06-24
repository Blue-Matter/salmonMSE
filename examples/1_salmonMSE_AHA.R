
# Compare salmonMSE and AHA outputs

# Install packages
remotes::install_github("Blue-Matter/MSEtool")
remotes::install_github("Blue-Matter/salmonMSE")


library(salmonMSE)

class?SOM # Definition of inputs

Bio <- new(
  "Bio",
  nsim = 3,
  maxage = 3,
  p_mature = c(0, 0, 1),
  SRrel = "BH",
  capacity_smolt = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
  kappa = 3,                  # Productivity in recruits per spawner
  Mocean_NOS = c(0, -log(0.01), 0),
  fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
  p_female = 0.49
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
  Mocean_HOS = Bio@Mocean_NOS,
  gamma = 0.8,
  pmax_NOB = 0.7,
  ptarget_NOB = 0.51,
  phatchery = 0.8,
  premove_HOS = 0,
  theta = c(100, 80),
  rel_loss = c(0.5, 0.4, 0.1),
  fec_brood = c(0, 0, 5040),
  fitness_type = "Ford",
  pbar_start = c(93.1, 92),
  fitness_variance = 10,
  selection_strength = 3,
  heritability = 0.5,
  fitness_floor = 0.5
)

Habitat <- new(
  "Habitat",
  capacity_smolt_improve = 1,    # Keep carrying capacity (SR alpha/beta) constant
  kappa_improve = 1              # Keep productivity (SR alpha) constant
)

Harvest <- new(
  "Harvest",
  u_preterminal = 0,             # No pre-terminal fishery
  u_terminal = 0.203,            # Specify fixed harvest rate of mature fish
  m = 0,                         # Mark rate of hatchery releases
  release_mort = c(0.1, 0.1),
  vulPT = c(0, 0, 0),
  vulT = c(1, 1, 1)
)

HistN <- array(0, c(3, 3, 2, 2))
HistN[, 1, 1, 1] <- 1000

Historical <- new(
  "Historical",
  HistN = HistN
)


# Stitched salmon operating model
SOM <- new("SOM",
           nyears = 2,
           proyears = 50,
           Bio, Habitat, Hatchery, Harvest, Historical = Historical)

# run salmonMSE
SMSE <- salmonMSE(SOM, convert = FALSE)
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

compare(SAHA$Fry_NOS, SMSE@Fry_NOS[, 1, ], "Fry_NOS", ylim = c(0, 150000))

#png("man/figures/example-NOS.png", height = 3, width = 7, res = 300, units = 'in')
compare(SAHA$NOS, SMSE@NOS[, 1, ], "NOS", ylim = c(0, 60), yline = 28.3)
#dev.off()
