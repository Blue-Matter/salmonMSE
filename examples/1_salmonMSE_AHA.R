
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
  gamma = 0.8,
  pmax_NOB = 0.7,
  ptarget_NOB = 0.51,
  premove_HOS = 0.8,
  theta = c(100, 80, 70),
  rel_loss = c(0.5, 0.4, 0.1),
  fec_brood = 5040,
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
  prod_smolt_improve = 1         # Keep productivity (SR alpha) constant
)

Harvest <- new(
  "Harvest",
  u = sum(c(0.038, 0.025, 0, 0.14)),        # Specify fixed harvest rate of mature fish
  m = 0                                     # Mark rate of hatchery releases (future feature)
)

# Stitched salmon operating model
SOM <- new("SOM",
           proyears = 50,
           Bio, Habitat, Hatchery, Harvest)

# run salmonMSE
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

compare(SAHA$Fry_NOS, SMSE@Fry_NOS[, 1, ], "Fry_NOS", ylim = c(0, 150000))

#png("man/figures/example-NOS.png", height = 3, width = 7, res = 300, units = 'in')
compare(SAHA$NOS, SMSE@NOS[, 1, ], "NOS", ylim = c(0, 60), yline = 28.3)
#dev.off()
