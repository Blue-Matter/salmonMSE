

#### Install packages ----

#remotes::install_github("Blue-Matter/MSEtool")
#remotes::install_github("Blue-Matter/salmonMSE")

#### Simple example, used to compare salmonMSE and AHA outputs ----
library(salmonMSE)

class?SOM # Definition of inputs

SAR <- 0.01
nsim <- 3
Bio <- new(
  "Bio",
  maxage = 3,
  p_mature = c(0, 0, 1),
  SRrel = "BH",
  capacity = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
  kappa = 3,                  # Productivity in recruits per spawner
  Mjuv_NOS = c(0, -log(SAR), 0),
  fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
  p_female = 0.49,
  s_enroute = 1
  #strays = 0
)

Hatchery <- new(
  "Hatchery",
  n_yearling = 10000,             # Management lever. No hatchery if both this line and next line are zero
  n_subyearling = 0,              # Management lever. No hatchery if both this line and previous line are zero
  s_prespawn = 1,                 # Survival prior to spawning
  s_egg_smolt = 0.92,             # Survival of eggs in hatchery
  s_egg_subyearling = 1,
  Mjuv_HOS = Bio@Mjuv_NOS,
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

SMSE <- salmonMSE(SOM)
saveRDS(SMSE, "examples/SMSE/SMSE_simple.rds")
class?SMSE

# run AHA - list of vectors by generation
SAHA <- AHA(SOM, ngen = 20)

# Compare AHA (x) and salmonMSE (y) output
if (FALSE) {

  compare <- function(x, y, ylab = "y", ylim, yline) {
    y[y < 1e-8] <- NA
    if (missing(ylim)) ylim <- range(x, y, na.rm = TRUE)
    par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
    matplot(t(x), xlab = "Generation", ylab = paste("AHA:", ylab), ylim = ylim, type = 'o', pch = 1, lwd = 1, col = 1,
            panel.first = graphics::grid())
    if (!missing(yline)) abline(h = yline, lty = 2)
    matplot(t(y), xlab = "Projection year", ylab = paste("salmonMSE:", ylab),
            ylim = ylim, typ = 'o', pch = 1, lwd = 1, col = 1,
            panel.first = graphics::grid())
    if (!missing(yline)) abline(h = yline, lty = 2)
  }

  compare(SAHA$Egg_NOS[, 1, ], SMSE@Egg_NOS[, 1, ], "Egg_NOS", ylim = c(0, 200000))

  png("man/figures/example-NOS.png", height = 3, width = 7, res = 300, units = 'in')
  NOS_AHA <- SAHA$NOS[, 1, ]
  NOS_SMSE <- apply(SMSE@NOS[, 1, , ], c(1, 3), sum)
  compare(NOS_AHA, NOS_SMSE, "NOS", ylim = c(0, 100), yline = SAHA$NOS[1, 1, 20])
  dev.off()

}

#### MSF example ----
Harvest_MSF <- new(
  "Harvest",
  u_preterminal = 0,             # No pre-terminal fishery
  u_terminal = 0.203,            # Specify fixed harvest rate of mature fish
  MSF_PT = FALSE,
  MSF_T = TRUE,
  release_mort = c(0.1, 0.1),
  vulPT = c(0, 0, 0),
  vulT = c(1, 1, 1)
)

# Stitched salmon operating model
SOM_MSF <- new("SOM",
               Bio, Habitat, Hatchery, Harvest_MSF, Historical,
               nsim = nsim, nyears = 2, proyears = 50)

SMSE_MSF <- salmonMSE(SOM_MSF)
saveRDS(SMSE_MSF, "examples/SMSE/SMSE_MSF.rds")



#### Stochastic example ----
SAR <- 0.01
nsim_stochastic <- 100

set.seed(100)
kappa_mean <- 3
kappa_sd <- 0.3
kappa <- rlnorm(nsim_stochastic, log(kappa_mean) - 0.5 * kappa_sd^2, kappa_sd)
Bio_stochastic <- new(
  "Bio",
  maxage = 3,
  p_mature = c(0, 0, 1),
  SRrel = "BH",
  capacity = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
  kappa = kappa,                  # Productivity in recruits per spawner
  Mjuv_NOS = c(0, -log(SAR), 0),
  fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
  p_female = 0.49,
  s_enroute = 1
  #strays = 0
)

# Run simulation
SOM_stochastic <- new("SOM",
                      Bio_stochastic, Habitat, Hatchery, Harvest, Historical,
                      nsim = nsim_stochastic, nyears = 2, proyears = 50)

tictoc::tic()
SMSE_stochastic <- salmonMSE(SOM_stochastic)
tictoc::toc()

# Save simulation object
saveRDS(SMSE_stochastic, "examples/SMSE/SMSE_stochastic.rds")

# Make markdown report
report(SMSE_stochastic, dir = "examples/reports", filename = "SMSE_stochastic")

if (FALSE) {

  # Plot productivity
  png("man/figures/example-kappa.png", height = 4, width = 6, res = 300, units = 'in')
  par(mfrow = c(1, 1), mar = c(5, 4, 1, 1))
  hist(kappa, main = NULL)
  dev.off()

  # Plot PNI time series
  png("man/figures/example-PNI-ts.png", height = 4, width = 6, res = 300, units = 'in')
  par(mar = c(5, 4, 1, 1))
  plot_statevar_ts(SMSE_stochastic, "PNI", quant = TRUE)
  dev.off()

  # Plot PNI histogram
  png("man/figures/example-PNI-hist.png", height = 4, width = 6, res = 300, units = 'in')
  par(mar = c(5, 4, 1, 1))
  plot_statevar_hist(SMSE_stochastic, "PNI", y = 49, breaks = 10, xlim = c(0.5, 0.9))
  dev.off()

  # Probability PNI > 0.8
  PNI_LT <- SMSE_stochastic@PNI[, 1, 49]
  mean(PNI_LT >= 0.8)
  quantile(PNI_LT, c(0.025, 0.5, 0.975))

  # PNI vs kappa
  png("man/figures/example-PNI-kappa.png", height = 4, width = 6, res = 300, units = 'in')
  par(mar = c(5, 4, 1, 1))
  plot(kappa, PNI_LT, xlim = c(1, 7), ylim = c(0.5, 1), panel.first = grid(), pch = 16)
  dev.off()

}
