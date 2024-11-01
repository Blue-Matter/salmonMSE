
# Simple example, used to compare salmonMSE and AHA outputs

# Install packages
remotes::install_github("Blue-Matter/MSEtool")
remotes::install_github("Blue-Matter/salmonMSE")


library(salmonMSE)

class?SOM # Definition of inputs

SAR <- 0.01
nsim <- 3
Bio <- new(
  "Bio",
  maxage = 3,
  p_mature = c(0, 0, 1),
  SRrel = "BH",
  capacity_smolt = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
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
  MSF = FALSE,
  release_mort = c(0.1, 0.1),
  vulPT = c(0, 0, 0),
  vulT = c(1, 1, 1)
)

# Return of 1000 natural and hatchery fish each for the first generation
nyears <- 2
HistN <- array(0, c(nsim, Bio@maxage, nyears+1, 2))
HistN[, 1, 1, ] <- HistN[, 2, 2, ] <- 1000/SAR

Historical <- new(
  "Historical",
  HistN = HistN
)

# Stitched salmon operating model
SOM <- new("SOM",
           Bio, Habitat, Hatchery, Harvest, Historical,
           nsim = nsim, nyears = 2, proyears = 50)

# run salmonMSE
#MOM <- SOM2MOM(SOM)
#H <- MSEtool::SimulateMOM(MOM, parallel = FALSE)
#MMSE <- salmonMSE(SOM, convert = FALSE)

#N <- apply(MMSE@N, c(1, 2, 3, 5), sum)
#N[1, 1, , ]

#MMSE@Misc$MICE$M_ageArray[1, 1, , 1, ]

SMSE <- salmonMSE(SOM)
class?SMSE # Definitions of arrays

# run AHA - list of vectors by generation
SAHA <- AHA(SOM, ngen = 20)

# Compare AHA (x) and salmonMSE (y) output
compare <- function(x, y, ylab = "y", ylim, yline) {
  y[y < 1e-8] <- NA
  if (missing(ylim)) ylim <- range(x, y, na.rm = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))
  matplot(t(x), xlab = "Generation", ylab = paste("AHA:", ylab), ylim = ylim, typ = 'o', pch = 1, lwd = 1, col = 1,
          panel.first = graphics::grid())
  if (!missing(yline)) abline(h = yline, lty = 2)
  matplot(t(y), xlab = "Projection year", ylab = paste("salmonMSE:", ylab),
          ylim = ylim, typ = 'o', pch = 1, lwd = 1, col = 1,
          panel.first = graphics::grid())
  if (!missing(yline)) abline(h = yline, lty = 2)
}

compare(SAHA$Egg_NOS, SMSE@Egg_NOS[, 1, ], "Egg_NOS", ylim = c(0, 200000))

png("man/figures/example-NOS.png", height = 3, width = 7, res = 300, units = 'in')
compare(SAHA$NOS, SMSE@NOS[, 1, ], "NOS", ylim = c(0, 100), yline = SAHA$NOS[1, 20])
dev.off()



# Stochastic example
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
  capacity_smolt = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
  kappa = kappa,                  # Productivity in recruits per spawner
  Mjuv_NOS = c(0, -log(SAR), 0),
  fec = c(0, 0, 5040),        # Spawning fecundity of NOS and HOS
  p_female = 0.49,
  s_enroute = 1
  #strays = 0
)

nyears <- 2
HistN <- array(0, c(nsim_stochastic, Bio@maxage, nyears+1, 2))
HistN[, 1, 1, ] <- HistN[, 2, 2, ] <- 1000/SAR

Historical <- new(
  "Historical",
  HistN = HistN
)

png("man/figures/example-kappa.png", height = 4, width = 6, res = 300, units = 'in')
par(mar = c(5, 4, 1, 1))
hist(kappa, main = NULL)
dev.off()

SOM_stochastic <- new("SOM",
                      Bio_stochastic, Habitat, Hatchery, Harvest, Historical,
                      nsim = nsim_stochastic, nyears = 2, proyears = 50)

SMSE_stochastic <- salmonMSE(SOM_stochastic)
saveRDS(SMSE_stochastic, "examples/SMSE_stochastic.rds")

report(SMSE_stochastic, dir = "scratch", filename = "SMSE_stochastic")

png("man/figures/example-PNI-ts.png", height = 4, width = 6, res = 300, units = 'in')
par(mar = c(5, 4, 1, 1))
plot_statevar_ts(SMSE_stochastic, "PNI", quant = TRUE)
dev.off()

png("man/figures/example-PNI-hist.png", height = 4, width = 6, res = 300, units = 'in')
par(mar = c(5, 4, 1, 1))
plot_statevar_hist(SMSE_stochastic, "PNI", y = 49, breaks = 10, xlim = c(0.5, 0.9))
dev.off()

PNI_LT <- SMSE_stochastic@PNI[, 1, 49]
mean(PNI_LT >= 0.8)
quantile(PNI_LT, c(0.025, 0.5, 0.975))

png("man/figures/example-PNI-kappa.png", height = 4, width = 6, res = 300, units = 'in')
par(mar = c(5, 4, 1, 1))
plot(kappa, PNI_LT, xlim = c(1, 7), ylim = c(0.5, 1), panel.first = grid(), pch = 16)
dev.off()

# Extra debugging stuff
SMSE@NOB[, 1, seq(1, 49, 3)]/
  SAHA$NOB[, 1:17]

SMSE@HOB[, 1, seq(1, 49, 3)]/SAHA$HOB[, 1:17]
SAHA$HOB

SMSE@NOB[, 1, ] + SMSE@HOB[, 1, ]
SAHA$NOB + SAHA$HOB

SMSE@NOS[, 1, ]
SAHA$NOS

SMSE@HOS[, 1, ]

SAHA$HOS

SMSE@Egg_HOS[, 1, ]
SAHA$Egg_HOS

SMSE@Egg_NOS[, 1, ]
SAHA$Egg_NOS

SMSE@Smolt_NOS[, 1, ]
SAHA$Smolt_NOS

SMSE@Smolt_HOS[, 1, ]
SAHA$Smolt_HOS

Smolt_all <- SMSE@Smolt_NOS[, 1, ] + SMSE@Smolt_HOS[, 1, ]

SMSE@Return_NOS[, 1, 3, ]
SAHA$Return_NOS

SMSE@Return_HOS[, 1, 3, ]
SAHA$Return_HOS

surv <- SAHA$Return_HOS[, 2:20]/SAHA$Smolt_Rel[, -1]
surv <- SMSE@Return_HOS[, 1, 3, seq(4, 49, 3)]/SMSE@Smolt_Rel[, 1, seq(2, 49, 3)]

SAHA$Smolt_NOS + SAHA$Smolt_HOS
SMSE@Smolt_NOS[, 1, ] + SMSE@Smolt_HOS[, 1, ]

surv <- SAHA$Return_NOS[, 2:20]/(SAHA$Smolt_NOS[, -1] + SAHA$Smolt_HOS[, -1])
surv <- SMSE@Return_NOS[, 1, 3, seq(4, 49, 3)]/(SMSE@Smolt_NOS[, 1, seq(2, 49, 3)] + SMSE@Smolt_HOS[, 1, seq(2, 49, 3)])

plot(surv[1, ], typ = 'o', ylim = c(0.009, 0.01))

exp(-SMSE@Mjuv_loss[, 1, 2, ])[1, ] %>% plot(typ = 'o', ylim = c(0.009, 0.01))
SMSE@fitness[, 1, ]
SAHA$fitness

SOM <- salmonMSE:::check_SOM(SOM)
MOM <- SOM2MOM(SOM)

SAHA$alpha
salmonMSE_env$state %>% reshape2::acast(list("x", "t"), value.var = "alpha")

SAHA$beta
salmonMSE_env$state %>% reshape2::acast(list("x", "t"), value.var = "beta")

