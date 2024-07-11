
# A simple example, same as example 1, with mark-selective fishing (see Harvest object)
library(salmonMSE)

class?SOM # Definition of inputs

SAR <- 0.01
Bio <- new(
  "Bio",
  nsim = 3,
  maxage = 3,
  p_mature = c(0, 0, 1),
  SRrel = "BH",
  capacity_smolt = 17250,     # Beverton-Holt asymptote. Not unfished capacity!!
  kappa = 3,                  # Productivity in recruits per spawner
  Mocean_NOS = c(0, -log(SAR), 0),
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
  m = 0.5,                       # Mark rate of hatchery releases
  release_mort = c(0.1, 0.1),
  vulPT = c(0, 0, 0),
  vulT = c(1, 1, 1)
)

# Return of 1000 natural and hatchery fish each for the first generation
nyears <- 2
HistN <- array(0, c(Bio@nsim, Bio@maxage, nyears, 2))
HistN[, 1, 1, ] <- HistN[, 2, 2, ] <- 1000/SAR

Historical <- new(
  "Historical",
  HistN = HistN
)

# Stitched salmon operating model
SOM <- new("SOM",
           nyears = 2,
           proyears = 50,
           Bio, Habitat, Hatchery, Harvest, Historical)

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

compare(SAHA$Fry_NOS, SMSE@Fry_NOS[, 1, ], "Fry_NOS", ylim = c(0, 150000))

png("man/figures/example-NOS.png", height = 3, width = 7, res = 300, units = 'in')
compare(SAHA$NOS, SMSE@NOS[, 1, ], "NOS", ylim = c(0, 60), yline = 28.9731)
dev.off()



#

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

SMSE@Fry_HOS[, 1, ]
SAHA$Fry_HOS

SMSE@Fry_NOS[, 1, ]
SAHA$Fry_NOS

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

exp(-SMSE@Mocean_loss[, 1, 2, ])[1, ] %>% plot(typ = 'o', ylim = c(0.009, 0.01))
SMSE@fitness[, 1, ]
SAHA$fitness

SOM <- salmonMSE:::check_SOM(SOM)
MOM <- SOM2MOM(SOM)

SAHA$alpha
salmonMSE_env$state %>% reshape2::acast(list("x", "t"), value.var = "alpha")

SAHA$beta
salmonMSE_env$state %>% reshape2::acast(list("x", "t"), value.var = "beta")

