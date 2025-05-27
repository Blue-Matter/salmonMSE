
library(salmonMSE)

# See Appendix A of Mahony et al. 2021. CSAS Research Document 2021/009.
# https://publications.gc.ca/site/eng/9.897509/publication.html

# Source code for original analysis is at: https://github.com/wchallenger/OK_CNPVA
# Use parameters in the markdown document: https://github.com/wchallenger/OK_CNPVA/blob/master/PVA-analysis.nb.html

maxage <- 5
nsim <- 100
nyears <- 2
proyears <- 45 # 2020 to 2065

Mjuv_NOS <- Mjuv_HOS <- array(0, c(nsim, maxage, nyears + proyears))

# For the PVA, we want to model:
# Deviations in the smolt production function to get to age 1
# To get to age 2, we need early juvenile survival (phi_ocean, stochastic), dam survival (fixed) and then SU_2 (fixed)
# Currently salmonMSE doesn't have Ricker deviations so we include it in the age 1 mortality rate
set.seed(1)
Ricker_sd <- 0.53
Ricker_dev <- exp(rnorm(nsim * proyears, 0, Ricker_sd)) |> matrix(nsim, proyears)

phi_ocean <- 0.084

# https://github.com/wchallenger/OK_CNPVA/blob/master/CNPVA/R/PopSim.fn.R#L41
D_juv <- rnorm(nsim * proyears, 0.35, 0.09) |> matrix(nsim, proyears)
D_juv[D_juv < 0] <- 0.01
D_juv[D_juv > 1] <- 1

SU_2 <- 1

# Age 1 survival
surv1 <- phi_ocean * D_juv * SU_2
surv24 <- c(0.7, 0.8, 0.9)  # SU_3, SU_4, SU_5

Mjuv_NOS[, 1, nyears + seq(1, proyears)] <- -log(surv1 * Ricker_dev)
Mjuv_NOS[, 1, seq(1, nyears)] <- array(Mjuv_NOS[, , nyears + 1], c(nsim, 1, nyears))

Mjuv_HOS[, 1, nyears + seq(1, proyears)] <- -log(surv1)
Mjuv_HOS[, 1, seq(1, nyears)] <- array(Mjuv_HOS[, , nyears + 1], c(nsim, 1, nyears))

# Age 2-4 survival
Mjuv_NOS[, 2:4, ] <- Mjuv_HOS[, 2:4, ] <- array(-log(surv24), c(3, nsim, nyears + proyears)) |> aperm(c(2, 1, 3))

# Arbitrary juvenile M for age 5
Mjuv_NOS[, 5, ] <- Mjuv_HOS[, 5, ] <- 0.01

Bio <- new(
  "Bio",
  maxage = maxage,
  p_mature = c(0, 0.04, 0.26, 0.72, 1),
  SRrel = "Ricker",
  kappa = 136, # set phi = 1 so that kappa = alpha = 136
  phi = 1,
  Smax = 2400,
  Mjuv_NOS = Mjuv_NOS,
  p_female = 1,         # fec = 1, p_female = 1 so egg production = number of spawners
  fec = rep(1, maxage), # fec = 1, p_female = 1 so egg production = number of spawners
  s_enroute = 0.94      # D_adult, not stochastic in salmonMSE
)

Harvest <- new(
  "Harvest",
  u_preterminal = 0.25,  # HR_ocean
  u_terminal = 0.42,     # HR_river
  vulPT = c(0, 1, 1, 1, 1),
  vulT = c(0, 1, 1, 1, 1)
)

Hatchery_base <- new(
  "Hatchery",
  n_yearling = 0,
  n_subyearling = 0
)

Hatchery_500k <- new(
  "Hatchery",
  n_yearling = 5e5,
  n_subyearling = 0,
  s_prespawn = 1,
  s_egg_smolt = 1,
  s_egg_subyearling = 1,
  brood_import = c(0, 0, 0, 0, 5e6),
  Mjuv_HOS = Mjuv_HOS,
  gamma = 1,
  m = 0,
  pmax_esc = 0,
  pmax_NOB = 0,
  ptarget_NOB = 0,
  phatchery = 0,
  premove_HOS = 0,
  fec_brood = Bio@fec,
  fitness_type = c("none", "none")
)

Habitat <- new(
  "Habitat",
  use_habitat = FALSE
)

# Calculate initial equilibrium juvenile abundance at age such that the number of spawners = 50
Nsp <- sapply(1:nsim, function(x) {
  nS <- 50
  # Per-smolt calculations
  #Fjuv <- -log(1 - Harvest@u_preterminal) * Harvest@vulPT
  #Fterm <- -log(1 - Harvest@u_terminal) * Harvest@vulT
  #juv_surv <- salmonMSE:::calc_survival(Mjuv_NOS[x, , nyears + 1] + Fjuv, Bio@p_mature)
  #esc <- juv_surv * Bio@p_mature * exp(-Fterm)

  juv_surv <- salmonMSE:::calc_survival(Mjuv_NOS[x, , nyears + 1], Bio@p_mature)
  esc <- juv_surv * Bio@p_mature
  spawners <- esc * Bio@s_enroute

  R <- nS/sum(spawners)
  Njuv <- R * juv_surv
  Nesc <- R * esc
  Nsp <- R * spawners
  Njuv
})

Njuv <- array(Nsp, c(maxage, nsim, nyears + 1)) |> aperm(c(2, 1, 3))

Historical <- new(
  "Historical",
  HistNjuv_NOS = Njuv,
  HistNjuv_HOS = array(0, c(nsim, maxage, nyears + 1))
)

SOM <- new("SOM",
           Name = "Okanagan Chinook PVA",
           nsim = nsim,
           nyears = nyears,
           proyears = proyears,
           seed = 1,
           Bio = Bio,
           Habitat = Habitat,
           Hatchery = Hatchery_base,
           Harvest = Harvest,
           Historical = Historical)

SMSE <- salmonMSE(SOM)

saveRDS(SMSE, file = "examples/Okanagan_PVA.rds")
report(SMSE, dir = "examples", filename = "Okanagan_PVA")

SOM_500k <- new("SOM",
                Name = "Okanagan Chinook PVA with 500k hatchery",
                nsim = nsim,
                nyears = nyears,
                proyears = proyears,
                seed = 1,
                Bio = Bio,
                Habitat = Habitat,
                Hatchery = Hatchery_500k,
                Harvest = Harvest,
                Historical = Historical)

SMSE_500k <- salmonMSE(SOM_500k)

saveRDS(SMSE_500k, file = "examples/Okanagan_PVA_500k.rds")
report(SMSE_500k, dir = "examples", filename = "Okanagan_PVA_500k")


# Plot distribution of spawners
Sp <- data.frame(
  NOS = apply(SMSE_500k@NOS[, 1, , proyears - 1], c(1, 3), sum),
  HOS = apply(SMSE_500k@HOS[, 1, , proyears - 1], c(1, 3), sum)
) %>%
  mutate(Total = NOS + HOS) %>%
  reshape2::melt()

quants <- Sp %>% dplyr::filter(variable == "Total") %>% pull(value) %>%
  quantile(c(0.025, 0.5, 0.975))

g <- ggplot(Sp, aes(value, colour = variable, fill = variable)) +
  geom_density(alpha = 0.5) +
  expand_limits(x = 0) +
  geom_vline(xintercept = quants["50%"], linetype = 2) +
  geom_vline(xintercept = quants[c("2.5%", "97.5%")], linetype = 3) +
  labs(x = "Spawners", y = "Density", colour = NULL, fill = NULL) +
  ggtitle("Long-term spawners with 500k hatchery releases")
#ggsave("examples/PVA_spawners.png", g, height = 3.5, width = 6)

# Plot annual spawners
Sp <- rbind(
  apply(SMSE_500k@NOS[, 1, , ], c(1, 3), sum) %>% reshape2::melt() %>% mutate(type = "NOS"),
  apply(SMSE_500k@HOS[, 1, , ], c(1, 3), sum) %>% reshape2::melt() %>% mutate(type = "HOS")
) %>%
  rename(Sim = Var1, Year = Var2) %>%
  dplyr::filter(Year < proyears)
Total <- Sp %>%
  summarise(value = sum(value), .by = c(Sim, Year)) %>%
  mutate(type = "Total")

g <- Total %>%
  summarise(p = mean(value > 1000), .by = Year) %>%
  ggplot(aes(Year, p)) +
  geom_line() +
  geom_point()

g <- rbind(Sp, Total) %>%
  summarise(med = median(value),
            min = quantile(value, 0.025),
            max = quantile(value, 0.975),
            .by = c(Year, type)) %>%
  mutate(type = factor(type, levels = c("NOS", "HOS", "Total"))) %>%
  ggplot(aes(Year, med, colour = type)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = min, ymax = max, fill = type), alpha = 0.25) +
  labs(x = "Year", y = "Spawners", fill = NULL, colour = NULL)
#ggsave("examples/PVA_spawners_ts.png", g, height = 3.5, width = 6)



