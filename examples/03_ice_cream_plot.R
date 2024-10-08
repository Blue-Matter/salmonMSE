

Design <- expand.grid(
  kappa = c(3, 6, 9),
  hatch = c(5, 10, 15) * 1000
)

wrapper <- function(x, Design) {

  SAR <- 0.01
  nsim <- 100

  set.seed(100)
  kappa_mean <- Design$kappa[x]
  kappa_sd <- 0.3
  kappa <- rlnorm(nsim, log(kappa_mean) - 0.5 * kappa_sd^2, kappa_sd)

  Bio <- new(
    "Bio",
    nsim = nsim,
    maxage = 3,
    p_mature = c(0, 0, 1),
    SRrel = "BH",
    capacity_smolt = 17250,
    kappa = kappa,
    Mjuv_NOS = c(0, -log(SAR), 0),
    fec = c(0, 0, 5040),
    p_female = 0.49,
    s_enroute = 1
  )

  Hatchery <- new(
    "Hatchery",
    n_yearling = Design$hatch[x],
    n_subyearling = 0,
    s_prespawn = 1,
    s_egg_smolt = 0.92,
    s_egg_subyearling = 1,
    Mjuv_HOS = Bio@Mjuv_NOS,
    gamma = 0.8,
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
    capacity_smolt_improve = 1,
    kappa_improve = 1
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

  # Return of 1000 natural and hatchery fish each for the first generation
  nyears <- 2
  HistN <- array(0, c(Bio@nsim, Bio@maxage, nyears+1, 2))
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

  SMSE <- salmonMSE(SOM)

  return(SMSE)
}

library(snowfall)
sfInit(parallel = TRUE, cpus = nrow(Design))
sfLibrary(salmonMSE)
sfExport("Design")

SMSE_list <- sfLapply(1:nrow(Design), wrapper, Design = Design)

saveRDS(SMSE_list, file = "examples/SMSE_list.rds")
sfStop()

# Make ice cream plot
pm_fn <- function(x, SMSE_list, Design, val = 0.75) {
  out <- Design[x, ]
  out$PNI_75 <- mean(SMSE_list[[x]]@PNI[, 1, 49] >= val)
  out$Catch <- mean(SMSE_list[[x]]@KT_NOS[, 1, 49] + SMSE_list[[x]]@KT_NOS[, 1, 49])
  return(out)
}

pm <- lapply(1:nrow(Design), pm_fn, SMSE_list, Design = Design) %>%
  bind_rows()

g <- plot_decision_table(pm$hatch, pm$kappa, pm$PNI_75, title = "Probability PNI > 0.75",
                         xlab =  "Hatchery releases", ylab = "Compensation ratio (productivity)")
ggsave("man/figures/decision_table_PNI75.png", g, height = 3, width = 3)

# Make tradeoff plot
g <- plot_tradeoff(pm$PNI_75, pm$Catch, factor(pm$kappa), factor(pm$hatch), "PNI_75", "Mean catch",
                   x1lab = "Compensation\nratio", x2lab = "Hatchery\nreleases") +
  scale_shape_manual(values = c(1, 4, 16))
ggsave("man/figures/tradeoff_plot.png", g, height = 3, width = 4.5)

# Make time series
PNI_ts <- lapply(1:nrow(Design), function(x) {
  out <- plot_statevar_ts(SMSE_list[[x]], "PNI", quant = TRUE, figure = FALSE) %>%
    reshape2::melt() %>%
    mutate(kappa = Design$kappa[x], hatch = Design$hatch[x])
  out
}) %>%
  bind_rows() %>%
  rename(Year = Var2) %>%
  dplyr::filter(!is.na(value)) %>%
  tidyr::pivot_wider(names_from = Var1) %>%
  mutate(hatch = paste("Hatchery production", hatch) %>% factor(levels = paste("Hatchery production", c(5, 10, 15) * 1000)))

g <- ggplot(PNI_ts, aes(Year)) +
  geom_line(aes(y = `50%`, colour = factor(kappa))) +
  geom_ribbon(linetype = 3, alpha = 0.1, aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(kappa))) +
  facet_wrap(vars(hatch)) +
  labs(x = "Projection Year", y = "PNI", colour = "Compensation\nratio", fill = "Compensation\nratio") +
  theme(legend.position = "bottom")
ggsave("man/figures/PNI_ts.png", g, height = 3, width = 6)

# Spawners
Design_txt <- Design
Design_txt[, 1] <- paste("Productivity =", Design[, 1])
Design_txt[, 2] <- factor(paste("Hatchery production", Design[, 2]),
                          levels = paste("Hatchery production", c(5000, 10000, 15000)))

g <- compare_spawners(SMSE_list, Design_txt) +
  coord_cartesian(ylim = c(0, 150), expand = FALSE)
ggsave("man/figures/compare_spawners.png", g, height = 5, width = 6)

g <- compare_spawners(SMSE_list, Design_txt, prop = TRUE)

# Fitness
g <- compare_fitness(SMSE_list, Design_txt)
ggsave("man/figures/compare_fitness.png", g, height = 5, width = 6)
