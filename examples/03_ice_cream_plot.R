

Design <- expand.grid(
  kappa = c(3, 6, 9),
  hatch = c(0, 5, 10, 15) * 1000
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
    maxage = 3,
    p_mature = c(0, 0, 1),
    SRrel = "BH",
    capacity = 17250,
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
    m = 0,
    pmax_esc = 0.50,
    pmax_NOB = 1,
    ptarget_NOB = 0.75,
    #m = 1,
    #pmax_esc = 1,
    #pmax_NOB = 0.7,
    #ptarget_NOB = 0.51,
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
    u_terminal = 0.5,              # Specify fixed harvest rate of mature fish
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

  return(SMSE)
}

library(snowfall)
sfInit(parallel = TRUE, cpus = min(0.5 * parallel::detectCores(), nrow(Design)))
sfLibrary(salmonMSE)
sfExport("Design")

SMSE_list <- sfLapply(1:nrow(Design), wrapper, Design = Design)

saveRDS(SMSE_list, file = "examples/SMSE/SMSE_list.rds")
sfStop()

if (FALSE) {

  # Make ice cream plot (PNI, catch, SMSY)
  pm_fn <- function(x, SMSE_list, Design) {
    out <- Design[x, ]
    out$PNI <- mean(SMSE_list[[x]]@PNI[, 1, 49])
    out$PNI_50 <- PNI50(SMSE_list[[x]], Yrs = c(49, 49))
    out$PNI_80 <- PNI80(SMSE_list[[x]], Yrs = c(49, 49))

    KNOS <- SMSE_list[[x]]@KT_NOS[, 1, 49] # Catch of natural fish
    KHOS <- SMSE_list[[x]]@KT_HOS[, 1, 49] # Catch of hatchery fish

    out$Catch <- mean(KNOS + KHOS)
    out$Catch60 <- mean((KNOS + KHOS) >= 60)
    out$`S/SMSY` <- SMSY85(SMSE_list[[x]], Yrs = c(49, 49))
    return(out)
  }

  pm <- lapply(1:nrow(Design), pm_fn, SMSE_list, Design = Design) %>%
    bind_rows()
  readr::write_csv(pm, "tests/tables/03_ice_cream_plot_pm_06.26.2025.csv")

  g <- plot_decision_table(pm$hatch, pm$kappa, pm$PNI_80, title = "Probability PNI > 0.80",
                           xlab =  "Hatchery releases", ylab = "Productivity")
  ggsave("man/figures/decision_table_PNI80.png", g, height = 3, width = 3)

  g <- plot_decision_table(pm$hatch, pm$kappa, pm$Catch60, title = "Probability Catch > 60",
                           xlab =  "Hatchery releases", ylab = "Productivity")
  ggsave("man/figures/decision_table_catch60.png", g, height = 3, width = 3)

  g <- plot_decision_table(pm$hatch, pm$kappa, pm$`S/SMSY`, title = "Probability NOS > SMSY",
                           xlab =  "Hatchery releases", ylab = "Productivity")
  #ggsave("man/figures/decision_table_SMSY.png", g, height = 3, width = 3)

  # Make tradeoff plot
  g <- plot_tradeoff(pm$PNI_80, pm$Catch60, factor(pm$kappa), factor(pm$hatch), "PNI_80", "Catch60",
                     x1lab = "Productivity", x2lab = "Hatchery\nreleases") +
    scale_shape_manual(values = c(1, 2, 4, 16))
  #ggsave("man/figures/tradeoff_plot_pm.png", g, height = 3, width = 4.5)

  g <- plot_tradeoff(pm$PNI, pm$Catch, factor(pm$kappa), factor(pm$hatch), "Mean PNI", "Mean Catch",
                     x1lab = "Productivity", x2lab = "Hatchery\nreleases") +
    scale_shape_manual(values = c(1, 2, 4, 16)) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = 2)
  ggsave("man/figures/tradeoff_plot_mean.png", g, height = 3, width = 4.5)

  # Make tradeoff figure with median and confidence intervals
  PNI_fn <- function(x, SMSE_list, Design) {
    out <- Design[x, ]

    val <- quantile(SMSE_list[[x]]@PNI[, 1, 49], c(0.025, 0.5, 0.975))

    out$lower <- val[1]
    out$median <- val[2]
    out$upper <- val[3]
    return(out)
  }
  PNI <- lapply(1:nrow(Design), PNI_fn, SMSE_list, Design = Design) %>%
    bind_rows()

  # Next calculate the median and bounds for catch for each scenario
  Catch_fn <- function(x, SMSE_list, Design) {
    out <- Design[x, ]

    KNOS <- SMSE_list[[x]]@KT_NOS[, 1, 49] # Catch of natural fish
    KHOS <- SMSE_list[[x]]@KT_HOS[, 1, 49] # Catch of hatchery fish
    val <- quantile(KNOS + KHOS, c(0.025, 0.5, 0.975))

    out$lower <- val[1]
    out$median <- val[2]
    out$upper <- val[3]
    return(out)
  }
  Catch <- lapply(1:nrow(Design), Catch_fn, SMSE_list, Design = Design) %>%
    bind_rows()

  # Provide the matrix of PNI and Catch to plot_tradeoff()
  g <- plot_tradeoff(as.matrix(PNI[, 3:5]), as.matrix(Catch[, 3:5]),
                     factor(PNI$kappa), factor(PNI$hatch), "PNI", "Catch",
                     x1lab = "Productivity", x2lab = "Hatchery\nreleases") +
    geom_vline(xintercept = c(0.5, 0.8), linetype = 2) +
    scale_shape_manual(values = c(1, 2, 4, 16))
  ggsave("man/figures/tradeoff_plot_median_ci.png", g, height = 3, width = 4.5)

  # Kobe figure
  Kobe_fn <- function(x, SMSE_list, Design) {
    out <- Design[x, ]

    Sp <- rowSums(SMSE_list[[x]]@NOS[, 1, , 49] + SMSE_list[[x]]@HOS[, 1, , 49])
    SMSY <- SMSE_list[[x]]@Misc$Ref[[1]]["Spawners_MSY", ]

    KNOS <- SMSE_list[[x]]@KT_NOS[, 1, 49] # Catch of natural fish
    KHOS <- SMSE_list[[x]]@KT_HOS[, 1, 49] # Catch of hatchery fish

    RNOS <- rowSums(SMSE_list[[x]]@Return_NOS[, 1, , 49])
    RHOS <- rowSums(SMSE_list[[x]]@Return_HOS[, 1, , 49])

    U <- (KNOS + KHOS)/(RNOS + RHOS)
    UMSY <- SMSE_list[[x]]@Misc$Ref[[1]]["UT_MSY", ]

    relS <- quantile(Sp/SMSY, c(0.025, 0.5, 0.975))
    relU <- quantile(U/UMSY, c(0.025, 0.5, 0.975))

    cbind(out, data.frame(t(relS)), data.frame(t(relU)))
    #out$lower <- val[1]
    #out$median <- val[2]
    #out$upper <- val[3]
    #return(out)
  }
  Kobe <- lapply(1:nrow(Design), Kobe_fn, SMSE_list, Design = Design) %>%
    bind_rows()


  g <- plot_tradeoff(as.matrix(Kobe[, 3:5]), as.matrix(Kobe[, 6:8]),
                     factor(PNI$kappa), factor(PNI$hatch), expression(S/S[MSY]), expression(U/U[MSY]),
                     x1lab = "Productivity", x2lab = "Hatchery\nreleases") +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_shape_manual(values = c(1, 2, 4, 16)) +
    expand_limits(x = 0, y = 0)
  #ggsave("man/figures/Kobe_plot_median_ci.png", g, height = 3, width = 4.5)

  # Make time series
  PNI_ts <- lapply(1:nrow(Design), function(x) {
    plot_statevar_ts(SMSE_list[[x]], "PNI", quant = TRUE, figure = FALSE) %>%
      reshape2::melt() %>%
      mutate(kappa = Design$kappa[x], hatch = Design$hatch[x])
  }) %>%
    bind_rows() %>%
    rename(Year = Var2) %>%
    dplyr::filter(!is.na(value)) %>%
    tidyr::pivot_wider(names_from = Var1) %>%
    mutate(hatch = paste("Hatchery production", hatch) %>% factor(levels = paste("Hatchery production", c(0, 5, 10, 15) * 1000)))

  g <- ggplot(PNI_ts, aes(Year)) +
    geom_line(aes(y = `50%`, colour = factor(kappa))) +
    geom_ribbon(linetype = 3, alpha = 0.1, aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(kappa))) +
    facet_wrap(vars(hatch)) +
    labs(x = "Projection Year", y = "PNI", colour = "Compensation\nratio", fill = "Compensation\nratio") +
    theme(legend.position = "bottom")
  ggsave("man/figures/PNI_ts.png", g, height = 4, width = 5)

  # Spawners
  Design_txt <- Design
  Design_txt[, 1] <- factor(paste("Productivity =", Design[, 1]),
                            levels = paste("Productivity =", c(9, 6, 3)))
  Design_txt[, 2] <- factor(paste("Hatchery production", Design[, 2]),
                            levels = paste("Hatchery production", c(0, 5000, 10000, 15000)))

  g <- compare_spawners(SMSE_list, Design_txt) +
    coord_cartesian(ylim = c(0, 75), expand = FALSE)
  ggsave("man/figures/compare_spawners.png", g, height = 5, width = 7)

  g <- compare_spawners(SMSE_list, Design_txt, prop = TRUE)

  # Fitness
  g <- compare_fitness(SMSE_list, Design_txt)
  ggsave("man/figures/compare_fitness.png", g, height = 5, width = 7)

  # Escapement
  g <- compare_escapement(SMSE_list, Design_txt)
  #ggsave("man/figures/compare_escapement.png", g, height = 5, width = 6)
}
