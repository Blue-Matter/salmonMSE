
kappa <- 3
capacity <- 1000
Smax <- 500
phi <- 1

i <- seq(1, 3, 0.5)

alpha <- kappa/phi

Spawners <- seq(0, 10000, 10)
Smolt_BH_a <- sapply(
  i, function(x) {
    alpha2 <- x * alpha
    beta <- alpha2/capacity
    alpha2 * Spawners / (1 + beta * Spawners)
  }
) %>%
  structure(dimnames = list(Spawners = Spawners, Improve = i)) %>%
  reshape2::melt() %>%
  mutate(SRR = "Beverton-Holt", par = "kappa")

Smolt_BH_b <- sapply(
  i, function(x) {
    beta <- alpha/capacity/x
    alpha * Spawners / (1 + beta * Spawners)
  }
) %>%
  structure(dimnames = list(Spawners = Spawners, Improve = i)) %>%
  reshape2::melt() %>%
  mutate(SRR = "Beverton-Holt", par = "Capacity")

Smolt_R_a <- sapply(
  i, function(x) {
    alpha2 <- x * alpha
    beta <- 1/Smax
    alpha2 * Spawners * exp(-beta * Spawners)
  }
) %>%
  structure(dimnames = list(Spawners = Spawners, Improve = i)) %>%
  reshape2::melt() %>%
  mutate(SRR = "Ricker", par = "kappa")

Smolt_R_b <- sapply(
  i, function(x) {
    beta <- 1/Smax/x
    alpha * Spawners * exp(-beta * Spawners)
  }
) %>%
  structure(dimnames = list(Spawners = Spawners, Improve = i)) %>%
  reshape2::melt() %>%
  mutate(SRR = "Ricker", par = "Capacity")

g <- rbind(
  Smolt_BH_a,
  Smolt_BH_b,
  Smolt_R_a,
  Smolt_R_b
) %>%
  dplyr::filter(Spawners <= 1000) %>%
  ggplot(aes(Spawners, value, colour = factor(Improve))) +
  geom_line() +
  labs(y = "Recruitment", colour = "Improvement") +
  facet_grid(vars(SRR), vars(par)) +
  geom_abline(slope = 1, intercept = 2, linetype = 2) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"))
ggsave("man/figures/SRR.png", g, height = 4, width = 5)


