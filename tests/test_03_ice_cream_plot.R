
library(testthat)

source("examples/02_ice_cream_plot.R")

test_that("Test decision table", {
  expect_true(is.list(SMSE_list))

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
    bind_rows() %>%
    round(2)

  pm_test <- readr::read_csv("tests/tables/03_ice_cream_plot_pm_06.26.2025.csv") %>%
    round(2)

  expect_true(all(!abs(pm - pm_test))) # Check decision table is unchanged from June 26, 2025
})
