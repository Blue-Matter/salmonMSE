
library(testthat)

source("examples/01_simple.R")

test_that("Simple SMSE", {
  expect_s4_class(SMSE, "SMSE")
  expect_true(is.list(SAHA))
})

test_that("Simple SMSE with MSF", {
  u_terminal <- SMSE_MSF@Misc$SOM@Harvest[[1]]@u_terminal

  expect_s4_class(SMSE_MSF, "SMSE")
  expect_lt(max(SMSE_MSF@KT_NOS), 1e-8)  # No natural origin retained catch
  expect_true(round(max(SMSE_MSF@ExT_HOS), 3) == u_terminal)  # HO ex. rate = Harvest rate because m = 1

  m <- SMSE_MSF@Misc$SOM@Hatchery[[1]]@m
  delta <- SMSE_MSF@Misc$SOM@Harvest[[1]]@release_mort[2]

  F_kept <- -log(1 - max(SMSE_MSF@ExT_HOS))
  E <- F_kept/m
  F_NO <- delta * E
  Ex_NO <- 1 - exp(-F_NO)

  expect_gt(max(SMSE_MSF@ExT_NOS), 0)  # NO ex. rate > 0
  expect_equal(round(max(SMSE_MSF@ExT_NOS), 4), round(Ex_NO, 4))
})




test_that("Stochastic SMSE", {
  expect_s4_class(SMSE_stochastic, "SMSE")
})

test_that("Stochastic SMSE: P[PNI > 0.8] = 0.13?", {
  expect_true({
    PNI_LT <- SMSE_stochastic@PNI[, 1, 48]
    mean(PNI_LT >= 0.8) == 0.13
  })
})
