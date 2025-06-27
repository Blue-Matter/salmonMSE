
library(testthat)

source("examples/04_MSOM.R")

test_that("Test MSOM", {
  expect_s4_class(out, "SMSE")

  # Numbers as of June 25, 2025
  expect_equal(round(max(out@p_wild[, 2, ], na.rm = TRUE), 3), 0.444)
  expect_equal(round(max(out@p_wild[, 1, ], na.rm = TRUE), 3), 0.626)
})
