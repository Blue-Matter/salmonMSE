
#' Plot output from salmonMSE
#'
#' @description
#' Various functions that plot the state variables from salmonMSE projections:
#'
#' - `plot_statevar_ts()` produces a time series for all simulations, or with confidence intervals
#' - `plot_statevar_hist()` produces a histogram across all simulations for a particular year
#'
#' @param SMSE Class \linkS4class{SMSE} object returned by [salmonMSE()]
#' @param var Character. Slot for the state variable in `SMSE` object. See `slotNames(SMSE)` for options.
#' @param s Integer. Population index for multi-population model (e.g., `s = 1` is the first population in the model)
#' @param xlab Character. Name of time variable for the figure
#' @param quant Logical, whether to plot individual simulations (FALSE) or the median with 95 percent confidence intervals (TRUE)
#' @param ylab Character. Name of the state variable for the figure
#' @param ylim Vector. Y-axis limits
#' @return Base graphics
#' @importFrom graphics matplot grid
#' @importFrom methods slot
#' @importFrom stats quantile
#' @export
plot_statevar_ts <- function(SMSE, var = "PNI", s = 1, xlab = "Projection Year", quant = FALSE, ylab = var, ylim) {
  x <- slot(SMSE, var)[, s, ]

  if (!quant) {
    xplot <- x
  } else {
    xplot <- apply(x, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
  }

  ind <- colSums(xplot, na.rm = TRUE) > 0
  Year <- 1:SMSE@proyears

  if (missing(ylim)) ylim <- c(0, 1.1) * range(xplot, na.rm = TRUE)

  if (!quant) {
    matplot(Year[ind], t(xplot[, ind]), typ = 'l', col = "grey40", ylim = ylim, lty = 1,
            xlab = "Projection Year", ylab = ylab, panel.first = graphics::grid())
  } else {
    matplot(Year[ind], t(xplot[, ind]), typ = 'o', pch = c(NA, 1, NA), col = 1, lty = c(2, 1, 2), ylim = ylim,
            xlab = "Projection Year", ylab = ylab, panel.first = graphics::grid())
  }
  invisible()
}

#' @name plot_statevar_ts
#' @param y Integer. Projection year for the state variable to plot the histogram
#' @importFrom graphics hist
plot_statevar_hist <- function(SMSE, var = "PNI", s = 1, y, xlab = var, ...) {
  x <- slot(SMSE, var)[, s, y]
  hist(x, xlab = xlab, main = NULL, ...)
  invisible()
}
