
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
#' @return Functions return the matrix of plotted values invisibly. Figure plotted from base graphics
#' @importFrom graphics matplot grid
#' @importFrom methods slot
#' @importFrom stats quantile
#' @seealso [plot_decision_table()]
#' @export
plot_statevar_ts <- function(SMSE, var = "PNI", s = 1, xlab = "Projection Year", quant = FALSE, ylab = var, ylim, ...) {
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
    matplot(Year[ind], t(xplot[, ind]), type = 'l', col = "grey40", ylim = ylim, lty = 1,
            xlab = "Projection Year", ylab = ylab, panel.first = graphics::grid(), ...)
  } else {
    matplot(Year[ind], t(xplot[, ind]), type = 'o', pch = c(NA, 1, NA), col = 1, lty = c(2, 1, 2), ylim = ylim,
            xlab = "Projection Year", ylab = ylab, panel.first = graphics::grid(), ...)
  }
  invisible(xplot)
}

#' @name plot_statevar_ts
#' @param y Integer. Projection year for the state variable to plot the histogram
#' @param ... Additional arguments to base plot function
#' @importFrom graphics hist
plot_statevar_hist <- function(SMSE, var = "PNI", s = 1, y, xlab = var, ...) {
  x <- slot(SMSE, var)[, s, y]
  hist(x, xlab = xlab, main = NULL, ...)
  invisible(x)
}

#' Decision table of performance metrics
#'
#' Generates a coloured table of a performance metric across two axes, which may be a population dynamics variable
#' (e.g., productivity) or a management action (e.g., hatchery production levels or harvest strategy).
#' See example at \url{https://docs.salmonmse.com/articles/decision-table.html}
#'
#' @param x Atomic, vector of values for the x axis (same length as z). Will be converted to factors
#' @param y Atomic, vector of values for the y axis (same length as z). Will be converted to factors
#' @param z Numeric, vector of values for the performance metric
#' @param title Character, optional title of figure
#' @param xlab Character, optional x-axis label
#' @param ylab Character, optional y-axis label
#' @return ggplot object
#' @seealso [plot_statevar_ts()] [plot_tradeoff()]
#' @import ggplot2
#' @export
plot_decision_table <- function(x, y, z, title, xlab, ylab) {
  dt <- data.frame(x = x, y = y, z = z)
  dt$txt <- txt <- format(round(z, 2))

  g <- ggplot(dt, aes(factor(x), factor(y), fill = z)) +
    geom_tile(colour = "grey20", alpha = 0.6) +
    geom_text(aes(label = txt)) +
    coord_cartesian(expand = FALSE) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "deeppink", high = "green4", mid = "white", midpoint = 0.5)

  if (!missing(xlab)) g <- g + labs(x = xlab)
  if (!missing(ylab)) g <- g + labs(y = ylab)
  if (!missing(title)) g <- g + ggtitle(title)

  return(g)
}

#' Tradeoff figure
#'
#' Generates a tradeoff figure, a comparison between two performance metrics, across two variables which may represent
#' a population dynamics variable (e.g., productivity) or a management action (e.g., hatchery production levels or harvest strategy).
#' See example at \url{https://docs.salmonmse.com/articles/decision-table.html}
#'
#' @param pm1 Numeric, vector of values for the first performance metric on the x-axis
#' @param pm2 Numeric, vector of values for the second performance metric on the y-axis (same length as pm1)
#' @param x1 Atomic, vector of values for the first grouping variable. Various levels are represented by colours. Same length as pm1.
#' @param x2 Numeric, vector of values for the second grouping variable. Various levels are represented by shapes. Same length as pm1.
#' @param xlab Character, optional x-axis label
#' @param ylab Character, optional y-axis label
#' @param x1lab Character, optional label for the first grouping variable
#' @param x2lab Character, optional label for the second grouping variable
#' @return ggplot object
#' @seealso [plot_statevar_ts()] [plot_decision_table()]
#' @import ggplot2
#' @export
plot_tradeoff <- function(pm1, pm2, x1, x2, xlab, ylab, x1lab, x2lab) {

  if (missing(x1)) x1 <- 0
  if (missing(x2)) x2 <- 1

  dt <- data.frame(
    pm1 = pm1,
    pm2 = pm2,
    x1 = x1,
    x2 = x2
  )

  g <- ggplot(dt, aes(pm1, pm2, colour = x1, shape = x2)) +
    geom_point() +
    theme_bw()

  if (length(x1) == 1) {
    g <- g +
      scale_colour_manual(values = GeomPoint$default_aes$colour) +
      guides(colour = "none")
  }
  if (length(x2) == 2) {
    g <- g +
      scale_shape_manual(values = GeomPoint$default_aes$shape) +
      guides(shape = "none")
  }

  if (!missing(xlab)) g <- g + labs(x = xlab)
  if (!missing(ylab)) g <- g + labs(y = ylab)
  if (!missing(x1lab)) g <- g + labs(colour = x1lab)
  if (!missing(x2lab)) g <- g + labs(shape = x2lab)

  return(g)
}

