
#' Plot core output from salmonMSE
#'
#' @description
#' Various functions that plot the state variables from salmonMSE projections:
#'
#' - `plot_statevar_ts()` produces a time series for all simulations, or with medians and 95th percentile intervals
#' - `plot_statevar_hist()` produces a histogram across all simulations for a particular year
#' - `plot_spawners()` produces a summary barplot of spawners, including NOS, HOS, and wild spawners
#' - `plot_escapement()` produces a summary figure of the proportion of spawners and broodtake to escapement
#' - `plot_fitness()` produces a summary figure of metrics (fitness, PNI, pHOS, and pWILD) related to hatchery production
#' - `plot_fishery()` produces a summary figure of metrics related to the fishery, e.g., median catch, exploitation rate or harvest rate
#'
#' @param SMSE Class \linkS4class{SMSE} object returned by [salmonMSE()]
#' @param var Character. Slot for the state variable in `SMSE` object. See `slotNames(SMSE)` for options. Additional supported options are:
#' `"ESS"` (egg-smolt survival), `"pbrood"` (broodtake to escapement ratio), `"pNOSesc"` (NOS/natural escapement), `"pHOSesc"` (HOS/hatchery escapement),
#' `Total Spawners` (NOS + HOS), `NOS/SMSY`, `S/SMSY`, and `NOS/Sgen`.
#' @param s Integer. Population index for multi-population model (e.g., `s = 1` is the first population in the model)
#' @param xlab Character. Name of time variable for the figure
#' @param figure Logical, whether to generate a figure (set to FALSE if only using the function to return the data matrix)
#' @param quant Logical, whether to plot individual simulations (FALSE) or the median with 95 percent confidence intervals (TRUE)
#' @param ylab Character. Name of the state variable for the figure
#' @param ylim Vector. Y-axis limits
#' @param agg.fun Function. Defines how to aggregate state variables that are reported by age. Typically, `sum` is used but `max` is also
#' possible for reporting apical exploitation rates.
#' @return Functions return the matrix of plotted values invisibly. Figure plotted from base graphics
#' @importFrom graphics matplot grid legend
#' @importFrom methods slot
#' @importFrom stats quantile
#' @importFrom grDevices hcl.colors
#' @seealso [plot_decision_table()] [plot_LHG()] [compare_statevar_ts()]
#' @export
plot_statevar_ts <- function(SMSE, var = "PNI", s = 1, figure = TRUE, xlab = "Projection Year", quant = FALSE, ylab = var, ylim,
                             agg.fun = sum, ...) {

  if (length(s) > 1) {
    x <- sapply(s, function(i) get_statevar(SMSE, var, i, agg.fun = agg.fun), simplify = "array")
    xplot <- apply(x, 2:3, quantile, 0.5, na.rm = TRUE)

    if (figure && any(xplot > 0, na.rm = TRUE)) {
      ind <- rowSums(xplot, na.rm = TRUE) > 0
      Year <- 1:SMSE@proyears

      if (missing(ylim)) ylim <- c(0, 1.1) * range(xplot, na.rm = TRUE)

      col <- grDevices::hcl.colors(length(s), palette = "Dark 3", alpha = 1)

      matplot(Year[ind], xplot[ind, ], type = 'o', ylim = ylim, col = col, lty = 1, pch = 1,
              xlab = "Projection Year", ylab = ylab, ...)
      legend("topleft", legend = SMSE@Snames, col = col, lty = 1, pch = 1, bty = "n")

    }

  } else {

    x <- get_statevar(SMSE, var, s, agg.fun = agg.fun)

    if (!quant) {
      xplot <- x
    } else {
      xplot <- apply(x, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
    }

    if (figure && any(xplot > 0, na.rm = TRUE)) {
      ind <- colSums(xplot, na.rm = TRUE) > 0
      Year <- 1:SMSE@proyears

      if (missing(ylim)) ylim <- c(0, 1.1) * range(xplot, na.rm = TRUE)

      if (!quant) {
        matplot(Year[ind], t(xplot[, ind]), type = 'l', col = "grey40", ylim = ylim, lty = 1,
                xlab = "Projection Year", ylab = ylab, ...)
      } else {
        matplot(Year[ind], t(xplot[, ind]), type = 'o', pch = c(NA, 1, NA), col = 1, lty = c(2, 1, 2), ylim = ylim,
                xlab = "Projection Year", ylab = ylab, ...)
      }
    }
  }


  invisible(xplot)
}

#' @name plot_statevar_ts
#' @param y Integer. Projection year for the state variable to plot the histogram. If missing, the last projection year is used.
#' @param ... Additional arguments to base plot function
#' @export
#' @importFrom graphics hist
plot_statevar_hist <- function(SMSE, var = "PNI", s = 1, y, figure = TRUE, xlab = var, ...) {
  x <- get_statevar(SMSE, var, s)

  if (missing(y)) {
    xvar <- colSums(x, na.rm = TRUE)
    if (all(!xvar)) {
      y <- length(xvar)
    } else {
      ymax <- max(which(!is.na(xvar) & xvar > 0))
      y <- ymax
    }
  }
  xplot <- x[, y]
  if (figure && any(xplot > 0, na.rm = TRUE)) hist(xplot, xlab = xlab, main = NULL, ...)
  invisible(xplot)
}

get_statevar <- function(SMSE, var, s, agg.fun = sum) {

  x <- matrix(NA, SMSE@nsim, SMSE@proyears)

  if (var %in% slotNames(SMSE)) {

    .x <- slot(SMSE, var)
    if (length(dim(.x)) == 4) {
      x[] <- apply(.x[, s, , ], c(1, 3), agg.fun)
    } else if (length(dim(.x)) == 3) {
      x[] <- .x[, s, ]
    } else {
      stop("Dimension of slot ", var, " is not 3 or 4")
    }
    return(x)
  }

  if (var == "ESS") {
    x[, seq(2, SMSE@proyears)] <- local({
      Smolt <- SMSE@Smolt_NOS[, s, ] + SMSE@Smolt_HOS[, s, ]
      Egg <- SMSE@Egg_NOS[, s, ] + SMSE@Egg_HOS[, s, ]

      Smolt[, seq(2, SMSE@proyears)]/Egg[, seq(2, SMSE@proyears) - 1]
    })
  }

  if (var == "pbrood") {
    x[] <- local({
      Esc_NOS <- apply(SMSE@Escapement_NOS[, s, , ], c(1, 3), sum)
      Esc_HOS <- apply(SMSE@Escapement_HOS[, s, , ], c(1, 3), sum)
      (SMSE@NOB[, s, ] + SMSE@HOB[, s, ])/(Esc_NOS + Esc_HOS)
    })
  }

  if (var == "pNOSesc") {
    x[] <- local({
      Esc_NOS <- apply(SMSE@Escapement_NOS[, s, , ], c(1, 3), sum)
      NOS <- apply(SMSE@NOS[, s, , ], c(1, 3), sum)
      NOS/Esc_NOS
    })
  }
  if (var == "pHOSesc") {
    x[] <- local({
      Esc_HOS <- apply(SMSE@Escapement_HOS[, s, , ], c(1, 3), sum)
      HOS <- apply(SMSE@HOS[, s, , ], c(1, 3), sum)
      HOS/Esc_HOS
    })
  }

  if (var == "Smolt") x[] <- SMSE@Smolt_NOS[, s, ] + SMSE@Smolt_HOS[, s, ]

  if (var == "Total Spawners") x[] <- apply(slot(SMSE, "NOS")[, s, , ] + slot(SMSE, "HOS")[, s, , ], c(1, 3), sum)

  if (var == "NOS/SMSY" && length(SMSE@Misc$Ref[[s]])) {
    x[] <- local({
      NOS <- apply(slot(SMSE, "NOS")[, s, , ], c(1, 3), sum)
      SMSY <- SMSE@Misc$Ref[[s]]["Spawners_MSY", ]
      NOS/SMSY
    })
  }

  if (var == "S/SMSY" && length(SMSE@Misc$Ref[[s]])) {
    x[] <- local({
      S <- apply(slot(SMSE, "NOS")[, s, , ] + slot(SMSE, "HOS")[, s, , ], c(1, 3), sum)
      SMSY <- SMSE@Misc$Ref[[s]]["Spawners_MSY", ]
      S/SMSY
    })
  }

  if (var == "NOS/Sgen" && length(SMSE@Misc$Ref[[s]])) {
    x[] <- local({
      NOS <- apply(slot(SMSE, "NOS")[, s, , ], c(1, 3), sum)
      Sgen <- SMSE@Misc$Ref[[s]]["Sgen", ]
      NOS/Sgen
    })
  }

  x[is.infinite(x)] <- NA

  return(x)
}

#' @name plot_statevar_ts
#' @param prop Logical, whether to plot proportions or absolute numbers
#' @param FUN Summarizing function across simulations, typically [median()] or [mean()]
#' @importFrom graphics barplot box
#' @importFrom stats median
#' @export
plot_spawners <- function(SMSE, s = 1, prop = TRUE, FUN = median, figure = TRUE, ylim) {
  Year <- 1:SMSE@proyears

  .HOS <- apply(SMSE@HOS[, s, , ], c(1, 3), sum)
  HOS <- apply(.HOS, 2, FUN)

  .NOS <- apply(SMSE@NOS[, s, , ], c(1, 3), sum)
  NOS <- apply(.NOS, 2, FUN)

  p_wild <- SMSE@p_wild[, s, ]
  p_wild[is.na(p_wild)] <- 0

  Spawners <- HOS + NOS
  WILD <- apply(p_wild * .NOS, 2, FUN)
  NOS_notWILD <- apply((1 - p_wild) * .NOS, 2, FUN)

  x <- rbind(WILD, NOS_notWILD, HOS)

  if (prop) {
    x <- apply(x, 2, function(i) i/sum(i))
    x[is.na(x)] <- 0
  }

  if (figure) {
    if (missing(ylim)) {
      if (prop) {
        ylim <- c(0, 1)
      } else {
        ylim <- c(0, 1.1) * range(Spawners)
      }
    }
    col <- c("#004533", "#76A6D0", "#FFF7FD") #grDevices::hcl.colors(3, palette = "PuBuGn")
    plot(Year, Spawners, xlim = range(Year) + c(-1, 0),
         xlab = "Projection Year", ylab = ifelse(prop, "Proportion", "Spawners"),
         type = "n", ylim = ylim, xaxs = "i", yaxs = "i")
    barplot(x, legend.text = rownames(x), space = 0, xlim = range(Year),
            col = col, border = "grey40", add = TRUE, xpd = FALSE)
    box()
  }

  invisible(x)
}

#' @name plot_statevar_ts
#' @importFrom graphics legend lines
#' @export
plot_fitness <- function(SMSE, s = 1, FUN = median, figure = TRUE, ylim) {
  Year <- 1:SMSE@proyears

  Fitness <- apply(SMSE@fitness[, s, 1, ], 2, FUN)
  PNI <- apply(SMSE@PNI[, s, ], 2, FUN)
  pHOSeff <- apply(SMSE@pHOS_effective[, s, ], 2, FUN)
  pWILD <- apply(SMSE@p_wild[, s, ], 2, FUN)

  x <- cbind(Fitness, PNI, pHOSeff, pWILD)

  if (figure) {
    if (missing(ylim)) ylim <- c(0, 1)
    matplot(Year, x, type = "n", xlab = "Projection Year", ylab = "Value",
            ylim = ylim)

    col <- 1:ncol(x)
    pch <- c(1, 4, 16, 18)
    for (i in 1:ncol(x)) {
      x_i <- x[, i]
      lines(Year[!is.na(x_i)], x_i[!is.na(x_i)], type = "o", col = col[i], pch = pch[i])
    }
    legend("bottomleft", legend = colnames(x), col = col, pch = pch, lwd = 1, bty = "n")
  }

  invisible(x)
}

#' @name plot_statevar_ts
#' @export
plot_escapement <- function(SMSE, s = 1, FUN = median, figure = TRUE, ylim) {

  Year <- 1:SMSE@proyears

  pNOSesc <- plot_statevar_ts(SMSE, "pNOSesc", s, figure = FALSE) %>%
    apply(2, FUN)
  pHOSesc <- plot_statevar_ts(SMSE, "pHOSesc", s, figure = FALSE) %>%
    apply(2, FUN)
  pbrood <- plot_statevar_ts(SMSE, "pbrood", s, figure = FALSE) %>%
    apply(2, FUN)

  x <- cbind(pNOSesc, pHOSesc, pbrood)

  if (any(x > 1, na.rm = TRUE)) warning("Spawners exceed escapement. Error in reporting?")

  if (figure) {
    if (missing(ylim)) ylim <- c(0, 1)
    matplot(Year, x, type = "n", xlab = "Projection Year", ylab = "Proportion",
            ylim = ylim)

    legnames <- c("Spawner/Escapement (natural)", "Spawner/Escapement (hatchery)", "Broodstock/Escapement (total)")

    col <- 1:ncol(x)
    pch <- c(1, 4, 16, 18)
    for (i in 1:ncol(x)) {
      x_i <- x[, i]
      lines(Year[!is.na(x_i)], x_i[!is.na(x_i)], type = "o", col = col[i], pch = pch[i])
    }
    legend("bottomleft", legend = legnames, col = col, pch = pch, lwd = 1, bty = "n")
  }

  invisible(x)
}


#' @name plot_statevar_ts
#' @param type Character. For `plot_fishery`, the fishery state variable to plot.
#' @export
plot_fishery <- function(SMSE, s = 1, type = c("catch", "exploit", "harvest"), FUN = median, figure = TRUE, ylim, ylab, ...) {
  type <- match.arg(type)

  var <- switch(
    type,
    "catch" = "K",
    "exploit" = "Ex",
    "harvest" = "U"
  )

  var_vec <- outer(paste0(c("PT", "T"), "_"), c("NOS", "HOS"), FUN = "paste0") %>%
    as.character()

  x <- sapply(var_vec, function(i) {
    plot_statevar_ts(SMSE, var = paste0(var, i), s = s, figure = FALSE, quant = FALSE, agg.fun = max)
  }, simplify = "array") %>%
    apply(2:3, FUN)

  if (figure) {
    Year <- 1:SMSE@proyears

    if (missing(ylim)) ylim <- c(0, 1.1) * range(x, na.rm = TRUE)

    if (missing(ylab)) {
      ylab <- switch(
        type,
        "catch" = "Kept catch",
        "exploit" = "Exploitation rate",
        "harvest" = "Harvest rate"
      )
    }

    plot(Year, NULL, type = 'n', col = "grey40", ylim = ylim, xlab = "Projection Year",
         ylab = ylab, ...)
    col <- 1:length(var_vec)

    for (i in 1:length(var_vec)) {
      x_i <- x[, i]
      ind <- x_i > 0
      lines(Year[ind], x_i[ind], type = 'o', pch = 1, col = col[i], lty = 1)
    }
    legend("topleft", legend = var_vec, col = 1:4, pch = 1, lty = 1, bty = "n")
  }

  invisible(x)
}


#' @name plot_statevar_ts
#'
#' @param xlim Vector. X-axis limits
#' @param type For `plot_Kobe`, the fishery state variable to plot. Whether to plot the exploitation rate for the terminal (T) or pre-terminal fishery (PT).
#' @export
plot_Kobe <- function(SMSE, s = 1, FUN = median, figure = TRUE, xlim, ylim,
                      xlab = expression(NOS/S[MSY]), ylab = expression(U/U[MSY]), type = c("T", "PT")) {

  type <- match.arg(type)

  if (!length(SMSE@Misc$Ref) || !length(SMSE@Misc$Ref[[s]])) return(invisible(data.frame()))

  NOS <- apply(SMSE@NOS[, s, , ], c(1, 3), sum)

  S_SMSY <- apply(NOS/SMSE@Misc$Ref[[s]]["Spawners_MSY", ], 2, FUN)
  if (type == "T") {
    Ex <- apply(SMSE@ExT_NOS[, s, , ], c(1, 3), max)
    Ex_ExMSY <- apply(Ex/SMSE@Misc$Ref[[s]]["UT_MSY", ], 2, FUN)
  } else {
    Ex <- apply(SMSE@ExPT_NOS[, s, , ], c(1, 3), max)
    Ex_ExMSY <- apply(Ex/SMSE@Misc$Ref[[s]]["UPT_MSY", ], 2, FUN)
  }

  if (figure) {
    Sind <- !is.na(S_SMSY) & is.finite(S_SMSY)
    Eind <- !is.na(Ex_ExMSY) & is.finite(Ex_ExMSY)
    ind <- Sind & Eind
    npoints <- sum(ind)

    if (npoints > 0) {
      Year <- seq(1, SMSE@proyears)
      Year_legend <- pretty(Year[ind], n = 4)
      Year_legend <- Year_legend[Year_legend %in% Year]

      cols <- grDevices::hcl.colors(length(Year), palette = "viridis", rev = FALSE)

      if (missing(xlim)) {
        xlim <- c(0, 1.1) * range(S_SMSY, na.rm = TRUE)
      }
      if (missing(ylim)) {
        ylim <- c(0, 1.1) * range(Ex_ExMSY, na.rm = TRUE)
      }

      plot(S_SMSY[ind], Ex_ExMSY[ind], type = "p", pch = 21, bg = cols[ind],
           xlab = xlab, ylab = ylab,
           xlim = xlim, ylim = ylim)
      abline(h = 1, v = 1, lty = 2, col = "grey40")
      legend("topright", legend = Year_legend, pt.bg = cols[Year %in% Year_legend], pch = 21, bty = "n", title = "Year")
    }
  }

  invisible(data.frame(S_SMSY, Ex_ExMSY))
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
#' @param scenario Atomic, vector of faceting variables (same length as z) used to generate a grid of decision tables
#' @param ncol Integer, number of columns for decision table grid, only used if `scenario` is provided
#' @param dir Character, either "h" or "v" to describe how the grid of tables should be organized (horizontally or vertically)
#' @return ggplot object
#' @seealso [plot_statevar_ts()] [plot_tradeoff()]
#' @import ggplot2
#' @export
plot_decision_table <- function(x, y, z, title, xlab, ylab, scenario, ncol = NULL, dir = "v") {
  dt <- data.frame(x = x, y = y, z = z)
  dt$txt <- format(round(dt$z, 2))

  if (!missing(scenario)) dt$scenario <- scenario

  g <- ggplot(dt, aes(factor(.data$x), factor(.data$y), fill = .data$z)) +
    geom_tile(colour = "grey20", alpha = 0.6) +
    geom_text(aes(label = .data$txt)) +
    coord_cartesian(expand = FALSE) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "deeppink", high = "green4", mid = "white", midpoint = 0.5)

  if (!missing(scenario)) {
    g <- g +
      facet_wrap(vars(.data$scenario), ncol = ncol, dir = dir) +
      theme(strip.background = element_blank())
  }

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
#' @param pm1 Numeric or matrix. A vector of values for the first performance metric on the x-axis. Alternatively, provide a three column matrix corresponding
#' to the lower bound, central tendency, and upper bound.
#' @param pm2 Numeric or matrix. A vector of values for the second performance metric on the y-axis (same length as pm1). Alternatively, provide a three column matrix corresponding
#' to the lower bound, central tendency, and upper bound.
#' @param x1 Atomic, vector of values for the first grouping variable. Various levels are represented by colours. Same length as pm1.
#' @param x2 Numeric, vector of values for the second grouping variable. Various levels are represented by shapes. Same length as pm1.
#' @param xlab Character, optional x-axis label
#' @param ylab Character, optional y-axis label
#' @param x1lab Character, optional label for the first grouping variable
#' @param x2lab Character, optional label for the second grouping variable
#' @param scenario Atomic, vector of faceting variables (same length as `pm1`, `pm2`) used to generate a grid of decision tables
#' @param ncol Integer, number of columns for decision table grid, only used if `scenario` is provided
#' @param dir Character, either "h" or "v" to describe how the grid of tables should be organized (horizontally or vertically) , only used if `scenario` is provided
#' @return ggplot object
#' @seealso [plot_statevar_ts()] [plot_decision_table()]
#' @import ggplot2
#' @export
plot_tradeoff <- function(pm1, pm2, x1, x2, xlab, ylab, x1lab, x2lab, scenario, ncol = NULL, dir = "v") {

  if (missing(x1)) x1 <- 0
  if (missing(x2)) x2 <- 1

  dt <- data.frame(
    pm1 = if (is.matrix(pm1)) pm1[, 2] else pm1,
    pm2 = if (is.matrix(pm1)) pm2[, 2] else pm2,
    pm1_lower = if (is.matrix(pm1)) pm1[, 1] else pm1,
    pm2_lower = if (is.matrix(pm1)) pm2[, 1] else pm2,
    pm1_upper = if (is.matrix(pm1)) pm1[, 3] else pm1,
    pm2_upper = if (is.matrix(pm1)) pm2[, 3] else pm2,
    x1 = x1,
    x2 = x2
  )

  if (!missing(scenario)) dt$scenario <- scenario

  g <- ggplot(dt, aes(.data$pm1, .data$pm2, colour = .data$x1, shape = .data$x2)) +
    geom_linerange(aes(xmin = .data$pm1_lower, xmax = .data$pm1_upper), linewidth = 0.25) +
    geom_linerange(aes(ymin = .data$pm2_lower, ymax = .data$pm2_upper), linewidth = 0.25) +
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

  if (!missing(scenario)) {
    g <- g +
      facet_wrap(vars(.data$scenario), ncol = ncol, dir = dir) +
      theme(strip.background = element_blank())
  }

  if (!missing(xlab)) g <- g + labs(x = xlab)
  if (!missing(ylab)) g <- g + labs(y = ylab)
  if (!missing(x1lab)) g <- g + labs(colour = x1lab)
  if (!missing(x2lab)) g <- g + labs(shape = x2lab)

  return(g)
}


#' Compare state variables from simulation runs
#'
#' @description
#' Compare outputs from multiple simulations to evaluate performance across states of nature and/or management levers (identified by colour):
#' - `compare_statevar_ts()` produces a time series for all simulations, or with medians and 95th percentile intervals
#' - `compare_statevar_hist()` produces a histogram or density plot across all simulations for a particular year
#'
#' @inheritParams plot_statevar_ts
#' @param SMSE_list List of SMSE objects for multiple model runs returned by [salmonMSE()]
#' @param names Character vector `length(SMSE_list)` to label individual model runs
#' @param col_vec Character vector `length(SMSE_list)` for custom colour schemes for comparing across model scenarios in figures
#' @returns An array invisibly. Also generates base graphics if `figure = TRUE`
#' @seealso [plot_statevar_ts()] [compare()]
#' @importFrom grDevices hcl.colors
#' @importFrom graphics polygon matlines
#' @export
compare_statevar_ts <- function(SMSE_list, var = "PNI", s = 1, figure = TRUE, xlab = "Projection Year", quant = FALSE, ylab = var,
                                ylim, agg.fun = sum,
                                names, col_vec, ...) {

  if (missing(names)) names <- paste("Scenario", 1:length(SMSE_list))
  if (missing(col_vec)) {
    col_vec <- hcl.colors(length(SMSE_list), palette = "Dark 3", alpha = 1)
    shade_vec <- hcl.colors(length(SMSE_list), palette = "Dark 3", alpha = 0.25)
  } else {
    if (requireNamespace("scales", quietly = TRUE)) {
      shade_vec <- scales::alpha(col_vec, 0.25)
    } else {
      shade_vec <- col_vec
    }
  }
  if (length(s) > 1) stop("s must be a single integer")

  # if quant = FALSE: Sim x year x scenario
  # if quant = TRUE: 3 x year x scenario
  output <- sapply(SMSE_list, plot_statevar_ts, var = var, s = s, figure = FALSE, quant = quant, agg.fun = agg.fun,
                   simplify = "array")

  if (figure) {
    if (missing(ylim)) ylim <- c(0, 1.1) * range(output, na.rm = TRUE)

    Year <- seq(1, dim(output)[2])

    matplot(Year, Year, type = "n", ylim = ylim, xlab = "Projection Year", ylab = ylab, ...)

    for (i in 1:length(SMSE_list)) {
      xplot <- output[, , i]
      if (any(xplot > 0, na.rm = TRUE)) {
        ind <- colSums(xplot, na.rm = TRUE) > 0

        if (!quant) {
          matlines(Year[ind], t(xplot[, ind]), type = 'l', col = col_vec, lty = 1, ...)
        } else {
          #matlines(Year[ind], t(xplot[, ind]), type = 'o', col = col_vec[i], pch = c(NA, 1, NA), lty = c(2, 1, 2), ...)
          polygon(c(Year[ind], rev(Year[ind])), c(xplot[1, ind], rev(xplot[3, ind])),
                  col = shade_vec[i], border = NA)
          lines(Year[ind], xplot[2, ind], type = 'o', col = col_vec[i], pch = 1, lty = 1, ...)
        }
      }
    }
    legend("topleft", legend = names, col = col_vec, lty = 1, pch = 1, bty = "n")
  }

  return(invisible(output))
}

#' @rdname compare_statevar_ts
#' @inheritParams plot_statevar_hist
#' @param type Character, whether to generate a density figure or histogram
#' @importFrom stats density
#' @export
compare_statevar_hist <- function(SMSE_list, var = "PNI", s = 1, y, figure = TRUE, xlab = var,
                                  names, col_vec, type = c("density", "hist"), ...) {

  type <- match.arg(type)

  if (missing(names)) names <- paste("Scenario", 1:length(SMSE_list))
  if (missing(col_vec)) {
    col_vec <- hcl.colors(length(SMSE_list), palette = "Dark 3", alpha = 1)
    shade_vec <- hcl.colors(length(SMSE_list), palette = "Dark 3", alpha = 0.25)
  } else {
    if (requireNamespace("scales", quietly = TRUE)) {
      shade_vec <- scales::alpha(col_vec, 0.25)
    } else {
      shade_vec <- col_vec
    }
  }
  if (length(s) > 1) stop("s must be a single integer")

  output <- lapply(SMSE_list, plot_statevar_hist, s = s, y = y, var = var, figure = FALSE)

  if (figure) {
    xlim <- range(unlist(output), na.rm = TRUE)

    if (type == "hist") {
      hist_all <- hist(unlist(output), plot = FALSE)
      hist_i <- lapply(output, hist, breaks = hist_all$breaks, plot = FALSE)
      ylim <- c(0, 1.1) * range(sapply(hist_i, getElement, "counts"))

      for (i in 1:length(SMSE_list)) {
        plot(hist_i[[i]], xlab = xlab, main = NULL,
             col = shade_vec[i], border = col_vec[i], add = i > 1, ylim = ylim, ...)
      }

    } else {
      output_dens <- lapply(output, function(i) try(density(i), silent = TRUE))

      is_dens <- sapply(output_dens, inherits, "density")
      xd <- sapply(output_dens[is_dens], getElement, "x")
      yd <- sapply(output_dens[is_dens], getElement, "y")
      ylim <- c(0, 1.1) * range(yd)
      col_vec2 <- col_vec[is_dens]
      shade_vec2 <- shade_vec[is_dens]

      if (any(is_dens)) {
        plot(NULL, NULL, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = "Density")
        for (i in 1:sum(is_dens)) {
          polygon(c(xd[, i], rev(xd[, i])), c(yd[, i], rep(0, nrow(yd))), border = col_vec2[i], col = shade_vec2[i])
        }
      }

    }
    legend("topleft", legend = names, col = col_vec, lty = 1, pch = 1, bty = "n")
  }

  return(invisible(output))
}
