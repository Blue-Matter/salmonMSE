
.onLoad <- function(libname, pkgname) {
  salmonMSE_env$Ford <- data.frame()
  salmonMSE_env$N <- data.frame()
  salmonMSE_env$stateN <- data.frame()
  salmonMSE_env$H <- data.frame()
  salmonMSE_env$stateH <- data.frame()
}

if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(".", "2.5%", "50%", "97.5%", "Age", "Chain", "Iteration", "Origin",
      "Var1", "Var2", "Var3", "Var4", "Year", "value", "variable",
      "p",
      "Design1", "Design2",
      "glossary", "pm1_lower", "pm1_upper", "pm2_lower", "pm2_upper",
      "Covariate"
    )
  )
}

#' @importFrom graphics grid matplot abline par
plot.default <- function(..., zero_line = FALSE) {
  if (zero_line) {
    graphics::plot.default(..., panel.first = {graphics::grid(); abline(h = 0, col = "grey60")})
  } else {
    graphics::plot.default(..., panel.first = graphics::grid())
  }
}

matplot <- function(..., zero_line = FALSE) {
  if (zero_line) {
    graphics::matplot(..., panel.first = {graphics::grid(); abline(h = 0, col = "grey60")})
  } else {
    graphics::matplot(..., panel.first = graphics::grid())
  }
}

#' @importFrom graphics hist legend
#' @importFrom stats sd
hist.numeric <- function(x, ...) {
  n <- length(x)

  # Remove NA's
  if (any(is.na(x))) {
    na_rate <- mean(is.na(x))
    x <- x[!is.na(x)]
    legend_na <- paste0(round(100 * na_rate), "% NA's")
  } else {
    legend_na <- NULL
  }

  # Calculate skewness NA's already removed, n is the original length of x
  skewness <- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^1.5
  cv <- sd(x)/mean(x)

  # Remove outliers
  if (!is.na(skewness) && cv > 0.1) {
    if (skewness > 3) {
      max_x_plot <- quantile(x, 0.95)
      p_outlier <- signif(100 * mean(x > max_x_plot), 2)
      x <- x[x <= max_x_plot]
      legend_skew1 <- paste0(p_outlier, "% > ", round(max_x_plot, 2))
    } else {
      legend_skew1 <- NULL
    }
    if (skewness < -3) {
      min_x_plot <- quantile(x, 0.05)
      x <- x[x >= min_x_plot]
      p_outlier <- signif(100 * mean(x < min_x_plot), 2)
      legend_skew2 <- paste0(p_outlier, "% < ", round(min_x_plot, 2))
    } else {
      legend_skew2 <- NULL
    }
  } else {
    legend_skew1 <- legend_skew2 <- NULL
  }

  # Plot histogram
  if (all(!diff(signif(x, 3)))) { # If all identical values
    x <- signif(x, 3)
    if (all(!x)) { # if x is all zeros
      breaks <- c(-0.1, 0.1)
      xlim <- c(-1, 1)
    } else {
      breaks <- c(0.99, 1.01) * x[1]
      xlim <- c(x[1] - 0.2 * abs(x[1]), x[1] + 0.2 * abs(x[1]))
    }
    r <- graphics::hist.default(x, breaks = breaks, xlim = xlim, ...)
  } else {
    r <- graphics::hist.default(x, ...)
  }

  # Make legend
  legend_text <- c(legend_na, legend_skew1, legend_skew2)
  if (!is.null(legend_text)) {
    legend("topright", legend = legend_text, bty = "n", text.col = "red")
  }
  invisible(r)
}

#' @name glossary
#' @title salmonMSE glossary
#' @description Glossary of terms and parameters used in salmonMSE
#' @examples
#' data(glossary)
#' glossary[1:2, ]
NULL
