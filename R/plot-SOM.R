
plot_SOM <- function(SOM, var = "kappa", figure = TRUE, xlab = var, ...) {
  x <- slot(SOM, var)
  xplot <- x

  if (figure) hist(xplot, xlab = xlab, main = NULL, ...)
  invisible(xplot)
}

#' @importFrom grDevices rainbow
#' @importFrom graphics axis
plot_stray <- function(stray, xlab = "Destination", ylab = "Origin") {

  ns <- nrow(stray)

  vcol <- rainbow(100, end = 0.45)

  graphics::plot.default(
    NULL, xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "i",
    xaxt = "n", yaxt = "n", xlim = c(1, ns+1), ylim = c(1, ns+1)
  )

  for(x in 1:ns) {
    for(y in 1:ns) {
      m_yx <- round(stray[y, x], 2)
      rect(xleft = x, ybottom = y, xright = x+1, ytop = y+1, col = vcol[100 * m_yx])
      text(x + 0.5, y + 0.5, m_yx)
    }
  }

  axis(2, at = 1:ns + 0.5, labels = 1:ns)
  axis(1, at = 1:ns + 0.5, labels = 1:ns)

  invisible(stray)
}

plot_SRR <- function(SOM, s = 1, quant = FALSE, surv = FALSE, figure = TRUE, check = TRUE, ylim, ...) {
  if (check) SOM <- check_SOM(SOM)

  smolt0 <- sapply(1:SOM@nsim, function(x) {
    salmonMSE:::calc_smolt(
      N1 = SOM@Bio[[s]]@phi[x],
      kappa = SOM@Bio[[s]]@kappa[x],
      capacity = SOM@Bio[[s]]@capacity_smolt[x],
      Smax = SOM@Bio[[s]]@Smax[x],
      phi = SOM@Bio[[s]]@phi[x],
      kappa_improve = SOM@Habitat[[s]]@kappa_improve,
      capacity_improve = SOM@Habitat[[s]]@capacity_smolt_improve,
      SRrel = SOM@Bio[[s]]@SRrel,
      per_recruit = TRUE
    )
  })
  E0 <- smolt0 * SOM@Bio[[s]]@phi

  Egg <- seq(0, 1.5 * max(E0), length.out = 100)
  Smolt <- sapply(1:SOM@nsim, function(x) {
    sapply(Egg, function(N) {
      salmonMSE:::calc_smolt(
        N1 = N,
        kappa = SOM@Bio[[s]]@kappa[x],
        capacity = SOM@Bio[[s]]@capacity_smolt[x],
        Smax = SOM@Bio[[s]]@Smax[x],
        phi = SOM@Bio[[s]]@phi[x],
        kappa_improve = SOM@Habitat[[s]]@kappa_improve,
        capacity_improve = SOM@Habitat[[s]]@capacity_smolt_improve,
        SRrel = SOM@Bio[[s]]@SRrel
      )
    })
  })

  if (surv) {
    x <- t(Smolt/Egg)
    x[, 1] <- SOM@Bio[[s]]@kappa/SOM@Bio[[s]]@phi
  } else {
    x <- t(Smolt)
  }

  if (!quant) {
    xplot <- x
  } else {
    xplot <- apply(x, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
  }

  if (figure && any(xplot > 0, na.rm = TRUE)) {
    xlab <- "Egg production"
    ylab <- ifelse(surv, "Egg-smolt survival", "Smolt production")

    if (missing(ylim)) ylim <- c(0, 1.1) * range(xplot, na.rm = TRUE)

    if (!quant) {
      matplot(Egg, t(xplot), type = 'l', col = "grey40", ylim = ylim, lty = 1,
              xlab = xlab, ylab = ylab, ...)
    } else {
      matplot(Egg, t(xplot), type = 'o', pch = c(NA, 1, NA), col = 1, lty = c(2, 1, 2), ylim = ylim,
              xlab = xlab, ylab = ylab, ...)
      abline(v = median(E0), lty = 3)
    }
  }

  invisible(xplot)
}


