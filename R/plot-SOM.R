
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

