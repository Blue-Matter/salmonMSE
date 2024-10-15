
plot_SOM <- function(SOM, var = "kappa", figure = TRUE, xlab = var, ...) {
  x <- slot(SOM, var)
  xplot <- x

  if (figure) hist(xplot, xlab = xlab, main = NULL, ...)
  invisible(xplot)
}
