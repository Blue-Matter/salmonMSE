
#' Plot life history groups
#'
#' Plot the annual proportions of life history groups at various life stages
#'
#' @param SMSE Class \linkS4class{SMSE} object returned by [salmonMSE()]
#' @param var Character. Slot for the state variables in `SMSE@Misc$LHG[[1]]`.
#' @param type Character to indicate whether to plot proportion or absolute numbers
#' @param s Integer. Population index for multi-population model (e.g., `s = 1` is the first population in the model)
#' @param FUN Summarizing function across simulations, typically [median()] or [mean()]
#' @param xlab Character. Name of time variable for the figure
#' @param figure Logical, whether to generate a figure (set to FALSE if only using the function to return the data matrix)
#' @param ylab Character. Name of the state variable for the figure
#' @param LHG_names Character. Vector of names for the life history groups
#' @importFrom grDevices hcl.colors
#' @seealso [plot_statevar_ts()]
#' @export
#' @importFrom graphics box
plot_LHG <- function(SMSE, var = "NOS", type = c("prop", "abs"), s = 1, FUN = median, figure = TRUE, xlab = "Projection Year",
                     ylab, LHG_names) {

  type <- match.arg(type)
  LHG <- SMSE@Misc$LHG[[s]]

  if (length(LHG)) {
    var <- match.arg(var, choices = names(LHG))
    x <- LHG[[var]]
  } else {
    x <- matrix(1, 1, SMSE@proyears)
  }

  if (length(dim(x)) == 4) { # Sum across age
    x <- apply(x, c(1, 2, 4), sum)
  }

  if (type == "abs") {
    xplot <- apply(x, 2:3, FUN)
  } else {
    xplot <- apply(x, c(1, 3), function(i) i/sum(i)) %>% apply(c(1, 3), FUN)
  }

  if (figure && any(xplot > 0, na.rm = TRUE)) {
    ind <- colSums(xplot, na.rm = TRUE) > 0
    Year <- 1:SMSE@proyears

    if (missing(LHG_names)) LHG_names <- paste("LHG", 1:nrow(xplot))

    col <- grDevices::hcl.colors(nrow(xplot), palette = "Dark 2")

    if (type == "abs") {
      ylim <- c(0, 1.1) * range(colSums(xplot, na.rm = TRUE))
      if (missing(ylab)) ylab <- var
    } else {
      ylim <- c(0, 1)
      if (missing(ylab)) ylab <- paste("Proportion", var)
    }
    plot(Year, NULL, xlim = range(Year) + c(-1, 0),
         xlab = "Projection Year", ylab = ylab,
         typ = "n", ylim = ylim, xaxs = "i", yaxs = "i")
    barplot(xplot, legend.text = LHG_names, space = 0, xlim = range(Year),
            col = col, border = "grey40", add = TRUE, xpd = FALSE)
    box()
  }

  invisible(xplot)
}

