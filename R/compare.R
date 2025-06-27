
compare <- function(..., SMSE_list, Design) {
  dots <- list(...)

  if (length(dots)) {
    SMSE_list <- dots
    class_check <- sapply(SMSE_list, inherits, "SMSE")
    if (!all(class_check)) stop("Not all objects in ... are SMSE objects")
  } else {
    class_check <- sapply(SMSE_list, inherits, "SMSE")
    if (!all(class_check)) stop("Not all objects in SMSE_list are SMSE objects")
  }
}

#' Compare simulation runs
#'
#' Create figures that compare results across two dimensions
#'
#' - [compare_spawners()] generates a time series of the composition of spawners
#' - [compare_fitness()] generates a time series of metrics (fitness, PNI, pHOS, and pWILD) related to hatchery production
#' - [compare_escapement()] generates a time series of the proportion of spawners and broodtake to escapement
#'
#' @param SMSE_list A list of SMSE objects returned by [salmonMSE()]
#' @param Design A data frame with two columns that describes the factorial design of the simulations. Used to label the figure.
#' Rows correspond to each object in `SMSE_list`. There two columns are variables against which to plot the result. See
#' example in \url{https://docs.salmonmse.com/articles/decision-table.html}.
#' @param prop Logical, whether to plot absolute numbers over proportions
#' @param FUN Summarizing function across simulations, typically [stats::median()] or [base::mean()]
#' @returns A ggplot object
#' @export
compare_spawners <- function(SMSE_list, Design, prop = FALSE, FUN = median) {
  Sp <- lapply(1:nrow(Design), function(x) {
    plot_spawners(SMSE_list[[x]], prop = FALSE, FUN = FUN, figure = FALSE) %>%
      reshape2::melt() %>% # Metric x Year
      mutate(Design1 = Design[x, 1], Design2 = Design[x, 2])
  }) %>%
    bind_rows() %>%
    rename(Year = Var2) %>%
    dplyr::filter(value > 0)
  Sp$Var1 <- factor(Sp$Var1, levels = c("HOS", "NOS_notWILD", "WILD"))

  if (prop) {
    Sp <- mutate(Sp, p = value/sum(value), .by = c(Year, Design1, Design2))

    g <- ggplot(Sp, aes(.data$Year, .data$p, fill = .data$Var1)) +
      geom_area() +
      facet_grid(vars(Design1), vars(Design2)) +
      labs(x = "Projection Year", y = "Proportion", fill = NULL) +
      coord_cartesian(ylim = c(0, 1), expand = FALSE) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "PuBuGn")
  } else {
    g <- ggplot(Sp, aes(.data$Year, .data$value, fill = .data$Var1)) +
      geom_area() +
      facet_grid(vars(Design1), vars(Design2)) +
      labs(x = "Projection Year", y = "Spawners", fill = NULL) +
      coord_cartesian(expand = FALSE) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "PuBuGn")
  }
  g
}

#' @rdname compare_spawners
#' @export
compare_fitness <- function(SMSE_list, Design, FUN = median) {

  fitness <- lapply(1:nrow(Design), function(x) {
    plot_fitness(SMSE_list[[x]], figure = FALSE, FUN = FUN) %>%
      reshape2::melt() %>% # Metric x Year
      mutate(Design1 = Design[x, 1], Design2 = Design[x, 2])
  }) %>%
    bind_rows() %>%
    rename(Year = Var1) %>%
    dplyr::filter(value > 0)

  g <- ggplot(fitness, aes(.data$Year, .data$value, colour = .data$Var2)) +
    geom_line(linewidth = 1) +
    facet_grid(vars(Design1), vars(Design2)) +
    labs(x = "Projection Year", y = "Median", colour = "Metric") +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    scale_colour_brewer(palette = "RdYlBu")

  g
}

#' @rdname compare_spawners
#' @export
compare_escapement <- function(SMSE_list, Design, FUN = median) {

  label_fn <- function(x) {
    switch(
      x,
      "pNOSesc" = "Spawner/Escapement (natural)",
      "pHOSesc" = "Spawner/Escapement (hatchery)",
      "pbrood" = "Broodstock/Escapement (total)"
    )
  }

  esc <- lapply(1:nrow(Design), function(x) {
    plot_escapement(SMSE_list[[x]], figure = FALSE, FUN = FUN) %>%
      reshape2::melt() %>% # Metric x Year
      mutate(Design1 = Design[x, 1], Design2 = Design[x, 2])
  }) %>%
    bind_rows() %>%
    rename(Year = Var1) %>%
    dplyr::filter(value > 0) %>%
    mutate(Var2 = sapply(as.character(Var2), label_fn))

  g <- ggplot(esc, aes(.data$Year, .data$value, colour = .data$Var2)) +
    geom_line(linewidth = 1) +
    facet_grid(vars(Design1), vars(Design2)) +
    labs(x = "Projection Year", y = "Median proportion", colour = "Metric") +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    guides(colour = guide_legend(ncol = 2)) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    scale_colour_brewer(palette = "Accent")

  g
}
