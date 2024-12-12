

make_hatchery_table <- function(SMSE, s = 1) {
  hatch_settings <- c("n_yearling", "n_subyearling", "pmax_esc", "pmax_NOB", "premove_HOS",
                      "s_prespawn", "s_egg_smolt", "s_egg_subyearling", "ptarget_NOB",
                      "phatchery")

  hatch_table <- glossary[match(hatch_settings, glossary$Slot), c("Definition", "Slot")]

  hatch_table$Value <- sapply(hatch_table$Slot, function(i) {
    SOM <- SMSE@Misc$SOM
    v <- slot(SOM@Hatchery[[s]], i)
    if (!length(v)) v <- NA
    return(format(v))
  })
  names(hatch_table)[2] <- "Parameter"
  return(hatch_table)
}


make_fitness_table <- function(SMSE, s = 1) {

  fitness_settings <- c("gamma", "fitness_variance", "selection_strength", "heritability", "theta")

  fitness_extra <- data.frame(
    Definition = c("Fitness function for the natural and hatchery environments", "Minimum population fitness value"),
    Slot = c("fitness_type", "fitness_floor")
  )
  fitness_table <- rbind(
    glossary[match(fitness_settings, glossary$Slot), c("Definition", "Slot")],
    fitness_extra
  )
  fitness_table$Value <- sapply(fitness_table$Slot, function(i) {
    SOM <- SMSE@Misc$SOM
    v <- slot(SOM@Hatchery[[s]], i)
    if (!length(v)) v <- NA
    if (length(v) > 1) {
      v <- paste(v, sep = ",")
    }
    return(format(v))
  })
  names(fitness_table)[2] <- "Parameter"
  return(fitness_table)
}

