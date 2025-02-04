

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


CMpar_key <- data.frame(
  Parameter = c("log_cr", "log_so", "moadd", "wt", "wto", "wt_sd", "wto_sd",
                "logit_matt", "sd_matt",
                "log_fanomalyPT", "log_FbasePT", "fanomalyPT_sd",
                "log_fanomalyT", "log_FbaseT", "fanomalyT_sd",
                "lnE_sd", "b1", "b", "logit_vulPT", "logit_vulT"),

  Description = c("Log productivity (compensation ratio)",
                  "Log unfished natural spawners",
                  "Additional age 1 natural mortality (M)",
                  "Annual lognormal deviation in egg-smolt mortality",
                  "Annual lognormal deviation in age 1 natural mortality",
                  "Prior standard deviation in egg-smolt mortality deviates",
                  "Prior standard deviation in annual age 1 M deviates",
                  "Logit annual maturity at age",
                  "Prior standard deviation in annual maturity at age",
                  "Annual lognormal deviation in pre-terminal instantaneous fishing mortality (from specified trend)",
                  "Lognormal scaling parameter for pre-terminal fishing mortality (from specified trend)",
                  "Prior standard deviation in pre-terminal fishing mortality deviates",
                  "Annual lognormal deviation in terminal instantaneous fishing mortality (from specified trend)",
                  "Lognormal scaling parameter for terminal fishing mortality (from specified trend)",
                  "Prior standard deviation in terminal fishing mortality deviates",
                  "Lognormal standard deviation in total escapement (observation error)",
                  "Linear coefficients for age 1 mortality covariates",
                  "Linear coefficients for age 2+ mortality covariates",
                  "Logit preterminal fishery vulnerability at age",
                  "Logit terminal fishery vulnerability at age")
)

make_CM_table <- function(fit) {
  model_names <- unique(names(fit$obj$par))
  CMpar_key[match(model_names, CMpar_key$Parameter), ]
}

