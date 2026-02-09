

#' Run salmonMSE
#'
#' @description
#' `salmonMSE()` runs a salmon management strategy evaluation through the following steps:
#' - Converts a salmon operating model (\linkS4class{SOM}) to a multi-stock operating model ([MSEtool::MOM-class]) via [SOM2MOM()]
#' - Creates a harvest management procedure specifying the harvest control rule
#' - Generates the historical reconstruction of the state variables
#' - Runs projection (if `Hist = FALSE`)
#' - Converts the openMSE output, along with additional state variables recorded in [salmonMSE_env], into a salmon MSE object (SMSE) via [MMSE2SMSE()]
#' @param SOM An object of class \linkS4class{SOM}
#' @param Hist Logical, whether to stop the function stop after historical simulations?
#' @param silent Logical, whether to report progress in console
#' @param trace Logical, whether to report additional messages from openMSE
#' @param convert Logical, whether to convert the output into a salmon MSE (SHist or SMSE, depending on `Hist`) object
#' @return
#' If `Hist = TRUE`: if `convert = TRUE`, a \linkS4class{SHist} object or if `convert = FALSE`, a multiHist object (list).
#'
#' If `Hist = FALSE`: if `convert = TRUE`, a \linkS4class{SMSE} object or if `convert = FALSE`, a [MSEtool::MMSE-class] object.
#' @export
salmonMSE <- function(SOM, Hist = FALSE, silent = FALSE, trace = FALSE, convert = TRUE) {

  if (!silent) message("Converting salmon operating model to MOM..")

  SOM <- check_SOM(SOM, silent = silent)
  MOM <- SOM2MOM(SOM, check = FALSE)

  HMMP <- make_Harvest_MMP(SOM, check = FALSE)

  salmonMSE_env$Ford <- data.frame()
  salmonMSE_env$N <- data.frame()
  salmonMSE_env$stateN <- data.frame()
  salmonMSE_env$H <- data.frame()
  salmonMSE_env$stateH <- data.frame()

  if (!silent) message("Generating historical dynamics..")

  H <- SimulateMOM(MOM, parallel = FALSE, silent = !trace)

  # Add numbers at age
  H <- initialize_population(H, SOM)

  if (Hist) {
    if (!silent) message("Returning historical simulations..")

    if (convert) {
      if (!silent) message("Converting to salmon Hist object..")
      SHist <- multiHist2SHist(H, SOM, check = FALSE)
      return(SHist)
    } else {
      return(H)
    }
  }

  if (!silent) message("Running forward projections..")

  # Initialize zbar in data frame
  salmonMSE_env$Ford <- initialize_zbar(SOM)

  M <- ProjectMOM(H, MPs = "HMMP", parallel = FALSE, silent = !trace, checkMPs = FALSE,
                  dropHist = TRUE, extended = FALSE)
  M@multiHist <- H

  if (convert) {
    if (!silent) message("Converting to salmon MSE object..")
    SMSE <- MMSE2SMSE(M, SOM, HMMP, N = salmonMSE_env$N, salmonMSE_env$stateN, Ford = salmonMSE_env$Ford, salmonMSE_env$H, salmonMSE_env$stateH)
    SMSE@Misc$SHist <- multiHist2SHist(H, SOM, check = FALSE)
    return(SMSE)

  } else {
    return(M)
  }

}

initialize_zbar <- function(SOM) {
  pindex <- make_stock_index(SOM)
  ns <- length(SOM@Bio)

  zbar_df <- lapply(1:ns, function(s) {
    do_hatchery <- sum(SOM@Hatchery[[s]]@n_subyearling, SOM@Hatchery[[s]]@n_yearling) > 0
    has_strays <- any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)

    if ((has_strays || do_hatchery) && any(SOM@Hatchery[[s]]@fitness_type == "Ford")) {

      zbar_start <- reshape2::melt(SOM@Hatchery[[s]]@zbar_start)
      nyears_real <- 2
      df <- data.frame(
        x = zbar_start$Var1,
        s = s,
        t = 2 * (nyears_real + 1 - zbar_start$Var2),  # Even time steps relative to first projection year (remember MICE predicts Perr_y for next time step)
        type = ifelse(zbar_start$Var3 == 1, "natural", "hatchery"),
        zbar = zbar_start$value
      )

      fitness_type <- SOM@Hatchery[[s]]@fitness_type
      theta <- SOM@Hatchery[[s]]@theta

      fitness_variance <- SOM@Hatchery[[s]]@fitness_variance
      phenotype_variance <- SOM@Hatchery[[s]]@phenotype_variance

      fitness_floor <- SOM@Hatchery[[s]]@fitness_floor

      df$fitness <- sapply(1:nrow(df), function(i) {
        theta <- switch(df$type[i], "natural" = theta[1], "hatchery" = theta[2])
        fitness_fn <- switch(df$type[i], "natural" = fitness_type[1], "hatchery" = fitness_type[2])

        switch(
          fitness_fn,
          "Ford" = calc_fitness(df$zbar[i], theta, fitness_variance, phenotype_variance, fitness_floor),
          "none" = 1
        )
      })

      return(df)

    } else {
      return(data.frame())
    }
  })

  return(bind_rows(zbar_df))
}


initialize_population <- function(H, SOM) {
  pindex <- make_stock_index(SOM)
  ns <- length(SOM@Bio)
  nyears <- 4
  nareas <- 2

  for (s in 1:ns) {

    Bio <- SOM@Bio[[s]]
    Hatchery <- SOM@Hatchery[[s]]
    a_juv <- seq(1, 2 * (Bio@maxage-1), 2) + 1 # Age classes in second semester of last historical year, leading into first projection year

    do_hatchery <- sum(SOM@Hatchery[[s]]@n_yearling, SOM@Hatchery[[s]]@n_subyearling) > 0
    has_strays <- any(SOM@stray[-s, s] > 0) || sum(SOM@Hatchery[[s]]@stray_external)

    for (g in 1:Bio@n_g) {

      # Add juvenile population (for real age class 2, 3, ...)
      p_NOjuv <- pindex$p[pindex$s == s & pindex$g == g &
                            pindex$stage == "juvenile" &
                            pindex$origin == "natural"]
      H[[p_NOjuv]][[1]]@AtAge$Number[, a_juv, nyears, ] <-
        array(SOM@Historical[[s]]@InitNjuv_NOS[, -1, g]/nareas, c(SOM@nsim, Bio@maxage-1, nareas))

      # Add arbitrary 100 natural-origin spawners in first openMSE year which corresponds to last projection year in real salmon dynamics
      # Will not show up in the results
      p_NOesc <- pindex$p[pindex$s == s & pindex$g == g &
                            pindex$stage == "escapement" &
                            pindex$origin == "natural"]
      H[[p_NOesc]][[1]]@AtAge$Number[, 2 * Bio@maxage, nyears, ] <- array(100/nareas, c(SOM@nsim, 1, nareas))

      # Update Perr_y parameter to parameterize age 1 abundance in first openMSE projection year
      Perr_new <- sapply(1:SOM@nsim, function(x) {
        Egg_openMSE <- sum(100 * Bio@p_female[Bio@maxage] * Bio@fec[x, Bio@maxage, 1])
        SRRpars <- H[[p_NOjuv]][[1]]@SampPars$Stock$SRRpars[[x]]

        Pred_N1 <- H[[p_NOjuv]][[1]]@SampPars$Stock$SRRfun(Egg_openMSE, SRRpars)
        Perr <- SOM@Historical[[s]]@InitNjuv_NOS[x, 1, g]/Pred_N1
        return(Perr)
      })
      nyears <- 4 # Hard coded in
      H[[p_NOjuv]][[1]]@SampPars$Stock$Perr_y[, 2 * Bio@maxage + nyears + 1] <- Perr_new
    }

    if (do_hatchery || has_strays) {
      for (r in 1:Hatchery@n_r) {

        # Add juvenile population
        p_HOjuv <- pindex$p[pindex$s == s & pindex$r == r &
                              pindex$stage == "juvenile" &
                              pindex$origin == "hatchery"]
        H[[p_HOjuv]][[1]]@AtAge$Number[, a_juv, nyears, ] <-
          array(SOM@Historical[[s]]@InitNjuv_HOS[, -1, r]/nareas, c(SOM@nsim, Bio@maxage-1, nareas))

        # Add arbitrary 100 hatchery-origin spawners in first openMSE year which corresponds to last projection year in real salmon dynamics
        # Will not show up in the results
        p_HOesc <- pindex$p[pindex$s == s & pindex$r == r &
                              pindex$stage == "escapement" &
                              pindex$origin == "hatchery"]
        H[[p_HOesc]][[1]]@AtAge$Number[, 2 * Bio@maxage, nyears, ] <- array(100/nareas, c(SOM@nsim, 1, nareas))

        # Update Perr_y parameter to parameterize age 1 abundance in first openMSE projection year
        H[[p_HOjuv]][[1]]@SampPars$Stock$Perr_y[, 2 * Bio@maxage + nyears + 1] <- SOM@Historical[[s]]@InitNjuv_HOS[, 1, r]
      }
    }
  }

  return(H)
}
