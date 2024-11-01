

#' Run salmonMSE
#'
#' @description
#' `salmonMSE()` runs a salmon management strategy evaluation through the following steps:
#' - Converts a salmon operating model (\linkS4class{SOM}) to a multi-stock operating model (\linkS4class{MOM}) via [SOM2MOM()]
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
#' If `Hist = FALSE`: if `convert = TRUE`, a \linkS4class{SMSE} object or if `convert = FALSE`, a \linkS4class{MMSE} object.
#' @export
salmonMSE <- function(SOM, Hist = FALSE, silent = FALSE, trace = FALSE, convert = TRUE) {

  if (!silent) message("Converting salmon operating model to MOM..")

  SOM <- check_SOM(SOM, silent = silent)
  MOM <- SOM2MOM(SOM, check = FALSE)

  do_hatchery <- SOM@Hatchery@n_subyearling > 0 || SOM@Hatchery@n_yearling > 0

  HMMP <- make_Harvest_MMP(
    SOM@Harvest@u_terminal,
    SOM@Harvest@u_preterminal,
    SOM@Harvest@MSF,
    SOM@Hatchery@m,
    ifelse(do_hatchery, 0, SOM@Harvest@release_mort)
  )

  salmonMSE_env$Ford <- data.frame()
  salmonMSE_env$N <- data.frame()
  salmonMSE_env$state <- data.frame()

  if (!silent) message("Generating historical dynamics..")
  H <- SimulateMOM(MOM, parallel = FALSE, silent = !trace)

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
  if (do_hatchery && any(SOM@Hatchery@fitness_type == "Ford")) {
    zbar_start <- reshape2::melt(SOM@Hatchery@zbar_start)
    salmonMSE_env$Ford <- data.frame(
      x = zbar_start$Var1,
      p_smolt = 1,
      t = 2 * (SOM@nyears + 1 - zbar_start$Var2),  # Even time steps relative to first projection year (remember MICE predicts Perr_y for next time step)
      type = ifelse(zbar_start$Var3 == 1, "natural", "hatchery"),
      zbar = zbar_start$value
    )
  }

  M <- ProjectMOM(H, MPs = "HMMP", parallel = FALSE, silent = !trace, checkMPs = FALSE,
                  dropHist = TRUE, extended = FALSE)
  M@multiHist <- H

  if (convert) {
    if (!silent) message("Converting to salmon MSE object..")
    SMSE <- MMSE2SMSE(M, SOM, HMMP, N = salmonMSE_env$N, Ford = salmonMSE_env$Ford, salmonMSE_env$state)
    SHist <- multiHist2SHist(H, SOM, check = FALSE)

    Ref <- calc_ref(SOM, check = FALSE)

    SMSE@Misc <- list(
      Ref = Ref,
      SHist = SHist,
      SOM = SOM
    )

    return(SMSE)

  } else {
    return(M)
  }

}
