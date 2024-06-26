

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
#' @param Hist Logical, whether to stop the function stop after historical simulations? Returns a list containing all historical data
#' @param silent Logical, whether to report progress in console bar
#' @param trace Logical, whether to report additional messages from openMSE
#' @param convert Logical, whether to convert the output into a salmon MSE (SMSE) object
#' @return
#' If `Hist = TRUE` a multiHist object (list). Otherwise, if `convert = TRUE`, a \linkS4class{SMSE} object or if `convert = FALSE`, a \linkS4class{MMSE} object.
#' @export
salmonMSE <- function(SOM, Hist = FALSE, silent = FALSE, trace = FALSE, convert = TRUE) {

  if (!silent) message("Converting salmon operating model to MOM..")
  MOM <- SOM2MOM(SOM)

  Harvest_MMP <- make_Harvest_MMP(SOM@u_terminal, SOM@u_preterminal, SOM@m, SOM@release_mort)

  salmonMSE_env$Ford <- data.frame()
  salmonMSE_env$N <- data.frame()
  salmonMSE_env$state <- data.frame()

  if (!silent) message("Generating historical dynamics..")
  H <- SimulateMOM(MOM, parallel = FALSE, silent = !trace)

  if (Hist) {
    if(!silent) message("Returning historical simulations..")
    return(H)
  }

  if (!silent) message("Running forward projections..")
  M <- ProjectMOM(H, MPs = "Harvest_MMP", parallel = FALSE, silent = !trace, checkMPs = FALSE,
                  dropHist = TRUE, extended = FALSE)
  M@multiHist <- H

  if (convert) {
    if (!silent) message("Converting to salmon MSE object..")
    SMSE <- MMSE2SMSE(M, SOM, Harvest_MMP, N = salmonMSE_env$N, Ford = salmonMSE_env$Ford, salmonMSE_env$state)
  } else {
    return(M)
  }

}
