

#' Run salmonMSE
#'
#' @description
#' Runs a salmon management strategy evaluation through the following steps:
#' - Converts a salmon operating model to a multi-stock operating model (MOM) via `SOM2MOM`
#' - Creates a harvest management procedure specifying harvest rate
#' - Generates the historical reconstruction of the state dynamics
#' - Runs projection (if `Hist = TRUE`)
#' @param SOM An object of class \linkS4class{SOM}
#' @param Hist Logical, whether to stop the function stop after historical simulations? Returns a list containing all historical data
#' @param silent Logical, whether to report progress in console bar
#' @param trace Logical, whether to report additional messages from openMSE
#' @export
salmonMSE <- function(SOM, Hist = FALSE, silent = FALSE, trace = FALSE) {

  if (!silent) message("Converting salmon operating model to MOM..")
  MOM <- SOM2MOM(SOM)

  Harvest_MMP <- make_Harvest_MMP(SOM@u)

  if (!silent) message("Generating historical dynamics..")
  H <- SimulateMOM(MOM, parallel = FALSE, silent = !trace)

  if (Hist) {
    if(!silent) message("Returning historical simulations..")
    return(H)
  }

  if (!silent) message("Running forward projections..")
  M <- ProjectMOM(H, MPs = "Harvest_MMP", parallel = FALSE, silent = !trace, checkMPs = FALSE,
                  dropHist = TRUE, extended = FALSE)

  return(M)
}
