

#' @import MSEtool
NULL

#' Class \code{'SOM'}
#'
#' An object containing all the parameters for a salmon operating model.
#'
#' @name SOM-class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('SOM', Bio, Habitat, Hatchery, Harvest)}.

#' @slot Name Name of the operating model
#' @slot nsim The number of simulations (necessarily stochastic, a minimum of 2)
#' @slot nyears The number of historical years
#' @slot proyears The number of projected years
#' @slot seed A random seed to ensure users can reproduce results exactly
#'
#' @template Bio_template
#' @template Habitat_template
#' @template Hatchery_template
#' @template Harvest_template
#' @keywords classes
#'
#' @export
SOM <- setClass(
  "SOM",
  slots = c(
    Name = "character",
    nsim = "numeric",
    ngen = "numeric",
    seed = "numeric"
  ),
  contains = c("Bio", "Habitat", "Hatchery", "Harvest")
)

#' @importFrom utils packageVersion
setMethod("initialize", "SOM", function(.Object, Bio, Habitat, Hatchery, Harvest,
                                        nsim = 2, ngen = 20, seed = 1, ...) {
  dots <- list(...)

  if (is.null(dots$Name)) .Object@Name <- "Salmon operating model"
  .Object@nsim <- nsim
  .Object@ngen <- ngen
  .Object@seed <- seed

  for(i in slotNames(Bio)) slot(.Object, i) <- slot(Bio, i)
  for(i in slotNames(Habitat)) slot(.Object, i) <- slot(Habitat, i)
  for(i in slotNames(Hatchery)) slot(.Object, i) <- slot(Hatchery, i)
  for(i in slotNames(Harvest)) slot(.Object, i) <- slot(Harvest, i)

  attr(.Object, "version") <- paste("salmonMSE", packageVersion("salmonMSE"), "with MSEtool", packageVersion("MSEtool"))
  attr(.Object, "date") <- date()
  attr(.Object, "R.version") <- getRversion()

  return(.Object)
})


