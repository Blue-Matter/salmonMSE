

#' @import MSEtool
NULL



# ---- Bio Class -----
#' Class \code{"Bio"}
#'
#' The component of the operating model that controls biological dynamics, i.e., natural production.
#'
#' Various parameters can be stochastic (length `nsim`) or input as a single numeric
#' (value identical across all simulations).
#'
#' @name Bio-class
#' @docType class
#' @template Bio_template
#'
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("Bio")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("Bio")
setClass(
  "Bio",
  slots = c(
    Name = "character",
    nsim = "numeric",
    maxage = "numeric",
    p_mature = "numeric",         # Age at which adults mature and return
    capacity_smolt = "numeric",   # Beverton-Holt asymptote. Not unfished capacity!!
    prod_smolt = "numeric",       # Productivity adult/SAR. At unfished, prod_smolt = 1/SAR
    SAR = "numeric",              # Future feature to allow for time-varying (PDO forcing)
    fec = "numeric",              # Spawning fecundity of NOS and HOS
    p_female = "numeric"
    #strays = 0
  )
)



# ---- Habitat Class -----
#' Class \code{"Habitat"}
#'
#' The component of the operating model that controls habitat management.
#'
#' @name Habitat-class
#' @docType class
#' @template Habitat_template
#'
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("Habitat")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("Habitat")
setClass(
  "Habitat",
  slots = c(
    Name = "character",
    capacity_smolt_improve = "numeric",    # Improves Beverton-Holt asymptote by 10% in projection
    prod_smolt_improve = "numeric"           # Keep productivity (SR alpha) constant
  )
)

# ---- Harvest Class -----
#' Class \code{"Harvest"}
#'
#' The component of the operating model that controls harvest.
#'
#' @name Harvest-class
#' @docType class
#' @template Harvest_template
#'
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("Harvest")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("Harvest")
setClass(
  "Harvest",
  slots = c(
    Name = "character",
    u = "numeric",
    m = "numeric"
  )
)



# ---- Hatchery Class -----
#' Class \code{"Hatchery"}
#'
#' The component of the operating model that controls the hatchery management.
#'
#' Various parameters can be stochastic (length `nsim`) or input as a single numeric
#' (value identical across all simulations).
#'
#' @name Hatchery-class
#' @docType class
#' @template Hatchery_template
#'
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("Hatchery")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("Hatchery")
setClass(
  "Hatchery",
  slots = c(
    Name = "character",
    n_yearling = "numeric",           # Management lever. No hatchery if both this line and next line are zero
    n_subyearling = "numeric",        # Management lever. No hatchery if both this line and previous line are zero
    s_prespawn = "numeric",           # Survival prior to spawning
    s_egg_smolt = "numeric",          # Survival of eggs in hatchery
    s_egg_subyearling = "numeric",
    gamma = "numeric",
    pmax_NOB = "numeric",
    ptarget_NOB = "numeric",
    premove_HOS = "numeric",
    theta = "numeric",
    rel_loss = "numeric",
    fec_brood = "numeric",
    fitness_type = "character",
    Zpop_start = "numeric",
    fitness_variance = "numeric",
    selection_strength = "numeric",
    heritability = "numeric",
    fitness_floor = "numeric"
  )
)


#' Class \code{"SOM"}
#'
#' An object containing all the parameters for a salmon operating model (SOM).
#'
#' @name SOM-class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("SOM", Bio, Habitat, Hatchery, Harvest)}.

#' @slot Name Name of the operating modelx
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
    nyears = "numeric",
    proyears = "numeric",
    seed = "numeric"
  ),
  contains = c("Bio", "Habitat", "Hatchery", "Harvest")
)

#' @importFrom utils packageVersion
setMethod("initialize", "SOM",
          function(.Object, Bio, Habitat, Hatchery, Harvest,
                   nyears = 2, proyears = 20, seed = 1, ...) {

  dots <- list(...)

  if (is.null(dots$Name)) .Object@Name <- "Salmon operating model"
  .Object@nyears <- nyears
  .Object@proyears <- proyears
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


