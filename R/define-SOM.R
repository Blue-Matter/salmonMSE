

#' @import MSEtool
#' @importFrom dplyr %>% pull filter
#' @importFrom methods new slot slotNames slot<-
#' @importFrom stats rlnorm
#' @importFrom rlang .data .env
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
#' @slot Name Character. Identifying name
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
#' @slot Name Character. Identifying name
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
#' @slot Name Character. Identifying name
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
    u_preterminal = "numeric",
    u_terminal = "numeric",
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
#' @slot Name Character. Identifying name
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
    fec_brood = "numeric",
    fitness_type = "character",
    theta = "numeric",
    rel_loss = "numeric",
    pbar_start = "numeric",
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
#'
#' @slot Name Character. Identifying name
#' @slot nyears Integer. The number of historical years
#' @slot proyears Integer. The number of projected years
#' @slot seed Integer. A random seed to ensure users can reproduce results exactly
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

            if (!missing(Bio)) for(i in slotNames(Bio)) slot(.Object, i) <- slot(Bio, i)
            if (!missing(Habitat)) for(i in slotNames(Habitat)) slot(.Object, i) <- slot(Habitat, i)
            if (!missing(Hatchery)) for(i in slotNames(Hatchery)) slot(.Object, i) <- slot(Hatchery, i)
            if (!missing(Harvest)) for(i in slotNames(Harvest)) slot(.Object, i) <- slot(Harvest, i)

            attr(.Object, "version") <- paste("salmonMSE", packageVersion("salmonMSE"), "with MSEtool", packageVersion("MSEtool"))
            attr(.Object, "date") <- date()
            attr(.Object, "R.version") <- getRversion()

            return(.Object)
          })


# ---- SMSE Class -----
#' Class \code{"SMSE"}
#'
#' Stores the outputs from the simulation of salmon operating models.
#'
#' @name SMSE-class
#' @docType class
#' @slot Name Character. Identifying name
#' @slot nyears Integer. The number of historical years
#' @slot proyears Integer. The number of projected years
#' @slot nsim Integer. The number of simulations
#' @slot nstocks Integer. The number of stocks
#' @slot Snames Character. Stock names
#' @slot Fry_NOS Array `[nsim, nstocks, proyears]`. Spawning output of natural origin spawners.
#' @slot Fry_HOS Array `[nsim, nstocks, proyears]`. Spawning output of hatchery origin spawners.
#' @slot Smolt_NOS Array `[nsim, nstocks, proyears]`. Smolts that are offspring of natural origin spawners.
#' @slot Smolt_HOS Array `[nsim, nstocks, proyears]`. Smolts that are offspring of hatchery origin spawners.
#' @slot Smolt_Rel Array `[nsim, nstocks, proyears]`. Smolts that are offspring of broodtake, i.e., hatchery releases.
#' @slot Return_NOS Array `[nsim, nstocks, nage, proyears]`. Mature fish that will be natural origin spawners.
#' @slot Return_HOS Array `[nsim, nstocks, nage, proyears]`. Mature fish that will be hatchery origin spawners.
#' @slot Escapement_NOS Array `[nsim, nstocks, nage, proyears]`. The escapement of mature fish that will be natural origin spawners.
#' @slot Escapement_HOS Array `[nsim, nstocks, nage, proyears]`. The escapement of mature fish that will be hatchery origin spawners.
#' @slot NOB Array `[nsim, nstocks, proyears]`. The broodtake of natural origin spawners.
#' @slot HOB Array `[nsim, nstocks, proyears]`. The broodtake of hatchery origin spawners.
#' @slot NOS Array `[nsim, nstocks, proyears]`. Natural origin spawners.
#' @slot HOS Array `[nsim, nstocks, proyears]`. Hatchery origin spawners.
#' @slot HOS_effective Array `[nsim, nstocks, proyears]`. Hatchery origin spawners discounted by `gamma`.
#' @slot CatchPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery catch of natural origin spawners.
#' @slot CatchT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery catch of natural origin spawners.
#' @slot CatchPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery catch of hatchery origin spawners.
#' @slot CatchT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery catch of hatchery origin spawners.
#' @slot UPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery harvest rate of natural origin spawners.
#' @slot UT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery harvest rate of natural origin spawners.
#' @slot UPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery harvest rate of hatchery origin spawners.
#' @slot UT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery harvest rate of hatchery origin spawners.
#' @slot fitness Array `[nsim, nstocks, proyears]`. Fitness.
#' @slot PNI Array `[nsim, nstocks, proyears]`. Proportionate natural influence of hatchery on mean phenotypic values.
#' @slot p_wild Array `[nsim, nstocks, proyears]`. Proportion of wild spawners, defined under Canada's Wild Salmon Policy.
#' @slot SAR_loss Array `[nsim, nstocks, proyears]`. Realized SAR due to fitness loss.
#' @slot Misc List. Miscellaneous output
#'
#' @details
#' In generation \eqn{t}, proportionate natural influence (PNI) is defined as:
#'
#' \deqn{\textrm{PNI}_t = \dfrac{p^{\textrm{NOB}}_t}{p^{\textrm{NOB}}_t + p^{\textrm{HOS}}_t}}
#'
#' The proportion of wild salmon is defined as:
#'
#' \deqn{p^{\textrm{WILD}}_t = q^{\textrm{HOS}}_t
#' \dfrac{(q^{\textrm{HOS}}_{t-1})^2}
#' {(q^{\textrm{HOS}}_{t-1})^2 + 2\gamma \times p^{\textrm{HOS}}_{t-1} q^{\textrm{HOS}}_{t-1} +
#' \gamma^2 (p^{\textrm{HOS}}_{t-1})^2}}
#'
#' where \eqn{q = 1-p}.
#'
#' @references
#' Withler et al. 2018. Genetically Based Targets for Enhanced Contributions to Canadian Pacific Chinook Salmon Populations.
#' DFO Can. Sci. Advis. Sec. Res. Doc. 2018/019. xii + 88 p.
#'
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("SMSE")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("SMSE")
setClass(
  "SMSE",
  slots = c(
    Name = "character",
    nyears = "numeric",
    proyears = "numeric",
    nsim = "numeric",
    nstocks = "numeric",
    Snames = "character",
    Fry_NOS = "array",
    Fry_HOS = "array",
    Smolt_NOS = "array",
    Smolt_HOS = "array",
    Smolt_Rel = "array",
    Return_NOS = "array",
    Return_HOS = "array",
    Escapement_NOS = "array",
    Escapement_HOS = "array",
    NOB = "array",
    HOB = "array",
    NOS = "array",
    HOS = "array",
    HOS_effective = "array",
    CatchPT_NOS = "array",
    CatchT_NOS = "array",
    CatchPT_HOS = "array",
    CatchT_HOS = "array",
    UPT_NOS = "array",
    UT_NOS = "array",
    UPT_HOS = "array",
    UT_HOS = "array",
    fitness = "array",
    SAR_loss = "array",
    PNI = "array",
    p_wild = "array",
    Misc = "list"
  )
)


setMethod("initialize", "SMSE",
          function(.Object, ...) {

            dots <- list(...)
            if (length(dots)) {
              for (i in names(dots)) slot(.Object, i) <- dots[[i]]
            }

            attr(.Object, "version") <- paste("salmonMSE", packageVersion("salmonMSE"), "with MSEtool", packageVersion("MSEtool"))
            attr(.Object, "date") <- date()
            attr(.Object, "R.version") <- getRversion()

            return(.Object)
          })

