

#' @import MSEtool
#' @importFrom dplyr %>% pull filter
#' @importFrom methods new slot slotNames slot<-
#' @importFrom stats rlnorm
#' @importFrom rlang .data .env
NULL

setClassUnion("num.array", c("numeric", "array"))
setClassUnion("num.matrix", c("numeric", "matrix"))

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
    p_mature = "num.array",         # Age at which adults mature and return
    SRrel = "character",
    capacity_smolt = "numeric",   # Beverton-Holt asymptote. Not unfished capacity!!
    kappa = "numeric",
    Smax = "numeric",
    phi = "numeric",
    Mocean_NOS = "num.array",     # Future feature to allow for time-varying (PDO forcing)
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
    kappa_improve = "numeric"
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
    m = "numeric",
    release_mort = "numeric",
    vulPT = "numeric",
    vulT = "numeric"
  )
)


# ---- Historical Class -----
#' Class \code{"Historical"}
#'
#' The component of the operating model that specifies the historical dynamics.
#'
#' @name Historical-class
#' @docType class
#' @slot Name Character. Identifying name
#' @template Historical_template
#'
#' @section Creating Object:
#' Objects can be created by calls of the form \code{new("Historical")}
#'
#' @export
#' @keywords classes
#' @examples
#' showClass("Historical")
setClass(
  "Historical",
  slots = c(
    Name = "character",
    HistSmolt = "num.matrix",
    HistSpawner = "array",
    HistN = "array",
    HistYearling = "numeric",
    HistFPT = "num.array",
    HistFT = "num.array"
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
    Mocean_HOS = "num.array",         # Future feature to allow for time-varying (PDO forcing)
    gamma = "numeric",
    pmax_NOB = "numeric",
    ptarget_NOB = "numeric",
    phatchery = "numeric",
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
#' \code{new("SOM", Bio, Habitat, Hatchery, Harvest, Historical)}.
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
#' @template Historical_template
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
  contains = c("Bio", "Habitat", "Hatchery", "Harvest", "Historical")
)

#' @importFrom utils packageVersion
setMethod("initialize", "SOM",
          function(.Object, Bio, Habitat, Hatchery, Harvest, Historical,
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
            if (!missing(Historical)) for(i in slotNames(Historical)) slot(.Object, i) <- slot(Historical, i)

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
#' @slot KPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery kept catch of natural origin spawners.
#' @slot KT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery kept catch of natural origin spawners.
#' @slot KPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery kept catch of hatchery origin spawners.
#' @slot KT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery kept catch of hatchery origin spawners.
#' @slot DPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery released catch (live and dead) of natural origin spawners.
#' @slot DT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery released catch (live and dead) of natural origin spawners.
#' @slot DPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery released catch (live and dead) of hatchery origin spawners.
#' @slot DT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery released catch (live and dead) hatchery origin spawners.
#' @slot UPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery harvest rate (from kept catch) of natural origin spawners.
#' @slot UT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery harvest rate of natural origin spawners.
#' @slot UPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery harvest rate of hatchery origin spawners.
#' @slot UT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery harvest rate of hatchery origin spawners.
#' @slot ExPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery exploitation rate (from kept catch and dead releases) of natural origin spawners.
#' @slot ExT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery exploitation rate of natural origin spawners.
#' @slot ExPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery exploitation rate of hatchery origin spawners.
#' @slot ExT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery exploitation rate of hatchery origin spawners.
#' @slot fitness Array `[nsim, nstocks, proyears]`. Fitness of the natural spawning population.
#' @slot PNI Array `[nsim, nstocks, proyears]`. Proportionate natural influence, index of gene flow from hatchery to the natural environment.
#' @slot p_wild Array `[nsim, nstocks, proyears]`. Proportion of wild spawners, natural spawners whose parents were also produced in the natural environment assuming
#' non-assortative mating, defined under Canada's Wild Salmon Policy.
#' @slot SAR_loss Array `[nsim, nstocks, proyears]`. Realized SAR due to fitness loss.
#' @slot Misc List. Miscellaneous output
#'
#' @details
#' In generation \eqn{t}, proportionate natural influence (PNI) is defined as:
#'
#' \deqn{\textrm{PNI}_t = \dfrac{p^\textrm{NOB}_t}{p^\textrm{NOB}_t + p^\textrm{HOSeff}_t}}
#'
#' with \eqn{p^\textrm{HOSeff} = \textrm{HOSeff}/(\textrm{NOS} + \textrm{HOSeff})}.
#'
#' The proportion of wild salmon is defined as:
#'
#' \deqn{p^{\textrm{WILD}}_t = q^\textrm{HOScen}_t
#' \dfrac{(q^\textrm{HOScen}_{t-1})^2}
#' {(q^\textrm{HOScen}_{t-1})^2 + 2\gamma \times p^\textrm{HOScen}_{t-1} q^\textrm{HOScen}_{t-1} +
#' \gamma^2 (p^\textrm{HOScen}_{t-1})^2}}
#'
#' where \eqn{q = 1-p} and \eqn{p^\textrm{HOScen} = \textrm{HOS}/(\textrm{NOS} + \textrm{HOS})}.
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
    KPT_NOS = "array",
    KT_NOS = "array",
    KPT_HOS = "array",
    KT_HOS = "array",
    DPT_NOS = "array",
    DT_NOS = "array",
    DPT_HOS = "array",
    DT_HOS = "array",
    UPT_NOS = "array",
    UT_NOS = "array",
    UPT_HOS = "array",
    UT_HOS = "array",
    ExPT_NOS = "array",
    ExT_NOS = "array",
    ExPT_HOS = "array",
    ExT_HOS = "array",
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

