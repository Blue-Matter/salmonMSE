

#' @import MSEtool
#' @importFrom dplyr %>% pull filter
#' @importFrom methods new slot slotNames slot<-
#' @importFrom stats rlnorm
#' @importFrom rlang .data .env
NULL

setClassUnion("num.array", c("numeric", "array"))
setClassUnion("num.logical", c("numeric", "logical"))
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
    maxage = "numeric",
    n_g = "numeric",
    p_LHG = "numeric",
    p_mature = "num.array",
    SRrel = "character",
    capacity = "numeric",
    kappa = "numeric",
    Smax = "numeric",
    phi = "numeric",
    Mjuv_NOS = "num.array",
    fec = "num.array",
    p_female = "numeric",
    s_enroute = "numeric"
    #strays = 0
  )
)



# ---- Habitat Class -----
#' Class \code{"Habitat"}
#'
#' The component of the operating model that controls survival in the freshwater environment. Includes changes in survival from either
#' environmental/climate effects or habitat mitigation.
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
    use_habitat = "logical",
    prespawn_rel = "character",
    prespawn_prod = "numeric",
    prespawn_capacity = "numeric",
    egg_rel = "character",
    egg_prod = "numeric",
    egg_capacity = "numeric",
    fry_rel = "character",
    fry_prod = "numeric",
    fry_capacity = "numeric",
    fry_sdev = "matrix",
    smolt_rel = "character",
    smolt_prod = "numeric",
    smolt_capacity = "numeric",
    smolt_sdev = "matrix"
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
    type_PT = "character",
    type_T = "character",
    u_preterminal = "numeric",
    u_terminal = "numeric",
    K_PT = "numeric",
    K_T = "numeric",
    MSF_PT = "logical",
    MSF_T = "logical",
    release_mort = "numeric",
    vulPT = "num.matrix",
    vulT = "num.matrix"
  )
)


# ---- Historical Class -----
#' Class \code{"Historical"}
#'
#' Optional component of the operating model that specifies the historical dynamics.
#'
#' Several approaches are possible:
#' - No set up. Default option sets 1000 natural-origin juveniles (age 1), and 1000 hatchery-origin juveniles (age 1) if there is hatchery production (otherwise, zero).
#' - *Recommended option*: specify the initial spawning abundance in the terminal age class.
#' - Detailed setup that reconstructs a historical population by specifying the juvenile abundance (at the beginning of the year),
#' annual fishing mortality rates, and spawner abundance. Typically used if there an estimation/conditioning model is used to inform
#' parameters of the operating model.
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
    HistSpawner_NOS = "num.array",
    HistSpawner_HOS = "num.array",
    HistNjuv_NOS = "array",
    HistNjuv_HOS = "array",
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
    n_r = "numeric",
    n_yearling = "numeric",           # Management lever. No hatchery if both this line and next line are zero
    n_subyearling = "numeric",        # Management lever. No hatchery if both this line and previous line are zero
    s_prespawn = "numeric",           # Survival prior to spawning
    s_egg_smolt = "numeric",          # Survival of eggs in hatchery
    s_egg_subyearling = "numeric",
    brood_import = "numeric",
    Mjuv_HOS = "num.array",
    p_mature_HOS = "num.array",
    stray_external = "matrix",
    gamma = "numeric",
    m = "numeric",
    pmax_esc = "numeric",
    pmax_NOB = "numeric",
    ptarget_NOB = "numeric",
    phatchery = "num.logical",
    premove_HOS = "numeric",
    fec_brood = "num.array",
    fitness_type = "character",
    theta = "numeric",
    rel_loss = "numeric",
    zbar_start = "num.array",
    phenotype_variance = "numeric",
    fitness_variance = "numeric",
    heritability = "numeric",
    fitness_floor = "numeric"
  )
)

#' @importFrom methods .hasSlot
setMethod("initialize", "Hatchery",
          function(.Object, ...) {
            dots <- list(...)
            for (i in names(dots)) {
              if (.hasSlot(.Object, i)) slot(.Object, i) <- dots[[i]]
            }
            if (!length(.Object@phatchery)) .Object@phatchery <- NA_real_
            return(.Object)
          })


setClassUnion("Bio.list", c("Bio", "list"))
setClassUnion("Habitat.list", c("Habitat", "list"))
setClassUnion("Hatchery.list", c("Hatchery", "list"))
setClassUnion("Harvest.list", c("Harvest", "list"))
setClassUnion("Historical.list", c("Historical", "list"))


#' Class \code{"SOM"}
#'
#' An object containing all the parameters for a salmon operating model (SOM).
#'
#' @name SOM-class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("SOM", Bio, Habitat, Hatchery, Harvest, Historical)}.
#'
#' @slot Name Character. Identifying name
#' @slot nsim Integer. Number of simulations
#' @slot nyears Integer. The number of historical years
#' @slot proyears Integer. The number of projected years
#' @slot seed Integer. A random seed to ensure users can reproduce results exactly
#' @slot Bio \linkS4class{Bio} object informing biological parameters, natural production, and habitat effects. Provide a list of Bio objects for multi-population models.
#' @slot Habitat \linkS4class{Habitat} object containing management levers for habitat mitigation. Provide a list of Habitat objects for multi-population models.
#' @slot Hatchery \linkS4class{Hatchery} object containing management levers for hatchery production. Provide a list of Hatchery objects for multi-population models.
#' @slot Harvest \linkS4class{Harvest} object containing management levers for harvest. Provide a list of Harvest objects for multi-population models.
#' @slot Historical \linkS4class{Historical} object to inform historical reconstruction and informing starting abundance for the projection. Provide a list of Historical objects for multi-population models.
#' @slot stray Matrix `[np, np]` where `np = length(Bio)` and row `p` indicates the re-assignment of hatchery fish to each population when they mature (at the recruitment life stage). For example,
#' `SOM@stray <- matrix(c(0.75, 0.25, 0.25, 0.75), 2, 2)` indicates that 75 percent of mature fish return to their natal river and 25 percent stray in both populations. By default, an identity matrix is used (no straying).
#' @keywords classes
#'
#' @export
SOM <- setClass(
  "SOM",
  slots = c(
    Name = "character",
    nsim = "numeric",
    nyears = "numeric",
    proyears = "numeric",
    seed = "numeric",
    Bio = "Bio.list",
    Habitat = "Habitat.list",
    Hatchery = "Hatchery.list",
    Harvest = "Harvest.list",
    Historical = "Historical.list",
    stray = "array"
  )
)

#' @importFrom utils packageVersion
setMethod("initialize", "SOM",
          function(.Object, Bio, Habitat, Hatchery, Harvest, Historical,
                   nsim = 3, nyears = 2, proyears = 20, seed = 1, Name = "Salmon operating model", ...) {

            dots <- list(...)

            .Object@nsim <- nsim
            .Object@nyears <- nyears
            .Object@proyears <- proyears
            .Object@seed <- seed

            if (!missing(Bio)) .Object@Bio <- Bio
            if (!missing(Habitat)) .Object@Habitat <- Habitat
            if (!missing(Hatchery)) .Object@Hatchery <- Hatchery
            if (!missing(Harvest)) .Object@Harvest <- Harvest
            if (!missing(Historical)) .Object@Historical <- Historical

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
#' @slot Egg_NOS Array `[nsim, nstocks, proyears]`. Spawning output, i.e., egg production, of natural origin spawners.
#' @slot Egg_HOS Array `[nsim, nstocks, proyears]`. Spawning output of hatchery origin spawners.
#' @slot Fry_NOS Array `[nsim, nstocks, proyears]`. Fry that are offspring of natural origin spawners.
#' @slot Fry_HOS Array `[nsim, nstocks, proyears]`. Fry that are offspring of hatchery origin spawners.
#' @slot Smolt_NOS Array `[nsim, nstocks, proyears]`. Smolts that are offspring of natural origin spawners.
#' @slot Smolt_HOS Array `[nsim, nstocks, proyears]`. Smolts that are offspring of hatchery origin spawners.
#' @slot Smolt_Rel Array `[nsim, nstocks, proyears]`. Smolts that are offspring of broodtake, i.e., hatchery releases.
#' @slot Njuv_NOS Array `[nsim, nstocks, nage, proyears]`. Abundance of juvenile natural origin fish at the beginning of the year.
#' @slot Njuv_HOS Array `[nsim, nstocks, nage, proyears]`. Abundance of juvenile hatchery origin fish at the beginning of the year.
#' @slot Return_NOS Array `[nsim, nstocks, nage, proyears]`. Mature fish that will be natural origin spawners.
#' @slot Return_HOS Array `[nsim, nstocks, nage, proyears]`. Mature fish that will be hatchery origin spawners.
#' @slot Escapement_NOS Array `[nsim, nstocks, nage, proyears]`. The escapement of mature fish that will be natural origin spawners.
#' @slot Escapement_HOS Array `[nsim, nstocks, nage, proyears]`. The escapement of mature fish that will be hatchery origin spawners.
#' @slot NOB Array `[nsim, nstocks, proyears]`. Natural origin broodtake.
#' @slot HOB Array `[nsim, nstocks, proyears]`. Hatchery origin broodtake (local + strays).
#' @slot HOB_stray Array `[nsim, nstocks, proyears]`. Hatchery origin broodtake (strays only).
#' @slot HOB_import Array `[nsim, nstocks, proyears]`. Imported hatchery origin broodtake used for hatchery production.
#' @slot NOS Array `[nsim, nstocks, nage, proyears]`. Natural origin spawners.
#' @slot HOS Array `[nsim, nstocks, nage, proyears]`. Hatchery origin spawners (local + strays).
#' @slot HOS_stray Array `[nsim, nstocks, nage, proyears]`. Hatchery origin spawners (strays only).
#' @slot HOS_effective Array `[nsim, nstocks, nage, proyears]`. Hatchery origin spawners (local + strays) discounted by `gamma`.
#' @slot KPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery kept catch of natural origin spawners.
#' @slot KT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery kept catch of natural origin spawners.
#' @slot KPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery kept catch of hatchery origin spawners.
#' @slot KT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery kept catch of hatchery origin spawners.
#' @slot DPT_NOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery released catch (live and dead) of natural origin spawners.
#' @slot DT_NOS Array `[nsim, nstocks, proyears]`. Terminal fishery released catch (live and dead) of natural origin spawners.
#' @slot DPT_HOS Array `[nsim, nstocks, proyears]`. Pre-terminal fishery released catch (live and dead) of hatchery origin spawners.
#' @slot DT_HOS Array `[nsim, nstocks, proyears]`. Terminal fishery released catch (live and dead) hatchery origin spawners.
#' @slot UPT_NOS Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery harvest rate (from kept catch) of natural origin spawners.
#' @slot UT_NOS Array `[nsim, nstocks, nage, proyears]`. Terminal fishery harvest rate of natural origin spawners.
#' @slot UPT_HOS Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery harvest rate of hatchery origin spawners.
#' @slot UT_HOS Array `[nsim, nstocks, nage, proyears]`. Terminal fishery harvest rate of hatchery origin spawners.
#' @slot ExPT_NOS Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery exploitation rate (from kept catch and dead releases) of natural origin spawners.
#' @slot ExT_NOS Array `[nsim, nstocks, nage, proyears]`. Terminal fishery exploitation rate of natural origin spawners.
#' @slot ExPT_HOS Array `[nsim, nstocks, nage, proyears]`. Pre-terminal fishery exploitation rate of hatchery origin spawners.
#' @slot ExT_HOS Array `[nsim, nstocks, nage, proyears]`. Terminal fishery exploitation rate of hatchery origin spawners.
#' @slot fitness Array `[nsim, nstocks, 2, proyears]`. Fitness of the population in the natural (1) and hatchery (2) environments.
#' @slot pNOB Array `[nsim, nstocks, proyears]`. Proportion of natural fish in the brood.
#' @slot pHOS_census Array `[nsim, nstocks, proyears]`. Proportion of spawners of hatchery origin, weighted by age class fecundity.
#' @slot pHOS_effective Array `[nsim, nstocks, proyears]`. Proportion of spawners of hatchery origin, discounted by `gamma`, weighted by age class fecundity.
#' @slot PNI Array `[nsim, nstocks, proyears]`. Proportionate natural influence, index of gene flow from hatchery to the natural environment.
#' @slot p_wild Array `[nsim, nstocks, proyears]`. Proportion of wild spawners, natural spawners whose parents were also produced in the natural environment assuming
#' non-assortative mating, defined under Canada's Wild Salmon Policy.
#' @slot Mjuv_loss Array `[nsim, nstocks, nage, proyears]`. Realized juvenile natural mortality, which may differ from inputs due to fitness loss.
#' @slot Misc List. Miscellaneous output:
#'
#' - `Ref` for reference points
#' - `SHist` for the [salmonMSE::SHist-class] object
#' - `SOM` for the [salmonMSE::SOM-class] object.
#' - `LHG` list `nstocks` long containing state variables by life history group
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
    Egg_NOS = "array",
    Egg_HOS = "array",
    Fry_NOS = "array",
    Fry_HOS = "array",
    Smolt_NOS = "array",
    Smolt_HOS = "array",
    Smolt_Rel = "array",
    Njuv_NOS = "array",
    Njuv_HOS = "array",
    Return_NOS = "array",
    Return_HOS = "array",
    Escapement_NOS = "array",
    Escapement_HOS = "array",
    NOB = "array",
    HOB = "array",
    HOB_stray = "array",
    HOB_import = "array",
    NOS = "array",
    HOS = "array",
    HOS_stray = "array",
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
    pNOB = "array",
    pHOS_census = "array",
    pHOS_effective = "array",
    Mjuv_loss = "array",
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



# ---- SHist Class -----
#' Class \code{"SHist"}
#'
#' Stores the outputs from the historical reconstruction of salmon operating models.
#'
#' @name SHist-class
#' @docType class
#' @slot Name Character. Identifying name
#' @slot nyears Integer. The number of historical years
#' @slot nsim Integer. The number of simulations
#' @slot nstocks Integer. The number of stocks
#' @slot Snames Character. Stock names
#' @slot Egg_NOS Array `[nsim, nstocks, nyears]`. Spawning output, i.e., egg production, of natural origin spawners.
#' @slot Egg_HOS Array `[nsim, nstocks, nyears]`. Spawning output of hatchery origin spawners.
#' @slot Smolt Array `[nsim, nstocks, nyears]`. Natural smolt production (sum of offspring of natural and hatchery spawners).
#' @slot Smolt_Rel Array `[nsim, nstocks, proyears]`. Smolts that are offspring of broodtake, i.e., hatchery releases.
#' @slot Njuv_NOS Array `[nsim, nstocks, nage, nyears]`. Abundance of juvenile natural origin fish at the beginning of the year.
#' @slot Njuv_HOS Array `[nsim, nstocks, nage, nyears]`. Abundance of juvenile hatchery origin fish at the beginning of the year.
#' @slot Return_NOS Array `[nsim, nstocks, nage, nyears]`. Mature fish that will be natural origin spawners.
#' @slot Return_HOS Array `[nsim, nstocks, nage, nyears]`. Mature fish that will be hatchery origin spawners.
#' @slot Escapement_NOS Array `[nsim, nstocks, nage, nyears]`. The escapement of mature fish that will be natural origin spawners.
#' @slot Escapement_HOS Array `[nsim, nstocks, nage, nyears]`. The escapement of mature fish that will be hatchery origin spawners.
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
#' @slot Misc List. Miscellaneous output
#' @export
#' @keywords classes
#' @examples
#' showClass("SHist")
setClass(
  "SHist",
  slots = c(
    Name = "character",
    nyears = "numeric",
    nsim = "numeric",
    nstocks = "numeric",
    Snames = "character",
    Egg_NOS = "array",
    Egg_HOS = "array",
    Smolt = "array",
    Smolt_Rel = "array",
    Njuv_NOS = "array",
    Njuv_HOS = "array",
    Return_NOS = "array",
    Return_HOS = "array",
    Escapement_NOS = "array",
    Escapement_HOS = "array",
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
    Misc = "list"
  )
)


setMethod("initialize", "SHist",
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



#' @export
setGeneric("report", function(object, ...) standardGeneric("report"))

#' @name report
#' @title Generate markdown reports
#'
#' @description Generate a markdown report for outcomes from a single operating model projection
#'
#' @param object \linkS4class{SMSE} object
#' @param name Character string for the model name to include in the report, e.g., model run number.
#' @param filename Character string for the name of the markdown and HTML files.
#' @param dir The directory in which the markdown and HTML files will be saved.
#' @param open_file Logical, whether the HTML document is opened after it is rendered.
#' @param render_args List of arguments to pass to [rmarkdown::render()].
#' @param ... Additional arguments (not used)
#' @importFrom utils browseURL
#' @aliases report,SMSE-method
#' @importFrom rmarkdown render
#' @return Returns invisibly the output of [rmarkdown::render()], typically the path of the output file
#' @export
setMethod("report", "SMSE",
          function(object, name = object@Name, filename = "SMSE", dir = tempdir(), open_file = TRUE, render_args = list(), ...) {

            if (missing(name)) name <- substitute(object) %>% as.character()

            dots <- list(...)
            SMSE <- object # Needed for markdown file

            ####### Function arguments for rmarkdown::render
            rmd <- system.file("include", "SMSEreport.Rmd", package = "salmonMSE") %>% readLines()
            rmd_split <- split(rmd, 1:length(rmd))

            name_ind <- grep("NAME", rmd)
            rmd_split[[name_ind]] <- paste("#", name, "{.tabset}")

            stock_ind <- grep("ADD RMD BY STOCK", rmd)
            n_g <- sapply(SMSE@Misc$SOM@Bio, slot, "n_g")
            n_r <- sapply(SMSE@Misc$SOM@Hatchery, slot, "n_r")
            rmd_split[[stock_ind]] <- Map(make_rmd_stock, s = 1:SMSE@nstocks, sname = SMSE@Snames, n_g = n_g, n_r = n_r) %>% unlist()

            if (SMSE@nstocks > 1) {
              rmd_stock_compare <- make_rmd_stock_comparison()
              rmd_split[[stock_ind]] <- c(rmd_split[[stock_ind]], rmd_stock_compare)
            }

            filename_rmd <- paste0(filename, ".Rmd")

            render_args$input <- file.path(dir, filename_rmd)
            if (is.null(render_args$quiet)) render_args$quiet <- TRUE

            # Generate markdown report
            if (!dir.exists(dir)) {
              message("Creating directory: ", dir)
              dir.create(dir)
            }
            write(unlist(rmd_split), file = file.path(dir, filename_rmd))

            # Rendering markdown file
            message("Rendering markdown file: ", file.path(dir, filename_rmd))
            output_filename <- do.call(rmarkdown::render, render_args)
            message("Rendered file: ", output_filename)

            if (open_file) browseURL(output_filename)
            invisible(output_filename)
          })
