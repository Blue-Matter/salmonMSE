% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define-SOM.R
\docType{class}
\name{Bio-class}
\alias{Bio-class}
\title{Class \code{"Bio"}}
\description{
The component of the operating model that controls biological dynamics, i.e., natural production.
}
\details{
Various parameters can be stochastic (length \code{nsim}) or input as a single numeric
(value identical across all simulations).
}
\section{Slots}{

\describe{
\item{\code{Name}}{Character. Identifying name}

\item{\code{maxage}}{Integer. The maximum age of the population age structure.}

\item{\code{n_g}}{Integer. Number of life history groups within a cohort.
Life history groups (LHGs) are sub-units of a cohort that have different biological parameters, e.g., survival,
but the egg production and smolt production in the next generation is calculated from the sum across life history groups.
Default is 1.}

\item{\code{p_LHG}}{The proportion of the total egg production assigned to each life history group within a cohort. For example,
if \code{Bio@n_g <- 2}, then \code{Bio@p_LHG <- c(0.9, 0.1)}, then 90 percent of the egg production in the first population is assigned to the first life history group and
ten percent to the second LHG. Default is \code{rep(1/Bio@n_g, Bio@n_g)}}

\item{\code{p_mature}}{Either vector by age (length \code{maxage}) or an array with dimension \verb{[nsim, maxage, nyears+proyears]}. The proportion mature by age.}

\item{\code{SRrel}}{Character, stock-recruit relationship for density-dependent smolt production from fry. Either "BH" (Beverton-Holt) or "Ricker".
Not used if habitat component is used. See \code{Habitat} object.}

\item{\code{capacity}}{Vector length \code{nsim}. The asymptote of the Beverton-Holt stock-recruit function, or the Ricker maximum for
density-dependent natural smolt production from egg production. \strong{Units of smolts.}
Not used if habitat component is used.}

\item{\code{kappa}}{Vector length \code{nsim}. The adult productivity ratio for the stock-recruit function. \strong{Units of recruits per spawner.}
Natural per-capita production of recruits as the population approaches zero (density-independent component).
In stage-based models, equivalent to the product of smolt productivity (smolts per spawner) and marine survival.
Not used if habitat component is used.}

\item{\code{Smax}}{Vector length \code{nsim}. The egg production that maximizes smolt production in the Ricker stock-recruit function. \strong{Units of eggs.}
Equivalent to units of spawners if \code{fec = 1} for all spawners. Not used if habitat component is used.}

\item{\code{phi}}{Optional parameter, vector length \code{nsim}. Unfished egg production per smolt. \strong{Units of egg per smolt}.
The \code{alpha} parameter of the stock-recruit function will be the ratio of \code{kappa} and \code{phi}. In stage-based models,
this is the product of marine survival, fecundity, and proportion female. If not provided, \code{phi} will be calculated from \code{Mjuv_NOS}, \code{p_mature}, \code{s_enroute}, \code{p_female}, \code{fec}, and \code{p_LHG} corresponding
to the first year and weighted by life history groups. Not used if habitat component is used.}

\item{\code{Mjuv_NOS}}{Either vector by age (length \code{maxage}) or an array with dimension \verb{[nsim, maxage, nyears+proyears, n_g]}.
Natural mortality of immature natural origin fish.
To replicate the SAR parameter of a stage-specific model, set \code{Mjuv_NOS[a] = -log(SAR)} where \code{a} is the age class prior to maturation (and zero for all other ages).}

\item{\code{fec}}{Vector by age (length \code{maxage}). Female fecundity of natural origin spawners.}

\item{\code{p_female}}{Numeric. The proportion of females in the spawning population. Default is 0.5.}

\item{\code{s_enroute}}{Numeric. Survival of escapement to the spawning grounds (for spawning and for broodtake). Default is 1.}
}}

\section{Creating Object}{

Objects can be created by calls of the form \code{new("Bio")}
}

\examples{
showClass("Bio")
}
\keyword{classes}
