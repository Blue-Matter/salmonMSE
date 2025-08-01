% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ref.R
\name{calc_ref}
\alias{calc_ref}
\alias{calc_MSY}
\alias{calc_Sgen}
\title{Reference points}
\usage{
calc_ref(SOM, rel_F, check = TRUE, maximize = c("MSY", "MER"))

calc_MSY(
  Mjuv,
  fec,
  p_female,
  rel_F,
  vulPT,
  vulT,
  p_mature,
  s_enroute,
  n_g = 1,
  p_LHG = 1,
  SRRpars,
  maximize = c("MSY", "MER"),
  F_search = c(1e-08, 5)
)

calc_Sgen(
  Mjuv,
  fec,
  p_female,
  rel_F,
  vulPT,
  vulT,
  p_mature,
  s_enroute,
  n_g = 1,
  p_LHG = 1,
  SRRpars,
  SMSY,
  F_search = c(1e-08, 100),
  nyears
)
}
\arguments{
\item{SOM}{An object of class \linkS4class{SOM}}

\item{rel_F}{Numeric length 2, indicates the relative effort in the preterminal and terminal fisheries, with a maximum value of 1.
The default is \code{c(0, 1)} which indicates a yield calculation with only the terminal fishery.}

\item{check}{Logical, whether to check the SOM object using \code{\link[=check_SOM]{check_SOM()}}}

\item{maximize}{Character, whether the MSY calculation is the optimum that maximizes catch (\code{"MSY"}) or excess recruitment (\code{"MER"}). The two
methods should be equivalent when \code{rel_F = c(0, 1)}.}

\item{Mjuv}{Numeric \code{maxage} for juvenile natural mortality. Can be a matrix \verb{[maxage, n_g]}.}

\item{fec}{Numeric \code{maxage} for fecundity. Can be a matrix \verb{[maxage, n_g]}.}

\item{p_female}{Numeric for proportion female spawners}

\item{vulPT}{Numeric \code{maxage} for preterminal vulnerability at age}

\item{vulT}{Numeric \code{maxage} for terminal vulnerability at age}

\item{p_mature}{Numeric \code{maxage} for maturity proportions at age. Can be a matrix \verb{[maxage, n_g]}.}

\item{s_enroute}{Numeric for en-route survival of escapement to spawning grounds}

\item{n_g}{Integer, number of life history groups within a cohort}

\item{p_LHG}{Numeric \code{n_g} for proportion of the total egg production assigned to each life history group within a cohort}

\item{SRRpars}{Data frame, one row, that contains the stock recruit parameters that predicts density-dependent survival at the egg-smolt life stage}

\item{F_search}{Numeric, length 2 for the range of F values to search for the instantaneous fishing mortality that produces MSY}

\item{SMSY}{Numeric, spawning abundance at MSY}

\item{nyears}{Integer, number of years to project the population with no fishing to reach \code{SMSY}. Default
is the minimum age of maturity.}
}
\value{
\itemize{
\item \code{calc_MSY} returns a vector of various state variables (catch, exploitation rate, egg production, spawners) at MSY
\item \code{calc_Sgen} returns a numeric
\item \code{calc_ref} returns a list by stock, each containing a matrix of MSY state variables and Sgen by simulation
}
}
\description{
Calculate MSY and Sgen reference points for the operating model. Uses the biological parameters (maturity, natural mortality)
in the last year of the projection.
\itemize{
\item \code{calc_MSY()} calculates the MSY reference points from a set of biological and fishery parameters
\item \code{calc_Sgen()} calculates the Sgen, the spawner abundance that would reach the spawner abundance at MSY after one generation without fishing
\item \code{calc_ref()} is a wrapper function that calculates MSY and Sgen for an operating model
}
}
\seealso{
\code{\link[=calc_Smsy_Ricker]{calc_Smsy_Ricker()}}
}
