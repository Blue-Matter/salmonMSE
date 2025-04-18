% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Harvest_MMP.R
\name{Harvest_MMP}
\alias{Harvest_MMP}
\title{Harvest component of operating model}
\usage{
Harvest_MMP(
  x = 1,
  DataList,
  reps = 1,
  u_terminal,
  u_preterminal,
  MSF = FALSE,
  m,
  release_mort,
  p_terminal = c(2, 5),
  p_preterminal = c(1, 4),
  pkey = data.frame(p = 1:6, s = 1),
  p_natural = 1:3,
  p_hatchery = 4:6,
  ...
)
}
\arguments{
\item{x}{Simulation number}

\item{DataList}{A nested list of \linkS4class{Data} objects by stock then fleet, generated by \code{\link[=multiMSE]{multiMSE()}}}

\item{reps}{The number of stochastic replicates to be returned by the function}

\item{u_terminal}{Numeric vector by population (s). Harvest rate of retained catch in the terminal fishery}

\item{u_preterminal}{Single numeric. Harvest rate of retained catch in the pre-terminal fishery}

\item{MSF}{Logical, whether to implement mark-selective fishing}

\item{m}{Numeric vector by population (s). Mark rate of hatchery origin fish, as a proxy for fishery retention. Only used to calculate the fishing effort.
Retention in the operating model is specified in the \linkS4class{MOM} object}

\item{release_mort}{Matrix \verb{[2, s]}. Release mortality of discarded fish in the pre-terminal (1st row) and terminal (2nd row) fishery. Only used
if \code{MSF = TRUE}. Only used to calculate the fishing effort.
Release mortality in the operating model is specified in the \linkS4class{MOM} object}

\item{p_terminal}{Numeric vector. Population index (p) for the recruitment that experiences the terminal fishing mortality}

\item{p_preterminal}{Numeric vector. Population index (p) for immature fish that experience the pre-terminal fishing mortality}

\item{pkey}{Data frame that assigns openMSE population by life stage and origin (p) to salmonMSE population (s)}

\item{p_natural}{Numeric vector. Population index (p) for natural origin fish}

\item{p_hatchery}{Numeric vector. Population index (p) for hatchery origin fish}

\item{...}{Not used}
}
\value{
A nested list of \linkS4class{Rec} objects, same dimension as \code{DataList}
}
\description{
A function used by openMSE to specify the fishing effort during the projections.
salmonMSE updates the arguments of this function from the salmon operating model.
}
\keyword{internal}
