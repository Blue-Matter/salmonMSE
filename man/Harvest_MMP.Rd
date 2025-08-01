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
  K_PT,
  K_T,
  type_PT = "u",
  type_T = "u",
  MSF_PT = FALSE,
  MSF_T = FALSE,
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

\item{DataList}{A nested list of \link[MSEtool:Data-class]{MSEtool::Data} objects by stock then fleet, generated by \code{\link[MSEtool:multiMSE]{MSEtool::multiMSE()}}}

\item{reps}{The number of stochastic replicates to be returned by the function}

\item{u_terminal}{Numeric vector by population (s). Harvest rate of retained catch in the terminal fishery}

\item{u_preterminal}{Single numeric. Harvest rate of retained catch in the pre-terminal fishery}

\item{K_PT}{Single numeric. Total retained catch in the pre-terminal fishery}

\item{K_T}{Numeric vector by population (s). Total retained catch in the terminal fishery}

\item{type_PT}{Single character, containing either "catch" or "u". Indicates whether to manage by harvest rate or total catch numbers for preterminal fisheries}

\item{type_T}{Character vector by population (s), containing either "catch" or "u". indicates whether to manage by harvest rate or total catch numbers for terminal fisheries}

\item{MSF_PT}{Single logical, whether to implement mark-selective fishing for the preterminal fishery}

\item{MSF_T}{Logical vector by population (s), whether to implement mark-selective fishing for the terminal fishery}

\item{m}{Numeric vector by population (s). Mark rate of hatchery origin fish, as a proxy for fishery retention. Only used to calculate the fishing effort.
Retention in the operating model is specified in the \link[MSEtool:MOM-class]{MSEtool::MOM} object}

\item{release_mort}{Matrix \verb{[2, s]}. Release mortality of discarded fish in the pre-terminal (1st row) and terminal (2nd row) fishery. Only used
if \code{MSF_PT = TRUE} or \code{MSF_T[s] = TRUE}. Only used to calculate the fishing effort.
Release mortality in the operating model is specified in the \link[MSEtool:MOM-class]{MSEtool::MOM} object}

\item{p_terminal}{Numeric vector. Population index (p) for the recruitment that experiences the terminal fishing mortality}

\item{p_preterminal}{Numeric vector. Population index (p) for immature fish that experience the pre-terminal fishing mortality}

\item{pkey}{Data frame that assigns openMSE population by life stage and origin (p) to salmonMSE population (s)}

\item{p_natural}{Numeric vector. Population index (p) for natural origin fish}

\item{p_hatchery}{Numeric vector. Population index (p) for hatchery origin fish}

\item{...}{Not used}
}
\value{
A nested list of \link[MSEtool:Rec-class]{MSEtool::Rec} objects, same dimension as \code{DataList}
}
\description{
A function used by openMSE to specify the fishing effort during the projections.
salmonMSE updates the arguments of this function from the salmon operating model.
}
\keyword{internal}
