% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define-SOM.R
\docType{class}
\name{Harvest-class}
\alias{Harvest-class}
\title{Class \code{"Harvest"}}
\description{
The component of the operating model that controls harvest.
}
\section{Slots}{

\describe{
\item{\code{Name}}{Character. Identifying name}

\item{\code{u}}{Numeric. The exploitation rate of the return.}

\item{\code{m}}{Numeric. The mark rate of hatchery origin fish, which affects fishery selectivity.}
}}

\section{Creating Object}{

Objects can be created by calls of the form \code{new("Harvest")}
}

\examples{
showClass("Harvest")
}
\keyword{classes}
