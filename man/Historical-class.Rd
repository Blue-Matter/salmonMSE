% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define-SOM.R
\docType{class}
\name{Historical-class}
\alias{Historical-class}
\title{Class \code{"Historical"}}
\description{
The component of the operating model that specifies the historical dynamics.
}
\section{Slots}{

\describe{
\item{\code{Name}}{Character. Identifying name}

\item{\code{HistSpawner_NOS}}{Array by \verb{[nsim, maxage, nyears, n_g]}. Natural origin spawners at age in the natural environment. \emph{Optional}}

\item{\code{HistSpawner_HOS}}{Array by \verb{[nsim, maxage, nyears, n_r]}. Hatchery origin spawners at age in the natural environment. \emph{Optional}}

\item{\code{HistNjuv_NOS}}{Array by \verb{[nsim, maxage, nyears+1, n_g]}. The abundance of immature natural origin fish at the beginning of the annual time step.
Default assumes 1000 smolts (age-1) fish annually.}

\item{\code{HistNjuv_HOS}}{Array by \verb{[nsim, maxage, nyears+1, n_r]}. The abundance of immature hatchery origin fish at the beginning of the annual time step.
Default assumes 1000 smolts (age-1) fish annually.}

\item{\code{HistFPT}}{Vector by historical years (\code{nyears}) or an array by dimension \verb{[nsim, nyears, 2]}. The instantaneous fishing mortality in the preterminal fishery.
The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish. Default is zero.}

\item{\code{HistFT}}{Vector by historical years (\code{nyears}) or an array by dimension \verb{[nsim, nyears, 2]}. The instantaneous fishing mortality in the terminal fishery.
The first array slice corresponds to F for natural origin fish and the second array slice corresponds to hatchery origin fish. Default is zero.}
}}

\section{Creating Object}{

Objects can be created by calls of the form \code{new("Historical")}
}

\examples{
showClass("Historical")
}
\keyword{classes}
