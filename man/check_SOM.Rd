% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SOM2MOM.R
\name{check_SOM}
\alias{check_SOM}
\title{Check inputs to SOM object}
\usage{
check_SOM(SOM, silent = FALSE)
}
\arguments{
\item{SOM}{\linkS4class{SOM} object}

\item{silent}{Logical, whether to report progress in console}
}
\value{
Updated \linkS4class{SOM} object with full dimensions in various slots
}
\description{
Ensures that the slots in the \linkS4class{SOM} object have the correct dimensions. Function will
update some slots to their full dimensions.
}
