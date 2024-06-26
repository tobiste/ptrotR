% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequences.R
\name{check.finite}
\alias{check.finite}
\title{Finite rotation object}
\usage{
check.finite(x)
}
\arguments{
\item{x}{data.frame containing columns lat, lon, angle, plate.rot, and plate.fix
\describe{
  \item{plate.rot}{ID of moving plate}
  \item{age}{Age for finite rotation}
  \item{lat}{Latitude of Euler pole of total reconstruction rotation}
  \item{lon}{Longitude}
  \item{angle}{Rotation angle in degree}
  \item{plate.fix}{ID of fixed/anchored plate}
  }}
}
\value{
data.frame
}
\description{
Check if object has columns of a finite rotation sequence
}
\examples{
data(pangea)
check.finite(pangea)
}
