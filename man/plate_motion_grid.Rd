% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequences.R
\name{plate_motion_grid}
\alias{plate_motion_grid}
\title{Plate motion grid}
\usage{
plate_motion_grid(euler, gridsize = 5)
}
\arguments{
\item{euler}{data.frame. containing lat, lon, and angle (optional) of Euler rotation}

\item{gridsize}{grid size in degree}
}
\value{
data.frame with plate motion direction at grid point
}
\description{
Create a grid of the plate motion direction and velocity
}
\examples{
data(pangea)
euler <- subset(pangea, pangea$plate.rot == 103 & pangea$age == 250)
plate_motion_grid(euler)
}