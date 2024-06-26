% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequences.R
\name{eulerpole_migration}
\alias{eulerpole_migration}
\title{Euler pole migration rate}
\usage{
eulerpole_migration(x)
}
\arguments{
\item{x}{data.frame. Sequence of stage rotations}
}
\value{
data.frame
}
\description{
Calculates the velocity and magnitude of a migrating Euler pole from a sequence of stage rotations
}
\examples{
data(pangea)
stages <- extract_stage_rotations(pangea, plate=103)
eulerpole_migration(stages)
}
