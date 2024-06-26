% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequences.R
\name{read.gplates}
\alias{read.gplates}
\title{Read GPLATES' rotation file}
\usage{
read.gplates(x, ...)
}
\arguments{
\item{x}{either a file name (.rot format) OR a matrix}

\item{...}{optional arguments to the read.table function}
}
\value{
object of class \code{"finite"}
\describe{
  \item{plate.rot}{ID of moving plate}
  \item{lat}{Latitude of Euler pole of total reconstruction rotation}
  \item{lon}{Longitude}
  \item{angle}{Rotation angle in degree}
  \item{plate.fix}{ID of fixed/anchored plate}
  \item{cmt}{Comments}
}
}
\description{
Imports a sequence of of total reconstruction rotations from a
GPLATES' .rot file
}
\details{
The comment column (last column) must not include **white space**. Use "_" to separate words instead.
}
\examples{
fname <- system.file("Pangea.rot", package="ptrotR")
Pangea <- read.gplates(fname)
print(Pangea)
}
\seealso{
\code{\link{check.finite}}
}
