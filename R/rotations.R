#' @title Antipodal Euler pole
#' @description Euler pole on the other side of the hemisphere
#' @param x data.frame containing the sequence of rotations or the rotation
#' @return Object of class \code{"euler.pole"}
#' @importFrom tectonicr longitude_modulo euler_pole
#' @export
#' @examples
#' euler.pole <- euler_pole(90, 0)
#' antipodal_euler_pole(euler.pole)
antipodal_euler_pole <- function(x) {
  lat.inv <- -x$lat
  lon.inv <- tectonicr::longitude_modulo(x$lon + 180)
  ep.inv <- tectonicr::euler_pole(lat.inv, lon.inv)
  return(ep.inv)
}


#' @title Euler rotation
#' @description Performs a Euler rotation of a vector
#' @param ep Euler pole object
#' @param psi Angle of rotation in degree
#' @param x Point to be rotated. Two-column Vector with  latitude and longitude
#' @return latitude and longitude of rotated vector
#' @export
#' @importFrom tectonicr cartesian_to_geographical geographical_to_cartesian
#' @examples
#' euler.pole <- euler_pole(90, 0)
#' euler_rotation(euler.pole, psi = 45, x = c(45, 45))
euler_rotation <- function(ep, psi, x) {
  x.rot <-
    tectonicr::cartesian_to_geographical(
      tectonicr::euler_rot(ep, psi) %*% tectonicr::geographical_to_cartesian(x)
    )
  return(x.rot)
}


#' @title Normalization of a vector
#' @description normalizes a vector to unit length
#' @param v vector
#' @return normalized vector
#' @export
#' @seealso \code{\link[ppls]{normalize.vector}}
#' @examples
#' normalize_vector(1:5)
normalize_vector <- function(v) {
  v / sqrt(sum(v^2))
}


#' @title Rotation
#' @description Rotates a vector around an axis and an angle
#' @param x vector in cartesian coordinates
#' @param n rotation axis in cartesian coordinates
#' @param alpha angle in degree
#' @details Rotation of a vector is the dot product of the rotation matrix and
#' the vector
#' @return rotated vector in cartesian coordinates
#' @importFrom tectonicr rotation_matrix
#' @export
rotation <- function(x, n, alpha) {
  c(tectonicr::rotation_matrix(n, alpha) %*% x)
}


#' @title Orthodromic distance
#' @description The great-circle distance or orthodromic distance is the
#' shortest distance between two points on the surface of a sphere,
#' measured along the surface of the sphere (as opposed to a straight line
#' through the sphere's interior)
#' @param a lat, lon coordinate of point 1
#' @param b lat, lon coordinate of point 2
#' @return distance in degree
#' @export
#' @examples
#' berlin <- c(52.517, 13.4)
#' tokyo <- c(35.7, 139.767)
#' orthodrome(berlin, tokyo)
orthodrome <- function(a, b){
  delta <- tectonicr::acosd(
    tectonicr::sind(a[1]) * tectonicr::sind(b[1]) +
      tectonicr::cosd(a[1]) * tectonicr::cosd(b[1]) * tectonicr::cosd(b[2]-a[2])
  )
  return(delta)
}

#' @title Great-circle distance distance
#' @description Distance along a great circle.
#' @param a lat, lon coordinates of point 1
#' @param b lat, lon coordinates of point 2
#' @return distance in degree
#' @export
#' @seealso \code{\link{orthodrome}}
#' @examples
#' berlin <- c(52.517, 13.4)
#' tokyo <- c(35.7, 139.767)
#' orthodrome(berlin, tokyo)
greatcircle_distance <- function(a, b){
  if(is.null(dim(a)) == T & is.null(dim(b))){  # a and b are vectors with one latitude/longitude
    delta <- orthodrome(a, b)
  }
  if(!is.null(dim(a)) == T & !is.null(dim(b))) { #a and b are data.frames or vectors with more latitude/longitude
    delta <- c()
    for (i in 1:nrow(a)) {
      delta[i] <- orthodrome(c(a$lat[i], a$lon[i]), c(b$lat[i], b$lon[i]))
    }
  }
  if(!is.null(dim(a)) == T & is.null(dim(b))) { #a is data.frames or vectors with more latitude/longitude and b is a vectors with one latitude/longitude
    delta <- c()
    for (i in 1:nrow(a)) {
      delta[i] <- orthodrome(c(a$lat[i], a$lon[i]), b)
    }
  }
  return(delta)
}
