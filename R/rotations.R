#' @title Euler pole object
#' @description Creates an object of the orientation of the Euler pole axis
#' @param x latitude or x coordinate of Euler pole axis
#' @param y longitude or y
#' @param z z coordinate
#' @param geo logical, TRUE (default) if Euler pole axis is given in
#' geographical coordinates (latitude, longitude). FALSE if given in cartesian
#' coordinates (x, y, z)
#' @return An object of class \code{"euler.pole"}, i.e., a \code{data.frame}
#' containing the Euler pole axis in geographical and cartesian coordinates.
#' @export
#' @importFrom PlateTectonicStressR geographical_to_cartesian cartesian_to_geographical
#' @examples
#' euler_pole(90, 0)
#' euler_pole(0, 0, 1, geo = FALSE)
euler_pole <- function(x, y, z = NA, geo = TRUE) {
  if (geo) {
    cart <- PlateTectonicStressR::geographical_to_cartesian(c(x, y))
    lat <- x
    lon <- y
    x <- cart[1]
    y <- cart[2]
    z <- cart[3]
  } else {
    pol <- PlateTectonicStressR::cartesian_to_geographical(c(x, y, z))
    lat <- pol[1]
    lon <- pol[2]
    x <- x
    y <- y
    z <- z
  }

  ep <- data.frame(lat, lon, x, y, z)
  class(ep) <- append(class(ep), "euler.pole")
  return(ep)
}

#' @title Antipodal Euler pole
#' @description Creates a Euler pole antidpodal location
#' @param ep Euler pole object
#' @return Object of class \code{"euler.pole"}
#' @export
#' @importFrom PlateTectonicStressR longitude_modulo
#' @examples
#' euler.pole <- euler_pole(90, 0)
#' antipodal_euler_pole(euler.pole)
antipodal_euler_pole <- function(ep) {
  lat.inv <- -ep$lat
  lon.inv <- PlateTectonicStressR::longitude_modulo(ep$lon + 180)
  ep.inv <- euler_pole(lat.inv, lon.inv)
  return(ep.inv)
}


#' @title Euler rotation matrix
#' @description Creates a matrix from the given set of values.
#' @param ep Euler pole object
#' @param psi Angle of rotation in degree
#' @return \code{matrix}
#' @references <div class="csl-entry">Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. <i>Computers and Geosciences</i>, <i>25</i>(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5</div>
#' @export
#' @importFrom pracma cosd sind
#' @examples
#' euler.pole <- euler_pole(90, 0)
#' euler_rot(euler.pole, psi = 45)
euler_rot <- function(ep, psi) {
  mat <- matrix(nrow = 3, ncol = 3)
  mat[1, 1] <- pracma::sind(ep$lat) * pracma::cosd(ep$lon)
  mat[1, 2] <- pracma::sind(ep$lat) * pracma::sind(ep$lon)
  mat[1, 3] <- -pracma::cosd(ep$lat)
  mat[2, 1] <- -pracma::sind(ep$lon)
  mat[2, 2] <- pracma::cosd(ep$lon)
  mat[2, 3] <- 0
  mat[3, 1] <- pracma::cosd(ep$lat) * pracma::cosd(ep$lon)
  mat[3, 2] <- pracma::cosd(ep$lat) * pracma::sind(ep$lon)
  mat[3, 3] <- pracma::sind(ep$lat)

  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- pracma::cosd(psi)
  R[1, 2] <- -pracma::sind(psi)
  R[1, 3] <- 0
  R[2, 1] <- pracma::sind(psi)
  R[2, 2] <- pracma::cosd(psi)
  R[2, 3] <- 0
  R[3, 1] <- 0
  R[3, 2] <- 0
  R[3, 3] <- 1

  A <- solve(mat) %*% R %*% mat
  return(A)
}


#' @title Euler rotation
#' @description Performs a Euler rotation of a vector
#' @param ep Euler pole object
#' @param psi Angle of rotation in degree
#' @param x Point to be rotated. Two-column Vector with  latitude and longitude
#' @return latitude and longitude of rotated vector
#' @export
#' @importFrom PlateTectonicStressR cartesian_to_geographical geographical_to_cartesian
#' @examples
#' euler.pole <- euler_pole(90, 0)
#' euler_rotation(euler.pole, psi = 45, x = c(45, 45))
euler_rotation <- function(ep, psi, x) {
  x.rot <-
    PlateTectonicStressR::cartesian_to_geographical(
      euler_rot(ep, psi) %*% PlateTectonicStressR::geographical_to_cartesian(x)
    )
  return(x.rot)
}


#' @title Normalization of a vector#'
#' @description normalizes a vector to unit length#'
#' @param v vector
#' @return normalized vector
#' @export
#' @examples
#' normalize_vector(1:5)
normalize_vector <- function(v) {
  v / sqrt(sum(v^2))
}

#' @title Rotation matrix
#' @description Creates rotation matrix from a rotation axis and angle
#' @param n vector in cartesian coordinates
#' @param alpha in degree
#' @return 3x3 matrix
#' @importFrom pracma cosd sind
#' @export
#' @examples
#' rotation_matrix(c(0, 0, 1), 45)
rotation_matrix <- function(n, alpha) {
  n <- normalize_vector(n) # unit vector
  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- n[1]^2 * (1 - pracma::cosd(alpha)) + pracma::cosd(alpha)
  R[1, 2] <- n[1] * n[2] * (1 - pracma::cosd(alpha)) - n[3] * pracma::sind(alpha)
  R[1, 3] <- n[1] * n[3] * (1 - pracma::cosd(alpha)) + n[2] * pracma::sind(alpha)
  R[2, 1] <- n[2] * n[1] * (1 - pracma::cosd(alpha)) + n[3] * pracma::sind(alpha)
  R[2, 2] <- n[2]^2 * (1 - pracma::cosd(alpha)) + pracma::cosd(alpha)
  R[2, 3] <- n[2] * n[3] * (1 - pracma::cosd(alpha)) - n[1] * pracma::sind(alpha)
  R[3, 1] <- n[3] * n[1] * (1 - pracma::cosd(alpha)) - n[2] * pracma::sind(alpha)
  R[3, 2] <- n[3] * n[2] * (1 - pracma::cosd(alpha)) + n[1] * pracma::sind(alpha)
  R[3, 3] <- n[3]^2 * (1 - pracma::cosd(alpha)) + pracma::cosd(alpha)
  return(R)
}

#' @title Rotation
#' @description Rotates a vector around an axis and an angle
#' @param x vector in cartesian coordinates
#' @param n rotation axis in cartesian coordinates
#' @param alpha angle in degree
#' @details Rotation of a vector is the dot product of the rotation matrix and
#' the vector
#' @return rotated vector in cartesian coordinates
#' @export
rotation <- function(x, n, alpha) {
  c(rotation_matrix(n, alpha) %*% x)
}


#' @title Rotation angle from rotation matrix
#' @description Extracts the rotation anglr from rotation matrix
#' @param A 3x3 matrix
#' @return angle in degree
#' @importFrom pracma rad2deg
#' @export
rotation_angle <- function(A) {
  psi <- acos((sum(diag(A)) - 1) / 2)
  return(pracma::rad2deg(psi))
}

#' @title Rotation axis from rotation matrix
#' @description Extracts the rotation axis from rotation matrix
#' @param A 3x3 matrix
#' @return vector
#' @importFrom pracma sind
#' @export
rotation_axis <- function(A) {
  psi <- rotation_angle(A)
  e1 <- (A[3, 2] - A[2, 3]) / 2 * pracma::sind(psi)
  e2 <- (A[1, 3] - A[3, 1]) / 2 * pracma::sind(psi)
  e3 <- (A[2, 1] - A[1, 2]) / 2 * pracma::sind(psi)
  return(c(e1, e2, e3))
}


#' @title Euler axis and angle from Euler matrix
#' @description Extracts the Euler pole and angle from a Euler matrix
#' @param A 3x3 matrix
#' @return list
#' \describe{
#' \item{pole}{Euler pole object}
#' \item{psi}{Euler angle in degree}
#' }
#' @export
euler_from_rot <- function(A) {
  psi <- rotation_angle(A)
  ra <- rotation_axis(A)
  ep <- euler_pole(ra[1], ra[2], ra[3], geo = FALSE)
  return(list(pole = ep, psi = psi))
}
