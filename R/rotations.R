#' @title Antipodal Euler pole
#' @description Euler pole on the other side of the hemisphere
#' @param x data.frame containing the sequence of rotations or the rotation
#' @return Object of class \code{"euler.pole"}
#' @importFrom tectonicr longitude_modulo euler_pole
#' @export
#' @examples
#' euler.pole <- c(90, 0)
#' antipodal_euler_pole(euler.pole)
antipodal_euler_pole <- function(x) {
  lat.inv <- -x[1]
  lon.inv <- tectonicr::longitude_modulo(x[2] + 180)
  c(lat.inv, lon.inv)
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
#' euler.pole <- tectonicr::euler_pole(90, 0, angle = 45)
#' euler_rotation(euler.pole, x = c(45, 45))
euler_rotation <- function(ep, x) {
  x.rot <-
    tectonicr::cartesian_to_geographical(
      euler::euler_matrix(ep %>% eulerpole_2_eulervec()) %*% tectonicr::geographical_to_cartesian(x)
    )
  return(x.rot)
}

#' @title Rotation
#' @description Rotates a vector around an axis and an angle
#' @param x vector in cartesian coordinates
#' @param n rotation axis in cartesian coordinates
#' @param alpha angle in degree
#' @details Rotation of a vector is the dot product of the rotation matrix and
#' the vector
#' @return rotated vector in cartesian coordinates
#' @importFrom euler rotation_matrix
#' @export
rotation <- function(x, n) {
  c(euler::rotation_matrix(n) %*% x)
}


#
# hav <- function(x){
#   sin(x/2)^2
# }
#
# ahav <- function(x){
#   2 * asin(sqrt(x))
# }

#' Angle along great circle on spherical surface
#'
#' Smallest angle between two points on the surface of a sphere, measured along
#' the surface of the sphere
#'
#' @param p1,p2 numeric vector; lat, lon coordinates of point p1 and p2
#' @details \describe{
#' \item{"orthodrome_haversine"}{uses Haversine formula (the default)}
#' \item{"orthodrome_haversine2"}{ uses Haversine formula optimized for 64-bit floating-point numbers}
#' \item{"orthodrome_vincenty"}{uses Vincenty formula for an ellipsoid with equal major and minor axes}
#' }
#' @return angle in radians
#' @name spherical_angle
#' @examples
#' berlin <- c(52.52, 13.41)
#' calgary <- c(51.04, -114.072)
#' orthodrome_haversine2(berlin, calgary)
#' orthodrome_haversine2_2(berlin, calgary)
#' orthodrome_vincenty2(berlin, calgary)
NULL

#' @rdname spherical_angle
#' @export
orthodrome_haversine2 <- function(p1, p2) {
  a <- p1 * pi / 180
  b <- p2 * pi / 180
  acos(sin(a[1]) * sin(b[1]) + cos(a[1]) * cos(b[1]) * cos(b[2] - a[2]))
}

#' @rdname spherical_angle
#' @export
orthodrome_haversine2_2 <- function(p1, p2) {
  a <- p1 * pi / 180
  b <- p2 * pi / 180
  havdlat <- tectonicr:::hav(b[1] - a[1])

  tectonicr:::ahav(
    havdlat + (1 - havdlat - hav(b[1] + a[1])) * hav(b[2] - a[2])
  )
}

#' @rdname spherical_angle
#' @export
orthodrome_vincenty2 <- function(p1, p2) {
  a <- p1 * pi / 180
  b <- p2 * pi / 180
  dlon <- b[2] - a[2]

  y <- sqrt(
    (cos(b[1]) * sin(dlon))^2 + (cos(a[1]) * sin(b[1]) - sin(a[1]) * cos(b[1]) * cos(dlon))^2
  )
  x <- sin(a[1]) * sin(b[1]) + cos(a[1]) * cos(b[1]) * cos(dlon)
  atan2(y, x)
}


#' @title Distances on Earth' surface
#' @description Distances on a sphere about the size of the Earth.
#' @param p1,p2 numeric vector; lat, lon coordinates of point p1 and p2
#' @param sm Angle between pole and small circle.
#' @param r radius of the sphere (default = 6371.0087714 km, i.e. the radius of
#' the Earth)
#' @param method Formula for calculating great circle distance:
#' \describe{
#' \item{"haversine"}{Haversine formula (the default)}
#' \item{"haversine2"}{Haversine formula optimized for 64-bit floating-point numbers}
#' \item{"vincenty"}{Vincenty formula for an ellipsoid with equal major and minor axes}
#' }
#' @return distance in units of r (default = kilometers)
#' @note The Earth is nearly spherical, so spherical formulas give the distance
#' between points on the surface of the Earth correct to within about 0.5%.
#' @details
#' The great-circle distance or orthodromic distance is the
#' shortest distance between two points on the surface of a sphere,
#' measured along the surface of the sphere.
#'
#' The loxodromic distance between two points gives the distance measured along
#' a line of constant bearing.
#'
#' The small circle distance between two points is measured along a small circle
#' given by the angle between the pole and the small circle.
#'
#' @references
#' Imboden, C., Imboden D. (1972) Orthodromic and loxodromic formula for the
#' calculation of distance and direction between ringing and finding place.
#' Vogelwarte 26: 336-346.
#' @name spherical_distance
#' @examples
#' berlin <- c(52.52, 13.41)
#' calgary <- c(51.04, -114.072)
#' dist_greatcircle2(berlin, calgary)
#' dist_loxodrome(berlin, calgary)
#' dist_smallcircle(berlin, calgary, sm = 45)
NULL

#' @rdname spherical_distance
#' @export
dist_greatcircle2 <- function(p1, p2, r = earth_radius(), method = c("haversine", "haversine2", "vincenty")) {
  method <- match.arg(method)
  stopifnot(is.numeric(r))

  if (method == "haversine") {
    d <- orthodrome_haversine2(p1, p2)
  }
  if (method == "haversine2") {
    d <- orthodrome_haversine2_2(p1, p2)
  }
  if (method == "vincenty") {
    d <- orthodrome_vincenty2(p1, p2)
  }
  d * r
}

#' @rdname spherical_distance
#' @export
dist_loxodrome <- function(p1, p2, r = earth_radius()) {
  stopifnot(is.numeric(r))

  p1 <- p1 * pi / 180
  p2 <- p2 * pi / 180

  dLat <- p2[1] - p1[1]
  dLon <- abs(p2[2] - p1[2])
  dPhi <- log(tan(p2[1] / 2 + pi / 4) / tan(p1[1] / 2 + pi / 4))
  i <- abs(dLat) > 1e-12
  q <- vector(length = length(i))
  q[i] <- dLat[i] / dPhi[i]
  q[!i] <- cos(p1[1][!i]) # E-W course becomes ill-conditioned with 0/0

  # if dLon over 180° take shorter rhumb line across the anti-meridian
  dLon[dLon > pi] <- 2 * pi - dLon[dLon > pi]

  r * sqrt(dLat * dLat + q * q * dLon * dLon)
}

#' @rdname spherical_distance
#' @export
dist_smallcircle <- function(p1, p2, sm, r = earth_radius(), method = c("haversine", "haversine2", "vincenty")) {
  stopifnot(is.numeric(r))
  dist_greatcircle2(p1, p2, r, method) * tectonicr::sind(sm) # sind(sm) == cosd(90-sm)
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
#' greatcircle_distance(berlin, tokyo)
greatcircle_distance <- function(a, b) {
  if (is.null(dim(a)) == T & is.null(dim(b))) { # a and b are vectors with one latitude/longitude
    delta <- orthodrome_haversine2(a, b)
  }
  if (!is.null(dim(a)) == T & !is.null(dim(b))) { # a and b are data.frames or vectors with more latitude/longitude
    delta <- c()
    for (i in 1:nrow(a)) {
      delta[i] <- orthodrome_haversine2(c(a$lat[i], a$lon[i]), c(b$lat[i], b$lon[i]))
    }
  }
  if (!is.null(dim(a)) == T & is.null(dim(b))) { # a is data.frames or vectors with more latitude/longitude and b is a vectors with one latitude/longitude
    delta <- c()
    for (i in 1:nrow(a)) {
      delta[i] <- orthodrome_haversine2(c(a$lat[i], a$lon[i]), b)
    }
  }
  return(delta)
}


gc_dist <- function(lat1, lon1, lat2, lon2) {
  acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1))
}



#' Misfit of Euler pole from azimuths of small-circle structures
#'
#' @param lat,lon,azi Latitude, longitude and azimuth of the small-circle structures (e.g. fracture zones along the crest of the sea-floor ridges), numeric in degrees
#' @param ep_lat,ep_lon estimated position of the Euler pole for relative plate motion, numeric in degrees
#' @details The actual search is made through a  grid centered on the original ep_lat and ep_lon
#' The process is repeated until the point is located to better than 1 degree.
#' @references Le Pichon, X. (1968), Sea-floor spreading and continental drift, J. Geophys. Res., 73( 12), 3661-3697, doi:10.1029/JB073i012p03661.
misfit_euler_from_azimuths <- function(lat, lon, azi, ep_lat, ep_lon) {
  data <- data.frame(lat, lon, azi)
  ep <- tectonicr::euler_pole(ep_lat, ep_lon, geo = TRUE)
  tazi <- tectonicr::model_shmax(data, ep) %>% pull(sc)

  sum(
    mapply(
      FUN = function(azi, azi_t) {
        (azi - azi_t)^2
      },
      azi = azi,
      azi_t = tazi
    )
  )
}

#' Misfit of Euler pole from spreading rates
#'
#' @param v spreading rate measured along a perpendicular to the crest of the ridge
#' @param lat,lon,azi Latitude, longitude and azimuth of the small-circle structures (e.g. fracture zones along the crest of the sea-floor ridges), numeric in degrees
#' @param ep_lat,ep_lon estimated position of the Euler pole for relative plate motion, numeric in degrees
#' @references Le Pichon, X. (1968), Sea-floor spreading and continental drift, J. Geophys. Res., 73( 12), 3661-3697, doi:10.1029/JB073i012p03661.
misfit_euler_from_rates <- function(lat, lon, azi, v, ep_lat, ep_lon) {
  # estimated euler pole
  ep <- tectonicr::euler_pole(ep_lat, ep_lon, geo = TRUE)

  # theoretical azimuths at e
  tazi <- tectonicr::model_shmax(data, ep) %>% pull(sc)

  p <- cbind(lat = lat, lon = lon, azi = azi, tazi = tazi) * pi / 180
  e <- c(lat = ep_lat, lon = ep_lon) * pi / 180

  # distances to e
  theta <- mapply(
    FUN = gc_dist,
    lat1 = p[, "lat"],
    lat2 = e["lat"],
    lon1 = p[, "lon"],
    lon2 = e["lon"]
  )

  tv_max <- v / sin(theta) # velocity at rotation's equator

  mapply(
    FUN = function(v, vmax, azi, tazi, sprn) {
      (v / (vmax * tectonicr::cosd(azi - tazi)) - sprn)^2
    },
    v = v,
    azi = p[, "azi"],
    tazi = p[, "tazi"],
    vmax = tv_max,
    sprn = sin(theta)
  )
}

morgan_alpha <- function(lon, eplon, eplat, theta) {
  asin(
    sin(eplon - lon) * cos(eplat) / sin(theta)
  )
}

#' Velocities at perpendicular along sea-floor ridge
#'
#' @param lat,lon,azi numeric. Latitude, longitude and strike of the great-circle structure (e.g. crest of the sea-floor ridges), in degrees
#' @param eplat,eplon numeric. Position of the Euler pole for relative plate motion (in degrees)
#' @param vmax numeric. Maximum velocity at equator of rotation around Euler pole
#' @references Morgan, W. J. (1968), Rises, trenches, great faults, and crustal blocks, J. Geophys. Res., 73( 6), 1959-1982, doi:10.1029/JB073i006p01959.
#' @examples
#' lat <- c(30, 18.5, 14.5, 10.8, 10.2, 9.4, 8.8, 7.6, 7.2, 4.0, 1.9, 1.1, -1.1, -0.1, -1.3, -1.9, -2.9, -7.5)
#' lon <- c(-42.3, -46.8, -46, -42.3, -40.9, -40, 38.7, -36.6, -34.3, -31.9, -30.6, -26, -24, -18, -14.5, -12.9, -12.5, -12.3)
#' strike <- c(99, 95, 91, 92, 94, 92, 92, 91, 91, 88, 86, 86, 81, 77, 75, 82, 73, 73)
#' ep_lat <- 62
#' ep_lon <- -36
#' vmax <- 1.8 # cm/yr
#' perpendicular_velocities(lat, lon, strike, ep_lat, ep_lon, vmax)
perpendicular_velocities <- function(lat, lon, strike, ep_lat, ep_lon, vmax) {
  p <- cbind(lat = lat, lon = lon, strike = strike) * pi / 180
  e <- c(lat = ep_lat, lon = ep_lon) * pi / 180

  theta <- mapply(
    FUN = gc_dist,
    lat1 = p[, "lat"],
    lat2 = e["lat"],
    lon1 = p[, "lon"],
    lon2 = e["lon"]
  )

  alpha <- mapply(
    FUN = morgan_alpha,
    lon = p[, "lon"],
    eplon = e["lon"],
    eplat = e["lat"],
    theta = theta
  )

  vmax * sin(theta) * cos(p[, "strike"] - alpha)
}

# Atlantis fracture zone
f1_lat = 31
f1_lon = -48
f2_lat = 29
f2_lon = -37
length = 1000

misfit_euler_from_line <- function(f1_lat, f1_lon, f2_lat, f2_lon, length, ep_lat, ep_lon, treshold=5){
  f1 <- c(lat = f1_lat, lon = f1_lon) * pi / 180
  f2 <- c(lat = f2_lat, lon = f2_lon) * pi / 180
  e <- c(lat = ep_lat, lon = ep_lon) * pi / 180

  theta1 <- gc_dist(
    lat1 = f1["lat"],
    lat2 = e["lat"],
    lon1 = f1["lon"],
    lon2 = e["lon"]
  )

  theta1 <- gc_dist(
    lat1 = f2["lat"],
    lat2 = e["lat"],
    lon1 = f2["lon"],
    lon2 = e["lon"]
  )

  c(theta1, theta2, theta1-theta2)
}

