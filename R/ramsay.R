#' The cone or plane best fit of conically or cylindrical disposed s-plane poles
#'
#' Finding the best fit pole of rotation for a given set of points that are
#' assumed to lie on a mutual small or great circle circle
#' @inheritParams mean_pole
#' @importFrom tectonicr deg2rad rad2deg
#' @importFrom dplyr mutate summarise
#' @name best_pole
#' @source Ramsay, 1967, p. 18-21
#' @examples
#' # example from Ramsay, 1967, p. 20
#' x <- rbind(
#'   c(-67, -31, -71),
#'   c(-62, -53, -50),
#'   c(-62, -75, -34),
#'   c(-58, 85, -34),
#'   c(-79, 40, -52),
#'   c(90, 14, -75),
#'   c(80, 10, 90)
#' ) |> data.frame()
#' colnames(x) <- c("a", "b", "c")
#' best_cone(x)
#' best_plane(x)
NULL

#' @rdname best_pole
#' @export
best_cone <- function(x) {
  # x <- tectonicr::rad2deg(x)
  x2 <- dplyr::mutate(x,
                      l = tectonicr:::cosd(a),
                      m = tectonicr:::cosd(b),
                      n = tectonicr:::cosd(c),
                      l = ifelse(a < 0, -l, l),
                      m = ifelse(b < 0, -m, m),
                      n = ifelse(c < 0, -n, n),
                      l2 = l^2,
                      m2 = m^2,
                      lm = l * m,
                      ln = l * n,
                      mn = m * n
  ) |>
    dplyr::summarise(
      sum_l = sum(l),
      sum_m = sum(m),
      sum_n = sum(n),
      sum_l2 = sum(l2),
      sum_m2 = sum(m2),
      sum_lm = sum(lm),
      sum_ln = sum(ln),
      sum_mn = sum(mn)
    )
  N <- nrow(x2)

  # Rule of Sarrus:
  D <- Da <- Db <- Dc <- matrix(nrow = 3, ncol = 3)
  D[1, 1] <- x2$sum_l2
  D[1, 2] <- D[2, 1] <- x2$sum_lm
  D[1, 3] <- D[3, 1] <- x2$sum_l
  D[2, 2] <- x2$sum_m2
  D[2, 3] <- D[3, 2] <- x2$sum_m
  D[3, 3] <- N

  Da[1, 1] <- -x2$sum_ln
  Da[1, 2] <- D[1, 2]
  Da[1, 3] <- D[1, 3]
  Da[2, 1] <- -x2$sum_mn
  Da[2, 2] <- D[2, 2]
  Da[2, 3] <- Da[3, 2] <- D[2, 3]
  Da[3, 1] <- -x2$sum_n
  Da[3, 3] <- N

  Db[1, 1] <- D[1, 1]
  Db[1, 2] <- Da[1, 1]
  Db[1, 3] <- D[1, 3]
  Db[2, 1] <- D[1, 2]
  Db[2, 2] <- Da[2, 1]
  Db[2, 3] <- D[2, 3]
  Db[3, 1] <- D[1, 3]
  Db[3, 2] <- Da[3, 1]
  Db[3, 3] <- N

  Dc[1, 1] <- D[1, 1]
  Dc[1, 2] <- D[1, 2]
  Dc[1, 3] <- Da[1, 1]
  Dc[2, 1] <- D[1, 2]
  Dc[2, 2] <- D[2, 2]
  Dc[2, 3] <- Da[2, 1]
  Dc[3, 1] <- D[1, 3]
  Dc[3, 2] <- D[2, 3]
  Dc[3, 3] <- Da[3, 1]
  # sum(l^2)*sum(m^2)*N + sum(l*m)*sum(m)*sum(l) + sum(l)*sum(l*m)*sum(m) - sum(l^2)*sum(m)^2 - N * sum(l*m)^2 - sum(m^2)*sum(l)^2

  A <- det(Da) / det(D)
  B <- det(Db) / det(D)
  C <- det(Dc) / det(D)

  gamma <- -tectonicr:::acosd((1 + A^2 + B^2)^(-1 / 2))
  alpha <- -tectonicr:::acosd(A * (1 + A^2 + B^2)^(-1 / 2))
  beta <- -tectonicr:::acosd(B * (1 + A^2 + B^2)^(-1 / 2))
  K <- tectonicr:::acosd(C * (1 + A^2 + B^2)^(-1 / 2)) |>
    tectonicr::deviation_norm()

  return(c("alpha" = alpha + 180, beta = beta, gamma = gamma, K = K))
}

#' @rdname best_pole
#' @export
best_plane <- function(x) {
  x2 <- dplyr::mutate(x,
                      l = tectonicr:::cosd(a),
                      m = tectonicr:::cosd(b),
                      n = tectonicr:::cosd(c),
                      l = ifelse(a < 0, -l, l),
                      m = ifelse(b < 0, -m, m),
                      n = ifelse(c < 0, -n, n),
                      l2 = l^2,
                      m2 = m^2,
                      lm = l * m,
                      ln = l * n,
                      mn = m * n
  )
  A <- sum(x2$lm) * sum(x2$mn) - sum(x2$ln) * sum(x2$m2) /
    (sum(x2$l2) * sum(x2$m2) - sum(x2$lm)^2)
  B <- sum(x2$lm) * sum(x2$ln) - sum(x2$mn) * sum(x2$l2) /
    (sum(x2$l2) * sum(x2$m2) - sum(x2$lm)^2)
  z <- 1 / sqrt((1 + A^2 + B^2))
  gamma <- -acos(z)
  alpha <- -acos(A * z)
  beta <- -acos(B * z)

  check <- cos(alpha)^2 + cos(beta)^2 + cos(gamma)^2

  return(
    c(
      tectonicr::rad2deg(
        c(alpha = alpha + pi, beta = beta, gamma = gamma)
      ),
      R = check
    )
  )
}

#' Unimodal pole distribution
#'
#' mean vector, variance and standard deviation
#'
#' @param x `data.frame` containing the three dimensions of stereographic
#' coordinates of the unit vectors in degree:
#' \describe{
#' \item{\code{a}}{angular deviation from x-axis E-W horizontal}
#' \item{\code{b}}{angular deviation from y-axis N-S horizontal}
#' \item{\code{c}}{angular deviation from z-axis vertical}
#' }
#' @importFrom tectonicr deg2rad rad2deg
#' @importFrom dplyr mutate
#' @source Ramsay, 1967, p. 15-17
#' @name pole_dist
#' @examples
#' # example from Ramsay, 1967, p. 17
#' df <- rbind(
#'   c(58, 78, 35),
#'   c(47, 72, 47),
#'   c(57, 68, 41),
#'   c(56, 64, 45),
#'   c(52, 65, 47),
#'   c(61, 62, 51),
#'   c(55, 60, 49),
#'   c(61, 62, 42),
#'   c(68, 52, 46),
#'   c(42, 58, 65)
#' ) |> as.data.frame()
#' colnames(df) <- c("a", "b", "c")
#' mean_pole(df)
#' var_pole(df)
#' sd_pole(df)
NULL

#' @rdname pole_dist
#' @export
mean_pole <- function(x) {
  cos_a <- tectonicr:::cosd(x$a)
  cos_b <- tectonicr:::cosd(x$b)
  cos_c <- tectonicr:::cosd(x$c)

  z <- sqrt(sum(cos_a)^2 + sum(cos_b)^2 + sum(cos_c)^2)
  cos_a <- sum(cos_a) / z
  cos_b <- sum(cos_b) / z
  cos_c <- sum(cos_c) / z

  tectonicr:::acosd(c(a = cos_a, b = cos_b, c = cos_c))
}
#' @rdname pole_dist
#' @export
var_pole <- function(x) {
  cos_a <- tectonicr:::cosd(x$a)
  cos_b <- tectonicr:::cosd(x$b)
  cos_c <- tectonicr:::cosd(x$c)

  sum(cos_a)^2 + sum(cos_b)^2 + sum(cos_c)^2
}
#' @rdname pole_dist
#' @export
sd_pole <- function(x) {
  sqrt(var_pole(x))
}


#' Spherical coordinate conversion
#'
#' @param x matrix of  coordinates
#' @name ramsay_coords
#' @examples
#' Stereographic coordinates (angle notation):
#' x <- rbind(
#'   c(58, 78, 35),
#'   c(47, 72, 47),
#'   c(57, 68, 41),
#'   c(56, 64, 45),
#'   c(52, 65, 47),
#'   c(61, 62, 51),
#'   c(55, 60, 49),
#'   c(61, 62, 42),
#'   c(68, 52, 46),
#'   c(42, 58, 65)
#' )
#' acoscartesian_to_cartesian(x)
#' acoscartesian_to_geographical(x)
#' acoscartesian_to_cartesian(x) |> cartesian_to_acoscartesian()
NULL

#' @rdname ramsay_coords
cartesian_to_acoscartesian <- function(x) {
  a <- acos(x[, 1])
  b <- acos(x[, 2])
  c <- acos(x[, 3])

  tectonicr::rad2deg(cbind(a, b, c))
}

#' @rdname ramsay_coords
acoscartesian_to_cartesian <- function(x) {
  x <- tectonicr::deg2rad(x)
  cx <- cos(x[, 1])
  cy <- cos(x[, 2])
  cz <- cos(x[, 3])
  cbind(cx, cy, cz)
}

#' @rdname ramsay_coords
geographical_to_cartesian2 <- function(x) {
  stopifnot(is.numeric(x))
  x <- tectonicr::deg2rad(x)
  cx <- cos(x[, 1]) * cos(x[, 2])
  cy <- cos(x[, 1]) * sin(x[, 2])
  cz <- sin(x[, 1])
  cbind(cx, cy, cz)
}

#' @rdname ramsay_coords
cartesian_to_geographical2 <- function(x) {
  stopifnot(is.numeric(x))
  r <- sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2)
  lat <- asin(x[, 3] / r)
  lon <- atan2(x[, 2], x[, 1])
  tectonicr::rad2deg(cbind(lat, lon))
}

#' @rdname ramsay_coords
geographical_to_acoscartesian <- function(x) {
  cartesian_to_acoscartesian(
    geographical_to_cartesian2(x)
  )
}

#' @rdname ramsay_coords
acoscartesian_to_geographical <- function(x) {
  cartesian_to_geographical2(
    acoscartesian_to_cartesian(x)
  )
}


ep_from_sc <- function(x) {
  x_sph <- geographical_to_acoscartesian(x)
  res <- best_cone(as.data.frame(x_sph))

  coords <- acoscartesian_to_geographical(t(as.matrix(res)))
  names(coords) <- c("lat", "lon")
  angle <- res[4]
  names(angle) <- ("angle")

  c(coords, angle)
}

ep_from_gc <- function(x) {
  x_sph <- geographical_to_acoscartesian(x)
  res <- best_plane(as.data.frame(x_sph))
  coords <- acoscartesian_to_geographical(t(as.matrix(res)))
  names(coords) <- c("lat", "lon")
  R <- res[4]
  c(coords, R)
}

#' Euler pole solution for geological structures
#'
#' find the best Euler pole solution for a set of geographic locations that
#' are aligned along a geological structure.
#'
#' @inheritParams to_geomat
#' @param sm logical. Are points `x` aligned on a small circle (`TRUE`) or great circle (`FALSE`)?
#' @importFrom sf st_coordinates
#' @return three-element vector given the latitude, longitude, and apical
#' angle of the best fit Euler pole.
#' @export
#' @examples
#' # Example from Price & Carmicheal (1986), GEOLOGY:
#'
#' # from matrix
#' rmt_mat <- rbind(
#'   "yukon" = c(66.1, -147.8),
#'   "bigbend" = c(52.25, -122.65),
#'   "washington" = c(47.85, -121.85)
#' )
#' # from data.frame
#' euler_solution(rmt_mat)
#'
#' rmt_df <- data.frame("lat" = rmt_mat[, 1], "lon" = rmt_mat[, 2])
#' euler_solution(rmt_df)
#'
#' # from sf object
#' rmt_sf <- sf::st_as_sf(rmt_df, coords = c("lon", "lat"), crs = "WGS84")
#' euler_solution(rmt_sf)
euler_solution <- function(x, sm = TRUE) {
  x_coords <- to_geomat(x)

  if (sm) {
    ep_from_sc(x_coords)
  } else {
    ep_from_gc(x_coords)
  }
}


#' Position statistics
#'
#' Statistics on the distribution of geographic locations
#'
#' @inheritParams to_geomat
#' @return three-element vector giving the mean geographic location, variance, and standard deviation
#' @export
geo_distribution <- function(x) {
  x_coords <- to_geomat(x)

  x_sph <- geographical_to_acoscartesian(x_coords)
  meanpole <- mean_pole(as.data.frame(x_sph))
  coords <- acoscartesian_to_geographical(t(as.matrix(meanpole)))
  names(coords) <- c("lat", "lon")

  varpole <- var_pole(as.data.frame(x_sph))
  sdpole <- sd_pole(as.data.frame(x_sph))
  c(coords, var = varpole, sd = sdpole)
}



#' Convert object to geographic coordinate matrix
#'
#' @param x a `matrix` or `data.frame` containing the geographical coordinates,
#' or a `sf` object
#'
#' @return matrix
to_geomat <- function(x){
  is.sf <- inherits(x, "sf")
  is.df <- is.data.frame(x)
  is.mat <- is.matrix(x)
  stopifnot(any(is.sf, is.df, is.mat))

  if(is.sf){
    x_coords <- sf::st_coordinates(x)
    cbind(x_coords[, 2], x_coords[, 1]) # switch columns
  } else if(is.df){
    cbind(x$lat, x$lon) # switch columns
  } else {
    x
  }
}

#
# library(ggplot2)
#
# rmt.res <- euler_solution(rmt)
# ep <- data.frame(lat = rmt.res[1], lon = rmt.res[2])
#
# ggplot()+
#   geom_sf(data = rmt) +
#   geom_sf(data = tectonicr::eulerpole_smallcircles(ep, 50)) +
#   coord_sf(xlim = c(sf::st_bbox(rmt)[1], sf::st_bbox(rmt)[3]), ylim = c(sf::st_bbox(rmt)[2], sf::st_bbox(rmt)[4]))
