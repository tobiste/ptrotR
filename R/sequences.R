#' @title Read GPLATES' rotation file
#' @description Imports a sequence of of total reconstruction rotations from a
#' GPLATES' .rot file
#' @param x either a file name (.rot format) OR a matrix
#' @param ... optional arguments to the read.table function
#' @return object of class \code{"finite"}
#' \describe{
#'   \item{plate.rot}{ID of moving plate}
#'   \item{lat}{Latitude of Euler pole of total reconstruction rotation}
#'   \item{lon}{Longitude}
#'   \item{angle}{Rotation angle in degree}
#'   \item{plate.fix}{ID of fixed/anchored plate}
#'   \item{cmt}{Comments}
#' }
#' @details The comment column (last column) must not include **white space**. Use "_" to separate words instead.
#' @seealso \code{\link{check.finite}}
#' @importFrom utils read.table
#' @export
#' @examples
#' fname <- system.file("Pangea.rot", package="ptrotR")
#' Pangea <- read.gplates(fname)
#' print(Pangea)
read.gplates <- function(x, ...) {
  data <- utils::read.table(file = x, header = FALSE, sep = "", dec = ".")
  colnames(data) <- c(
    "plate.rot", "age", "lat", "lon", "angle", "plate.fix",
    "sep", "cmt"
  )
  data$sep <- NULL

  return(data)
}

#' @title Finite rotation object
#'
#' @description Check if object has columns of a finite rotation sequence
#'
#' @param x data.frame containing columns lat, lon, angle, plate.rot, and plate.fix
#' \describe{
#'   \item{plate.rot}{ID of moving plate}
#'   \item{age}{Age for finite rotation}
#'   \item{lat}{Latitude of Euler pole of total reconstruction rotation}
#'   \item{lon}{Longitude}
#'   \item{angle}{Rotation angle in degree}
#'   \item{plate.fix}{ID of fixed/anchored plate}
#'   }
#' @return data.frame
#' @export
#' @examples
#' data(pangea)
#' check.finite(pangea)
check.finite <- function(x){
  if(!('plate.fix' %in% colnames(x))){
    stop("column 'plate.fix'is missing")
  }
  if(!('plate.rot' %in% colnames(x))){
    stop("column 'plate.rot' is missing")
  }
  if(!('angle' %in% colnames(x))){
    stop("column 'angle' is missing")
  }
  if(!('lat' %in% colnames(x))){
    stop("column 'lat' is missing")
  }
  if(!('lon' %in% colnames(x))){
    stop("column 'lon' is missing")
  }
  if(!('age' %in% colnames(x))){
    stop("column 'age' is missing")
  }
  return(TRUE)
}

#' @title Stage rotation object
#'
#' @description  Check if object has columns of a 'stage' rotation sequence
#'
#' @param x data.frame containing columns lat, lon, angle, plate.rot, and plate.fix
#' \describe{
#'   \item{plate.rot}{ID of moving plate}
#'   \item{lat}{Latitude of Euler pole of total reconstruction rotation}
#'   \item{lon}{Longitude}
#'   \item{min.age}{End of stage rotation}
#'   \item{max.age}{Start of stage rotation}
#'   \item{angle}{Rotation angle in degree}
#'   \item{plate.fix}{ID of fixed/anchored plate}
#'   }
#' @return data.frame
#' @export
#' @examples
#' data(pangea)
#' check.stage(extract_stage_rotations(pangea, plate=103))
check.stage <- function(x){
  if(!('plate.fix' %in% colnames(x))){
    stop("column 'plate.fix'is missing")
  }
  if(!('plate.rot' %in% colnames(x))){
    stop("column 'plate.rot' is missing")
  }
  if(!('angle' %in% colnames(x))){
    stop("column 'angle' is missing")
  }
  if(!('lat' %in% colnames(x))){
    stop("column 'lat' is missing")
  }
  if(!('lon' %in% colnames(x))){
    stop("column 'lon' is missing")
  }
  if(!('max.age' %in% colnames(x))){
    stop("column 'max.age' is missing")
  }
  if(!('min.age' %in% colnames(x))){
    stop("column 'min.age' is missing")
  }
  return(TRUE)
}


#' @title Antipodal rotation
#' @description Euler pole on the other side of the hemisphere
#' @param x data.frame containing the sequence of rotations or the rotation
#' @return object of with same class like x
#' @importFrom tectonicr longitude_modulo
#' @export
#' @examples
#' euler.pole <- data.frame(plate.rot = c("A", "B"), lat=c(27, -21), lon=c(17, -151), angle=c(0.4, 0.42), plate.fix = "C")
#' antipodal_rotation(euler.pole)
antipodal_rotation <- function(x) {
  for(i in seq_along(x$lat)) {
    x$lat[i] <- -x$lat[i]
    x$lon[i] <- tectonicr::longitude_modulo(x$lon[i] + 180)
    x$angle[i] <- -x$angle[i]
  }
  return(x)
}



#' @title Invert rotation
#' @description Changes plate motion (A relative to B) to (B relative to A)
#' @param x data.frame containing the sequence of rotations or the rotation
#' @return object of with same class like x
#' @export
#' @seealso \code{\link{check.finite}}, \code{\link{check.stage}}
#' @examples
#' x <- data.frame(plate.rot = "A", lat = 10, lon = -70, angle = 20, plate.fix = "B")
#' inverse_rotation(x)
inverse_rotation <- function(x){
  x.rev <- x
  x.rev$plate.rot <- x$plate.fix
  x.rev$plate.fix <- x$plate.rot
  x.rev$angle <- -x$angle
  return(x.rev)
}


eulerpole_2_eulervec <- function(x){
  v <- c(x$x, x$y, x$z,
         x$angle)
  class(v) <- append(class(v), "euler")
  return(v)
}

inv_euler <- function(x){
  x[4] <- -x[4]
  return(x)
}


#' @title Stage rotation extraction of rotation matrices
#' @description extract stage rotation of two finite rotation matrices
#' @param a1,a2 rotation matrix of finite rotations (t(a2)<t(a1))
#' @return list
#' @details x must be all equivalent total rotations.
#' @references Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. Computers and Geosciences, 25(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5
#' @export
#' @seealso \code{\link{extract_stage_rotations}}, \code{\link{extract_stage_rotation}}
extract_stage_rotation_matrix <- function(a1, a2) {
  a12 <- euler::matrix_2_angles(a1 %*% solve(a2))
  return(a12)
}

#' @title Stage rotation extraction of rotations
#' @description extract stage rotation of two finite rotations
#' @param r1,r2 rotation matrix of finite rotations (t(r1)<t(r2))
#' @return list
#' @details x must be all equivalent total rotations.
#' @references
#' Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of Tectonic
#' Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg (Eds.),
#' *Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series*
#' (pp. 1--7). Springer Nature Switzerland AG 2021.
#' \doi{10.1007/978-3-030-26050-7_435-1}
#' @seealso \code{\link{extract_stage_rotations}}, \code{\link{extract_stage_rotation_matrix}}
#' @examples
#' x <- data.frame(plate.rot = 101, age = c(400, 530), lat = c(72.72, 68.3319), lon = c(150, 11.5198), angle = c(25.59, -31.9464), plate.fix = 301)
#' extract_stage_rotation(x[1,], x[2, ])
extract_stage_rotation <- function(r1, r2){
  r01 <- euler::to_euler(c(r1$lat, r1$lon, r1$angle))
  r02 <- euler::to_euler(c(r2$lat, r2$lon, r2$angle))
  euler::relative_euler_schaeben2(r01, r02)
}


#' @title Stage rotations extraction
#' @description extract all stage rotations from a sequence of  total
#' reconstruction rotations
#' @param x data.frame. Sequence of total reconstruction rotations
#' @param plate ID of plate
#' @return data.frame. Sequence of stage rotations
#' @details x must ba all equivalent total rotations.
#' @references Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. Computers and Geosciences, 25(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5
#' @export
#' @seealso \code{\link{extract_stage_rotation}}, \code{\link{check.finite}}, \code{\link{check.stage}}
#' @examples
#' data(pangea)
#' extract_stage_rotations(pangea, plate=137)
extract_stage_rotations <- function(x, plate) {
  check.finite(x)

  data <- subset(x, x$plate.rot == plate)
  age.list <- unique(data$age)
  for (time in 2:length(age.list)) {
    rot.i <- subset(data, data$age == age.list[time - 1])
    rot.j <- subset(data, data$age == age.list[time])
    ep.ij <- extract_stage_rotation(rot.i, rot.j)
    if(is.nan(ep.ij$axis[1])){
      ep.ij$axis <- c(90, 0)
      ep.ij$angle <- 0
    }

    # all northern hemisphere poles:
    if (ep.ij$axis[1] < 0) {
      ep.ij$axis <- antipodal_euler_pole(c(ep.ij$axis[1], ep.ij$axis[2]))
      ep.ij$angle <- -ep.ij$angle
    }

    ep.ij <- data.frame(
      "lat" = ep.ij$axis[1],
      "lon" = ep.ij$axis[2],
      "angle" = ep.ij$angle,
      "max.age" = age.list[time],
      "min.age" = age.list[time - 1],
      "plate.rot" = data$plate.rot[1],
      "plate.fix" = data$plate.fix[1]
    )

    if (time == 2) {
      stage_poles <- ep.ij
    } else {
      stage_poles <- rbind(stage_poles, ep.ij)
    }
  }

  return(stage_poles)
}


#' @title Find missing rotations in a rotation sequence
#' @description Identifies gaps in a sequence of rotations with
#' different reference systems
#' @param x data.frame. Sequence of total
#' reconstruction rotations
#' @return data.frame. Sequence of total reconstruction rotations
#' @details x must ba all equivalent total rotations.
#' @importFrom dplyr first
#' @export
#' @seealso \code{\link{check.finite}}, \code{\link{check.stage}}
find_missing_rotations <- function(x) {
  age.list <- unique(x$age) # all unique ages

  for (i in unique(x$plate.rot)) { # loop through all moving plates
    x.i <- subset(x, x$plate.rot == i)
    for (age in age.list) { # loop through the age list
      if (!(age %in% x.i$age)) {
        missing.x.i <- data.frame("plate.rot" = i, "age" = age, missing = T)
      } else {
        missing.x.i <- data.frame("plate.rot" = i, "age" = age, missing = F)
      }
      if (age == first(age.list)) {
        missing.i.age <- missing.x.i
      } else {
        missing.i.age <- rbind(missing.i.age, missing.x.i)
      }
    }
    if (i == first(unique(x$plate.rot))) {
      missing.i <- missing.i.age
    } else {
      missing.i <- rbind(missing.i, missing.i.age)
    }
  }

  return(missing.i)
}




#' @title Last finite rotation in sequence
#' @description Adds last finite rotation in a sequence assuming that the plate motion does not change
#' @param x data.frame. Sequence of total reconstruction rotations
#' @return data.frame.
#' @export
#' @importFrom plyr rbind.fill
#' @examples
#' data(pangea)
#' add_last_finite_rotations(pangea)
add_last_finite_rotations <- function(x){
  t <- max(x$age)
  for(id in x$plate.rot){ # loop through all plates
    x.id <- subset(x, x$plate.rot == id)
    if(max(x.id$age) != t){
        x.id.oldest1 <- subset(x.id, x.id$age==max(x.id$age))
        x.id.oldest2 <- data.frame(
          plate.rot = x.id.oldest1$plate.rot,
          age = t,
          lat = x.id.oldest1$lat,
          lon = x.id.oldest1$lon,
          angle = x.id.oldest1$angle,
          plate.fix = x.id.oldest1$plate.fix,
          cmt = 'added_oldest_finite_rotation'
        )
        x <- plyr::rbind.fill(x, x.id.oldest2)
      }
    }
  return(x)
}


#' @title Interpolate gaps in the sequence of total reconstruction rotations
#' @description Interpolate missing rotations in a sequence of rotations with a
#' different reference systems
#' @param df data.frame. Sequence of total
#' reconstruction rotations
#' @return data.frame. Sequence of total
#' reconstruction rotations with filled gaps
#' @importFrom plyr rbind.fill
#' @export
#' @seealso \code{\link{check.finite}}, \code{\link{check.stage}}
#' @examples
#' data(pangea)
#' interpolate_missing_finite_poles(pangea)
interpolate_missing_finite_poles <- function(df) {
  check.finite(df)
  df <- add_last_finite_rotations(df)
  missing.df <- find_missing_rotations(df)
  missing.df <- missing.df[order(missing.df$plate.rot, missing.df$age),]

  missing.rot <- data.frame(
    plate.rot = as.character(),
    age = as.numeric(),
    lat = as.numeric(),
    lon = as.numeric(),
    angle = as.numeric(),
    plate.fix = as.character(),
    cmt = as.character()
  )

  df <- df[order(df$plate.rot, df$age),]

  for (id in unique(missing.df$plate.rot)) { # loop through all plates
    x <- subset(missing.df, missing.df$plate.rot == id)

    for (i in seq_along(x$missing)) { # loop through all finite rotations
      if (x$missing[i]) {
        missing.age <- x$age[i]

        if(missing.age==0){
          fixed <- subset(df, df$plate.rot == id)
          missing.rot <- rbind(
            missing.rot,
            data.frame(
              plate.rot = id,
              age = missing.age,
              lat = 90,
              lon = 0,
              angle = 0,
              plate.fix = unique(fixed$plate.fix),
              cmt = 'added_present-day'
            )
          )
          df <- unique(plyr::rbind.fill(df, missing.rot)) # add new rotation to rot file
          x$missing[i] <- FALSE # add tell that the rotation is not missing anymore

        } else {
          older.age <- subset(x, x$age > missing.age & !x$missing)
          # x %>% dplyr::filter(age > missing.age & missing == F)
          suppressWarnings({
            older.age <- min(older.age$age)
          })


          if (!is.infinite(older.age)) {
            younger.age <- subset(x, x$age < missing.age & !x$missing)
            # younger.age <- x %>% dplyr::filter(age < missing.age & missing == F)
            younger.age <- max(younger.age$age)

            fixed <- subset(df, df$plate.rot == id)
            fixed <- unique(fixed$plate.fix)

            # df.young <- df %>% dplyr::filter(plate.rot==id, age == younger.age)
            # df.old <- df %>% dplyr::filter(plate.rot==id, age == older.age)
            df.young <- subset(df, df$plate.rot == id & df$age == younger.age)
            df.old <- subset(df, df$plate.rot == id & df$age == older.age)

            if (df.young$lat == df.old$lat & df.young$lon == df.old$lon & df.young$angle == df.old$angle) {
              missing.rot <- rbind(
                missing.rot,
                data.frame(
                  plate.rot = id,
                  age = missing.age,
                  lat = df.young$lat,
                  lon = df.young$lon,
                  angle = df.young$angle,
                  plate.fix = fixed,
                  cmt = paste0(
                    "inteprolated_between_", older.age, "_Ma_and_",
                    younger.age, "_Ma"
                  )
                )
              )
              df <- unique(plyr::rbind.fill(df, missing.rot)) # add new rotation to rot file
              x$missing[i] <- FALSE # add tell that the rotation is not missing anymore


            } else {
              z <- finite_pole_interpolation(
                df.young,
                df.old,
                missing.age
              )

              missing.rot <- rbind(
                missing.rot,
                data.frame(
                  plate.rot = id,
                  age = missing.age,
                  lat = z$axis[1],
                  lon = z$axis[2],
                  angle = z$angle,
                  plate.fix = fixed,
                  cmt = paste0(
                    "inteprolated_between_", older.age, "_Ma_and_",
                    younger.age, "_Ma"
                  )
                )
              )
              df <- unique(plyr::rbind.fill(df, missing.rot)) # add new rotation to rot file
              x$missing[i] <- FALSE # add tell that the rotation is not missing anymore

            }
          }
        }
      }
    }
  }
  return(df)
}


#' @title Interpolate interjacent rotation between two total reconstruction rotations
#'
#' @param rot1,rot2 data.frame. \code{rot1$age} < \code{rot2$age}. Both have same fixed plate.
#' @param tx number. Age of the desired interjacent finite rotation
#' Must be in between \code{rot1$age} and \code{rot2$age}
#' @return list. Sequence of total
#' reconstruction rotations with filled gaps
#' @references Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. Computers and Geosciences, 25(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5
#' @export
#' @seealso \code{\link{check.finite}}
#' @importFrom dplyr between
#' @importFrom tectonicr relative_rotation euler_pole rad2deg
#' @importFrom euler euler_concatenation
#' @examples
#' x <- data.frame(plate.rot = 101, age = c(400, 530), lat = c(72.72, 68.3319), lon = c(150, 11.5198), angle = c(25.59, -31.9464), plate.fix = 301)
#' finite_pole_interpolation(x[1, ], x[2, ], 410)
#' finite_pole_interpolation(x[1, ], x[2, ], 430)
#' all.equal(finite_pole_interpolation(x[1, ], x[2, ], 430), finite_pole_interpolation(x[2, ], x[1, ], 430)) # mixing up the inputs doesn't matter!
#'
#' y = data.frame(plate.rot = 101, age = c(0, 10), lat = c(90, 0), lon = c(0, 0), angle = c(1, 2), plate.fix = 301)
#' finite_pole_interpolation(y[1, ], y[2, ], 9)
finite_pole_interpolation <- function(rot1, rot2, tx) {
  check.finite(rot1)
  check.finite(rot2)
  if (rot1$plate.fix != rot2$plate.fix) {
    stop("Anchored/fixed plate of both total reconstruction rotations must be identical")
  }

  if (rot1$age == rot2$age) {
    stop("Age of finite rotations is identical")
  }

  if (rot1$age == rot2$age) {
    stop("Age of finite rotations is identical")
  }

  # rot1 must younger than rot2
  if (rot1$age > rot2$age) {
    rot1.copy <- rot1
    rot1 <- rot2
    rot2 <- rot1.copy
  }

  if (!dplyr::between(tx, rot1$age, rot2$age)) {
    stop("Interjacent time lies not in between the finite rotations")
  }

  # rotation from t=0 to t=1
  ROT_t01 <- euler::to_euler(c(rot1$lat, rot1$lon, rot1$angle))

  # rotation from t=0 to t=2
  ROT_t02 <- euler::to_euler(c(rot2$lat, rot2$lon, rot2$angle))

  # stage pole between rot1 and rot2
  rot12 <- euler::relative_euler_schaeben2(ROT_t02, ROT_t01) # R01 * (R02)^-1

  # rotation from t=1 to t=tx
  ROT_t1x <- euler::to_euler(c(
          rot12$axis[1],
          rot12$axis[2],
          rot12$angle / (rot2$age - rot1$age) * (tx - rot1$age)
          ))

  # interjacent finite rotation is rotation composition of rot1 and stage rotation rot12
  euler::euler_concatenation(ROT_t01, inv_euler(ROT_t1x)) # (R1x)^-1 * R01
  }


#' @title Equivalent rotations of different reference system
#' @description Transforms a sequence of rotations with different reference
#' systems into one with a common reference system
#' @param x sequence of plate rotations. An object of \code{"data.frame"}, \code{"tibble"}, or \code{"matrix"}
#' @param fixed ID of fixed plate. Must be same object class as x
#' @return sequence of plate rotations. Same object class as x
#' @importFrom dplyr %>% mutate
#' @importFrom tectonicr equivalent_rotation
#' @export
#' @seealso \code{\link{check.finite}}, \code{\link{check.stage}}
#' @examples
#' data(pangea)
#' equivalent_rotations(pangea, fixed=103)
equivalent_rotations <- function(x, fixed) {
  if(missing(fixed)){
    stop('Argument "fixed" is missing')
  }

  x.compl <- interpolate_missing_finite_poles(x)
  x.eq <- x[0, ] # create blank

  v <- unique(x.compl$age) # list of all ages
  v <- v[v != 0] # exclude age = 0
  for (i in v) {
    x.age <- subset(x.compl, x.compl$age == i)
    if (nrow(x.age) > 1) {
      x.age.eq <- tectonicr::equivalent_rotation(x.age, fixed) %>%
        dplyr::mutate(age = i, cmt = paste0(x.age$cmt, "_+_transformed"))
      x.eq <- rbind(x.eq, x.age.eq)
    }
  }
  return(x.eq)
}


#' @title Euler pole migration rate
#' @description Calculates the velocity and magnitude of a migrating Euler pole from a sequence of stage rotations
#' @param x data.frame. Sequence of stage rotations
#' @return data.frame
#' @importFrom dplyr "%>%" first mutate lag group_by
#' @export
#' @examples
#' data(pangea)
#' stages <- extract_stage_rotations(pangea, plate=103)
#' eulerpole_migration(stages)
eulerpole_migration <- function(x) {
  plate.rot <- ep.migration.vel <- min.age <- NULL
  check.stage(x)

  x$antipodal <- ifelse(x$lat > 0, TRUE, FALSE)

  for (p in unique(x$plate.rot)) {
    df.p <- subset(x, x$plate.rot == p | x$antipodal == TRUE)
    t <- df.p$max.age - df.p$min.age
    for (i in 1:nrow(df.p)) {
      if (i == nrow(df.p)) {
        df.p$ep.migration[i] <- NA
        df.p$ep.migration.vel[i] <- NA
      } else {
        ep.jump.p1 <-
          greatcircle_distance(c(df.p$lat[i], df.p$lon[i]), c(df.p$lat[i + 1], df.p$lon[i +
                                                                                1]))
        antipodal <-
          antipodal_euler_pole(c(df.p$lat[i + 1], df.p$lon[i + 1]))
        ep.jump.p2 <-
          greatcircle_distance(c(df.p$lat[i], df.p$lon[i]),
                     c(antipodal[1], antipodal[2]))
        df.p$ep.migration[i] <-
          ifelse(ep.jump.p1 <= ep.jump.p2, ep.jump.p1, ep.jump.p2)
        df.p$ep.migration.vel[i] <- df.p$ep.migration[i] / t[i]
      }
    }
    if (p == dplyr::first(unique(x$plate.rot))) {
      data2 <- df.p
    } else {
      data2 <- rbind(data2, df.p)
    }
    data2 <- data2 %>%
      dplyr::group_by(plate.rot) %>%
      dplyr::mutate(ep.migration.dvel = ep.migration.vel - dplyr::lag(
        ep.migration.vel,
        default = dplyr::first(ep.migration.vel),
        order_by = min.age
      ))
  }
  return(data2)
}


#' @title Create a grid
#' @description Create a grid of equally spaced points
#' @param gridsize  grid size in degree
#' @param lat.lim vector, range of latitudes
#' @param lon.lim vector, range of longitudes
#' @return data.frame
#' @export
grid_points <- function(gridsize, lat.lim, lon.lim){

  lats <- seq(min(lat.lim), max(lat.lim), gridsize)
  lons <- seq(min(lon.lim), max(lon.lim), gridsize)

  grid.df <- data.frame("lon" = as.numeric(),
                        "lat" = as.numeric())

  for(i in lons) {
    grid.df.i <- data.frame(
      "lat" = lats,
      "lon" = i)
    grid.df <- rbind(grid.df, grid.df.i)
  }
  return(grid.df)
}


#' @title Plate motion grid
#' @description Create a grid of the plate motion direction and velocity
#' @param euler  data.frame. containing lat, lon, and angle (optional) of Euler rotation
#' @param gridsize  grid size in degree
#' @param lat.lim vector, range of latitudes
#' @param lon.lim vector, range of longitudes
#' @return data.frame with plate motion direction at grid point
#' @importFrom tectonicr get_azimuth abs_vel
#' @export
#' @examples
#' data(pangea)
#' euler <- subset(pangea, pangea$plate.rot == 103 & pangea$age == 250)
#' plate_motion_grid(euler)
plate_motion_grid <- function(euler, gridsize = 5, lat.lim = c(-90, 90), lon.lim = c(-180, 180)){
  grid <- grid_points(gridsize, lat.lim, lon.lim)
  #r <-  6371.00887714

  for(i in 1:nrow(grid)){
    grid$azimuth[i] <- (
      tectonicr::get_azimuth(
        c(grid$lat[i], grid$lon[i]),
        c(euler$lat[1], euler$lon[1])) - 90
      ) %% 180

    grid$sm[i] <- greatcircle_distance(
      c(grid$lat[i], grid$lon[i]),
      c(euler$lat[1], euler$lon[1])
      )

    if(!is.null(euler$angle)){
      grid$velocity.abs[i] <- tectonicr::abs_vel(euler$angle[1], grid$sm[i])
      #grid$velocity.rel[i] <- pracma::rad2deg(grid$velocity.abs[i] / (r * pracma::sind(grid$sm[i])))
    }
  }
  return(grid)
}
