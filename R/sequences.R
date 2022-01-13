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
#' @details The comment column (last column) must not include
#' **white space**. Use "_" to separate words instead.
#' @importFrom utils read.table
#' @export
#' @examples
#' fname <- system.file("Pangea.rot", package="PlateTectonicMotionR")
#' Pangea <- read.gplates(fname)
#' print(Pangea)
read.gplates <- function(x, ...) {
  data <- utils::read.table(file = x, header = FALSE, sep = "", dec = ".")
  colnames(data) <- c(
    "plate.rot", "age", "lat", "lon", "angle", "plate.fix",
    "sep", "cmt"
  )
  data$sep <- NULL

  class(data) <- append(class(data), "finite")
  return(data)
}

#' @title Stage rotation extraction
#' @description extract stage rotation of two finite rotations
#' @param a1 rotation matrix of finite rotation 1
#' @param a2 rotation matrix of finite rotation 2
#' @return list
#' @details x must ba all equivalent total rotations.
#' @references <div class="csl-entry">Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. <i>Computers and Geosciences</i>, <i>25</i>(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5</div>
#' @export
#' @seealso \code{\link{extract_stage_rotations}}
extract_stage_rotation <- function(a1, a2) {
  a12 <- euler_from_rot(a1 %*% solve(a2))
  return(a12)
}


#' @title Stage rotations extraction
#' @description extract all stage rotations from a sequence of  total
#' reconstruction rotations
#' @param x object of class \code{"finite"}. Sequence of total
#' reconstruction rotations
#' @param plate ID of plate
#' @return object of class \code{"stage"}. Sequence of stage rotations
#' @details x must ba all equivalent total rotations.
#' @references <div class="csl-entry">Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. <i>Computers and Geosciences</i>, <i>25</i>(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5</div>
#' @export
#' @seealso \code{\link{extract_stage_rotation}}
#' @examples
#' data(pangea)
#' extract_stage_rotations(pangea, plate=103)
extract_stage_rotations <- function(x, plate) {
  data <- subset(x, x$plate.rot == plate)
  age.list <- unique(data$age)
  for (time in 2:length(age.list)) {
    rot.i <- subset(data, data$age == age.list[time - 1])
    finite.i <-
      euler_rot(euler_pole(rot.i$lat, rot.i$lon), rot.i$angle)

    rot.j <- subset(data, data$age == age.list[time])
    finite.j <-
      euler_rot(euler_pole(rot.j$lat, rot.j$lon), rot.j$angle)

    euler_pole_error <- NA
    ep.i <- tryCatch(
      expr = extract_stage_rotation(finite.j, finite.i),
      error = function(cond) {
        pole <- data.frame(
          lat = 0,
          lon = 0,
          X = 0,
          Y = 0,
          Z = 0
        )
        psi <- 0
        error.pole <- list(pole = pole, psi = psi)
        return(error.pole)
      }
    )

    # all northern hemisphere poles:
    if (ep.i$pole$lat < 0) {
      ep.i$pole <- antipodal_euler_pole(ep.i$pole)
      ep.i$psi <- ep.i$psi * -1
    }

    ep.i <- data.frame(
      "lat" = ep.i$pole$lat,
      "lon" = ep.i$pole$lon,
      "angle" = ep.i$psi,
      "max.age" = age.list[time],
      "min.age" = age.list[time - 1],
      "plate.rot" = data$plate[1],
      "plate.fix" = data$fixed[1]
    )

    if (time == 2) {
      stage_poles <- ep.i
    } else {
      stage_poles <- rbind(stage_poles, ep.i)
    }
  }

  class(stage_poles)[2] <- "stage"
  return(stage_poles)
}




#' @title Find missing rotations in a rotatin sequence
#' @description Identifies gaps in a sequence of rotations with
#' different reference systems
#' @param x object of class \code{"finite"}. Sequence of total
#' reconstruction rotations
#' @return object of class \code{"finite"}. Sequence of total reconstruction rotations
#' @details x must ba all equivalent total rotations.
#' @importFrom dplyr first
#' @export
find_missing_rotations <- function(x) {
  # identify missing rotations in a sequence of rotations with different reference systems
  # input df: dataframe consisting of lat, lon, ID, angle, and fixed
  # output df: dataframe showing the missing rotations

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


#' @title Interpolate gaps in the sequence of total reconstruction rotations
#' @description Interpolate missing rotations in a sequence of rotations with a
#' different reference systems
#' @param df object of class \code{"finite"}. Sequence of total
#' reconstruction rotations
#' @return object of class \code{"finite"}. Sequence of total
#' reconstruction rotations with filled gaps
#' @importFrom plyr rbind.fill
#' @export
#' #' @examples
#' data(pangea)
#' interpolate_missing_finite_poles(pangea)
interpolate_missing_finite_poles <- function(df) {
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

  for (id in unique(missing.df$plate.rot)) {
    x <- subset(missing.df, missing.df$plate.rot == id)

    for (i in seq_along(x$missing)) {
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
              cmt = 'present-day'
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
                  lat = z$pole$lat,
                  lon = z$pole$lon,
                  angle = z$psi,
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


#' @title Interpolate gaps in the sequence of total reconstruction rotations
#' @description Interpolate missing rotations in a sequence of rotations with a
#' different reference systems
#' @param rot1 object of class \code{"finite"}.
#' @param rot2 object of class \code{"finite"}. Must have the same fixed plate as rot1
#' @param tx number. Age of the requested intermediate finite rotation
#' Must be in between \code{rot1$age} and \code{rot2$age}
#' @return object of class \code{"finite"}. Sequence of total
#' reconstruction rotations with filled gaps
#' @references <div class="csl-entry">Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. <i>Computers and Geosciences</i>, <i>25</i>(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5</div>
#' @export
#' @importFrom dplyr between
finite_pole_interpolation <- function(rot1, rot2, tx) {
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
    stop("Intermediate time lies not in between the finite rotations")
  }

  # rotation matrix for rotation from t=0 to t=1
  ROT_t01 <-
    euler_rot(euler_pole(rot1$lat, rot1$lon), psi = rot1$angle)

  # rotation matrix for rotation from t=0 to t=2
  ROT_t02 <-
    euler_rot(euler_pole(rot2$lat, rot2$lon), psi = rot2$angle)

  # stage pole between rot1 and rot2
  rot12 <- extract_stage_rotation(ROT_t02, ROT_t01)

  # rotation matrix for rotation from t=1 to t=tx
  ROT_t1x <-
    euler_rot(
      rot12$pole,
      rot12$psi / (rot2$age - rot1$age) * (rot2$age - tx)
    )

  # intermediate finite rotation is rotation composition of rot1 and stage rotation rot12
  ROT0tx <- ROT_t1x %*% ROT_t01
  return(euler_from_rot(ROT0tx))
}


#' @title Equivalent rotation
#' @description Transform a sequence of rotations into a new reference system
#' @param x sequence of plate rotations. Object of class \code{"finite"} of \code{"stage"}
#' @param fixed ID of new fixed plate. Has to be one out of \code{x$plate.fix}
#' @return sequence of plate rotations in new reference system. Same object class as x
#' @export
equivalent_rotation <- function(x, fixed) {
  if(!(fixed %in% x$plate.rot)){
    stop("'fixed' has to be one out of x$plate.rot")
  }
  lat.eq <- c()
  lon.eq <- c()
  angle.eq <- c()

  fixed.plate <- subset(x, x$plate.rot == fixed) # 'fixed' has to be one out of x$plate.rot
  fixed.ep <- euler_pole(fixed.plate$lat, fixed.plate$lon)

  if (exists(paste0("fixed.plate$plate.fix"))) {
    fixed0.plate <- subset(x, x$plate.rot == fixed.plate$plate.fix)
    fixed0.ep <- euler_pole(fixed0.plate$lat, fixed0.plate$lon)
    fixed0.rot <- euler_rot(fixed0.ep, -fixed0.plate$angle)

    fixed.rot <-
      euler_rot(fixed.ep, fixed.plate$angle) %*% fixed0.rot
  } else {
    fixed.rot <- euler_rot(fixed.ep, -fixed.plate$angle)
  }

  for (i in seq_along(x$plate.rot)) {
    if (x$plate.rot[i] == fixed) {
      # fixed plate has no rotation
      lat.eq[i] <- 90
      lon.eq[i] <- 0
      angle.eq[i] <- 0
    } else {
      # Composition of finite rotations
      equivalent.rot <- fixed.rot %*%
        euler_rot(euler_pole(x$lat[i], x$lon[i]), x$angle[i])

      equivalent.ep <- euler_from_rot(equivalent.rot)

      lat.eq[i] <- equivalent.ep$pole$lat
      lon.eq[i] <- equivalent.ep$pole$lon
      angle.eq[i] <- equivalent.ep$psi
    }
  }
  x.eq <- data.frame(
    plate.rot = x$plate.rot,
    lat = lat.eq,
    lon = lon.eq,
    angle = angle.eq,
    plate.fix = fixed
  )
  return(x.eq)
}


#' @title Equivalent rotations of different reference system
#' @description Transform a sequence of rotations with different reference
#' systems into one with a common reference system
#' @param x sequence of plate rotations. An object of \code{"finite"} of \code{"stage"}
#' @param fixed ID of fixed plate. Must be same object class as x
#' @return sequence of plate rotations. Same object class as x
#' @importFrom dplyr %>% mutate
#' @export
#' @examples
#' data(pangea)
#' equivalent_rotations(pangea, fixed=103)
equivalent_rotations <- function(x, fixed) {
  x.compl <- interpolate_missing_finite_poles(x)
  x.eq <- x[0, ] # create blank

  v <- unique(x.compl$age) # list of all ages
  v <- v[v != 0] # exclude age = 0
  for (i in v) {
    x.age <- subset(x.compl, x.compl$age == i)
    if (nrow(x.age) > 1) {
      x.age.eq <- equivalent_rotation(x.age, fixed) %>%
        mutate(age = i, cmt = paste0(x.age$cmt, "_+_transformed"))
      x.eq <- rbind(x.eq, x.age.eq)
    }
  }
  return(x.eq)
}
