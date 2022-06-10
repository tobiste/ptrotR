

#' Function to test if an R object is a Spatial* or a sf data-type.
#'
#' Helper functions to
#'
#' @param x Object to be tested.
#'
#' @return returns \code{TRUE} if its argument is a Spatial* or a sf
#'  (that is, has \code{"Spatial*"} or \code{"sf"} amongst its classes) and
#'  FALSE otherwise.
#'
#' @name is_spatial
NULL

#' @rdname is_spatial
#' @export
is.sp <- function(x) {
  #if (grepl("spatial", class(x)[1], ignore.case = TRUE)) TRUE else FALSE
  grepl("spatial", class(x)[1], ignore.case = TRUE)
}

#' @rdname is_spatial
#' @export
is.sf <- function(x) {
  inherits(x, "sf")
}



#' Clip grid along plate boundaries
#'
#' Crops a grid along plate boundaries
#'
#' @param grid \code{"data.frame"}
#' @param plate an object of \code{"SpatialLinesDataFrame"} or \code{"sf"}
#' @return \code{"data.frame"}
#' @importFrom sf st_as_sf sf_use_s2 st_intersection
#' @export
clip_grid_by_plate <- function(grid, plate){
  grid.sf <- sf::st_as_sf(grid, coords = c("lon", "lat"), crs = "WGS84")

  if(is.sp(plate)){
    plate <- sf::st_as_sf(plate)
  }

  suppressMessages(sf::sf_use_s2(FALSE))
  grid.clipped <- sf::st_intersection(grid.sf, plate)
  suppressMessages(sf::sf_use_s2(TRUE))

  sf.as.data.frame(grid.clipped)
}



#' Clip Euler circles
#'
#' Clip small circles, great circles, or loxodromes at plate boundaries
#'
#' @param x Euler circle object. Either an object of \code{"sf"} or
#' \code{"SpatialLinesDataFrame"}
#' @param crop_by an object of \code{"SpatialLinesDataFrame"}, \code{"sf"} or a
#' \code{"data.frame"} with the range of longitudes and latitudes (lon, lat).
#' @return An object of \code{class(x)}
#' @importFrom sf st_intersection st_as_sf sf_use_s2 as_Spatial st_bbox st_crs
#' @export
clip_eulercircles <- function(x, crop_by){
  stopifnot(is.sp(x) | is.sf(x))
  stopifnot(is.sp(crop_by) | any(class(crop_by) == c("sf", "data.frame")))

  if (is.sp(x)) {
    x <- sf::st_as_sf(x)
  }
  if (is.sp(crop_by)) {
    crop_by <- sf::st_as_sf(crop_by)

  } else if (!is.sf(crop_by)) {
    crop_by <- sf::st_bbox(c(
      xmin = min(crop_by$lon),
      xmax = max(crop_by$lon),
      ymin = min(crop_by$lat),
      ymax = max(crop_by$lat)
    ),
    crs = sf::st_crs("WGS84")
    )
  }

  suppressMessages(sf::sf_use_s2(FALSE))
  x.clipped <- sf::st_intersection(x, crop_by)
  suppressMessages(sf::sf_use_s2(TRUE))

  if(is.sp(x)) {
    x.clipped <- x.clipped %>% sf::as_Spatial()
  }

  return(x.clipped)
}

#' Simple feature to data.frame
#'
#' Converts an \code{sf} object into a \code{data.frame}
#'
#' @param x an object of class \code{sf}
#' @return an object of class \code{data.frame}
#' @importFrom sf st_coordinates st_drop_geometry
#' @export
sf.as.data.frame <- function (x) {
  coords <- sf::st_coordinates(x)
  colnames(coords) <- c('lon', 'lat')
  cbind(coords, sf::st_drop_geometry(x))
}





#' @title Rose Diagram
#'
#' @description Plots a rose diagram (rose of directions), the analogue of a histogram or
#' density plot for angular data.
#'
#' @author Copyright (C) 2021 Tobias Stephan
#'
#' @param x Data to be plotted. A numeric vector containing angles, or a
#' histogram object containing a histogram of angular values, or a density
#' object containing a smooth density estimate for angular data, or an fv object
#' giving a function of an angular argument.
#' @param bins number of arcs to partition the circle with.
#' @param axial Logical value indicating whether data are uniaxial (axial=FALSE)
#' or biaxial (axial=TRUE, the default).
#' @param clockwise Logical value indicating whether angles increase in the
#' clockwise direction (clockwise=TRUE, the default) or anti-clockwise,
#' counter-clockwise direction (clockwise=FALSE).
#' @param unit The unit in which the angles are expressed. "degree" the default,
#' or "radian".
#' @param ... Additional arguments passed to \code{spatstat::rose()}.
#' @return A window (class "owin") containing the plotted region.
#' @importFrom spatstat.core rose
#' @export
#' @examples
#' x <- runif(100, 60, 210)
#' rose(x, col="grey", axial=TRUE)
rose <- function(x, bins=NULL, axial=FALSE, clockwise=TRUE, start="N", unit="degree", sub,...){
  if(class(x) == "density"){
    bw <- x$bw
    spatstat.core::rose(
      x%%360, clockwise=clockwise, start = start, unit = unit, ...)
    if(missing(sub)) sub=paste0("Bin width: ", bw)
    title(sub=sub)
  } else {
    if (is.null(bins)==TRUE) {
      bins <- length(x)
    } else {
      bins <- round(bins)
      if (bins <= 0) {
        stop("bins must be non-negative")
      }
    }

    bw <- 360/bins # bin width
    x <- base::as.vector(x%%360)
    breaks=seq(0, 360, bw)

    if(axial==TRUE){
      x2 <- (x+180) %% 360
      x <- hist(x=c(x, x2), plot=FALSE, breaks = breaks)
    }
    spatstat.core::rose(
      x,
      breaks=breaks,
      clockwise=clockwise, start = start, unit = unit, ...)
    if(missing(sub)) sub=paste0("Bin width: ", bw)
    title(sub=sub)
  }
}

#' Inverse Distance Weighting with Directional Data
#'
#' Function for inverse distance weighted interpolation with directional data.
#' Useful for when you are working with data whose unit of measurement is
#' degrees (i.e. the average of 35 degrees and 355 degrees should be 15
#' degrees). It works by finding the shortest distance between two degree marks
#' on a circle.
#'
#' @param values azimuth in degree
#' @param coords An spatial object of locations where the values were measured.
#' @param grid An spatial objector with the locations to predict.
#' @param idp inverse distance weighting power
#' @return \code{data.frame} with the grid coordinates, the interpolated
#' mean values, and the standard deviation for each of the grid
#' @importFrom sp spDists
#' @export
#' @examples
#' values <- c(55, 355)
#' coords <- data.frame(lon = c(1, 2), lat = c(1, 2))
#' sp::coordinates(coords) <- ~lon+lat
#' sp::proj4string(coords) <- "+proj=longlat +datum=WGS84"
#'
#' grid <- data.frame(lon = c(1, 2, 1, 2), lat = c(1, 2, 2, 1))
#' sp::coordinates(grid) <- ~lon+lat
#' sp::proj4string(grid) <- "+proj=longlat +datum=WGS84"
#'
#' ## Perform the inverse distance weighted interpolation
#' res <- idw_circ(values, coords, grid)
#' head(res)
idw_circ <- function (values, coords, grid, idp = 2) {
  stopifnot(length(values) == nrow(coords))
  stopifnot(is.numeric(idp))

  distance <- t(sp::spDists(coords, grid, longlat = TRUE))
  w <- 1/(distance^idp)

  for (i in seq_len(nrow(w))) {
    if (sum(is.infinite(w[i, ])) > 0) {
      w[i, !is.infinite(w[i, ])] <- 0
      w[i, is.infinite(w[i, ])] <- 1
    }
  }
  Z <- apply(w, 1, sum, na.rm = TRUE)

  x <- values %% 180

  sin2 <- (tectonicr:::sind(2 * x))
  cos2 <- (tectonicr:::cosd(2 * x))
  sin2.w <- w %*% diag(sin2)
  cos2.w <- w %*% diag(cos2)

  meansin2 <- apply(sin2.w/Z, 1, sum, na.rm = TRUE)
  meancos2 <- apply(cos2.w/Z, 1, sum, na.rm = TRUE)
  shmax.mean <- (tectonicr:::atan2d(meansin2, meancos2)/2) %% 360

  meanR <- sqrt(meancos2^2 + meansin2^2)
  sd.rad <- sqrt(-2 * log(meanR))/2
  sd <- tectonicr::rad2deg(sd.rad)

  data.frame(sp::coordinates(grid), shmax.mean, sd)
}


#' @title Inverse Distance Weighting with Directional Data
#' @description Function for inverse distance weighted interpolation with directional data. Useful for when you are working with data whose unit of measurement is degrees (i.e. the average of 35 degrees and 355 degrees should be 15 degrees). It works by finding the shortest distance between two degree marks on a circle.
#' @param values A table of points to be interpolated. Table must contain x and y locations, and a column of values to be interpolated.
#' @param locations A table with lat and lon coordinates of the samples.
#' @param x.range Latitude range of the grid
#' @param y.range Longitude range of the grid
#' @param spacing Spacing of the grid. default is 1
#' @param idp The power to use in weight calculation. default is 2
#' @export
#' @importFrom sp coordinates proj4string gridded
interpolate_azimuth <- function(values, locations, x.range, y.range, spacing=1, idp=2){
  # set up the grid
  grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = spacing), y = seq(from = y.range[1], to = y.range[2], by = spacing)) # expand points to grid
  sp::coordinates(grd) <- ~ x + y
  sp::proj4string(grd) <- wgs84
  sp::gridded(grd) <- TRUE

  locations <- as.data.frame(locations)
  sp::coordinates(locations) <- ~X+Y
  sp::proj4string(locations) <- wgs84

  idw_circ(
    values=values,
    coords=locations,
    grid=grd,
    idp=idp
  )
}



#' @title Inverse Distance Weighting interpolation
#' @description This function interpolates a list of samples with location and a value to a table of coordinates, that generally represent a spatial grid. The interpolation is based on inverse distance weighting algoritm with three different methods available for weight calculation.
#' @param values A table of points to be interpolated. Table must contain x and y locations, and a column of values to be interpolated.
#' @param locations A table with lat and lon coordinates of the samples.
#' @param x.range Latitude range of the grid
#' @param y.range Longitude range of the grid
#' @param spacing Spacing of the grid
#' @param idp The power to use in weight calculation.
#' @param ... Other arguments to be passed to \code{\link[gstat]{idw}}
#' @export
#' @importFrom gstat idw
#' @importFrom sp coordinates proj4string CRS gridded
interpolate_idw <-
  function(values,
           locations,
           x.range,
           y.range,
           spacing,
           idp,
           ...) {
    grd <-
      expand.grid(
        x = seq(from = x.range[1], to = x.range[2], by = spacing),
        y = seq(from = y.range[1], to = y.range[2], by = spacing)
      ) # expand points to grid
    sp::coordinates(grd) <- ~ x + y
    sp::proj4string(grd) <-
      sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
    sp::gridded(grd) <- TRUE

    locations <- as.data.frame(locations)
    sp::coordinates(locations) <- ~ X + Y
    sp::proj4string(locations) <-
      sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

    idw <-
      gstat::idw(
        formula = values ~ 1,
        newdata = grd,
        locations =  locations,
        idp = idp
      ) %>%
      as.data.frame()

    names(idw)[1:3] <- c("lon", "lat", "var.pred")
    return(idw$var.pred)
  }

interpolate_krige <-
  function(values,
           locations,
           x.range,
           y.range,
           spacing,
           psill,
           nugget,
           range,
           model,
           ...) {
    grd <-
      expand.grid(
        x = seq(from = x.range[1], to = x.range[2], by = spacing),
        y = seq(from = y.range[1], to = y.range[2], by = spacing)
      ) # expand points to grid
    sp::coordinates(grd) <- ~ x + y
    sp::proj4string(grd) <-
      sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
    sp::gridded(grd) <- TRUE

    locations <- as.data.frame(locations)
    sp::coordinates(locations) <- ~ X + Y
    sp::proj4string(locations) <-
      sp::CRS("+proj=longlat +datum=WGS84 +no_defs")


    if(missing(psill) | missing(nugget) | missing(range) | missing(model)){
      semivariog <-
        gstat::variogram(values ~ 1, locations = locations, data = locations)
      print(plot(semivariog))

      print(semivariog)

      OK <- FALSE
      while (!OK) {
        psill <- as.numeric(readline("Enter psill value (min gamma):"))
        nugget <- as.numeric(readline("Enter nugget value (max gamma):"))
        range <- as.numeric(readline("Enter range value (max dist):"))
        model.i <-
          as.integer(readline("Enter model:\n  1 - Exp\n  2 - Sph\n  3 - Gau\n  4 - Mat"))
        if (model.i == 1) {
          model <- "Exp"
        } else if (model.i == 2) {
          model <- "Sph"
        } else if (model.i == 3) {
          model <- "Gau"
        } else if (model.i == 4) {
          model <- "Mat"
        } else {
          model <- "Exp"
          message("Exp")
        }
        model.variog <-
          gstat::vgm(
            psill = psill,
            model = model,
            nugget = nugget,
            range = range
          )
        fit.variog <- gstat::fit.variogram(semivariog, model.variog)
        print(plot(semivariog, fit.variog))

        decision <- readline("Satisfied? (Y or N):")
        if (decision == "Y" | decision == "y") {
          OK <- TRUE
        } else if (decision == "N" | decision == "n") {
          OK <- FALSE
        }
      }
    }
    model.variog <-
      gstat::vgm(
        psill = psill,
        model = model,
        nugget = nugget,
        range = range
      )

    krige.res <- gstat::krige(
      formula = values ~ 1,
      newdata = grd,
      locations =  locations,
      model = model.variog
    ) %>%
      as.data.frame()

    names(krige.res)[1:3] <- c("lon", "lat", "var.pred")
    return(krige.res$var.pred)
  }





#' Stress2Grid
#'
#' Stress pattern and wavelength analysis
#'
#' @param azi SHmax in degree
#' @param unc Uncertainties of SHmax in degree
#' @param type  methods used for the determination of the orientation of SHmax
#' @param lon,lat
#' @param lon.range,lat.range
#' @param gridsize Specify the spacing of the regular grid in decimal degree.
#' @param min_data minimum number of data per bin. Default is 3
#' @param threshold threshold for deviation of orientation. Default is 25
#' @param arte_thres Maximum distance (in km) of the gridpoint to the next datapoint. Default is 200
#' @param apply_mw Logical. If a method weighting should be applied: Default is FALSE
#' @param apply_qw Logical. If a quality weighting should be applied: Default is TRUE
#' @param dist_weight Distance weighting method which should be used. Specify ‘none’, ‘linear’, or ‘inverse’
#' @param dist_threshold distance weight to prevent overweight of data nearby (0 to 1). Default is 0.1
#' @param R_range
#' @export
stress2grid <- function( azi,
                         unc,
                         type,
                         lon, lat,
                         lon.range, lat.range,
                         gridsize = 2.5,
                         min_data = 3,
                         threshold = 25,
                         arte_thres = 200,
                         apply_mw = FALSE,
                         apply_qw = TRUE,
                         dist_weight = c("linear", "inverse", "none"),
                         dist_threshold = 0.1,
                         R_range = seq(50, 1000, 50)
) {

  # Several fixed variables are assigned and data is loaded:
  deg2rad = pi / 180
  num_r = length(R_range)

  locations <- cbind("X" = lon, "Y" = lat)
  sp::coordinates(locations) <- ~ X + Y
  sp::proj4string(locations) <- "+proj=longlat +ellps=WGS84"

  # WSM method weighting (from 0 to 5)
  if(apply_mw){
    methd <- c(FMS = 4, FMF = 5, BO = 5, DIF = 5,  HF = 4, GF = 5, GV = 4, OC = 2, rest = 1)

    method.weights <- data.frame(
      type = c("FMS", "FMF", "BO", "DIF", "HF", "GF", "GV", "OC", NA),
      w.method = c(4, 5, 5, 5,  4, 5, 4, 2, 1)
    )

  } else {
    w.method <- rep(1, length(lat))
  }

  datas <- data.frame(
    lon = lon,
    lat = lat,
    sin2 = sin(2 * azi * deg2rad),
    cos2 = cos(2 * azi * deg2rad),
    w.method = w.method,
    w.qual = 1 / unc,
    pb_dist = rep(1, length(lat))
  )

  # Regular grid
  if(missing(lon.range) | missing(lat.range)) {
    lon.range <- range(lon, na.rm = TRUE)
    lat.range <- range(lat, na.rm = TRUE)
  }

  G <-
    expand.grid(X = seq(lon.range[1], lon.range[2], gridsize),
                Y = seq(lat.range[1], lat.range[2], gridsize))
  XG  <- G$X
  YG  <- G$Y

  sp::coordinates(G) <- ~ X + Y
  sp::proj4string(G) <- "+proj=longlat +ellps=WGS84"
  sp::gridded(G) <- TRUE
  n_G <- length(XG)

  SH <- matrix()
  gmt_coord <- matrix()
  gmt_azi <- matrix()

  for (i in 1:n_G) {
    distij <-
      t(spDists(locations, matrix(XG[i], YG[i]), longlat = TRUE))
    #distij = ddistance(YG[i],XG[i],datas[, 2],datas[, 1])

    if (arte_thres > 0 && min(distij) > arte_thres) {

    }

    for (k in 1:length(R_range)) {
      R_search <-  R_range[k]
      ids_R <- which(distij < R_search)

      N_in_R  <-  length(ids_R)
      #sq_R  <-  sum(datas$pb_dist[ids_R])

      if (N_in_R < min_data) {
        # not enough data within search radius
        #continue;
      } else if (N_in_R == 1) {
        s0 <- 0
        meanSH <-  datas$w.qual[ids_R]
        mdr <-  distij[ids_R] / R_search
      } else {
        mdr  <-  mean(distij[ids_R]) / R_search

        dist_threshold_scal <-  R_search * dist_threshold

        if (dist_weight == 'inverse') {
          wd_R  <-  1 / max(dist_threshold_scal, distij[ids_R])
        } else if (dist_weight == 'linear') {
          wd_R <-  R_search + 1 - max(dist_threshold_scal, distij[ids_R])
        } else {
          wd_R <- rep(1, length(ids_R))
        }

        sw_R <- sum(wd_R * datas$w.qual[ids_R] * datas$w.method[ids_R])
        array_sin2  <-
          wd_R * datas$pb_dist[ids_R] * datas$w.method[ids_R] * datas$sin2[ids_R]
        array_cos2 <-
          wd_R * datas$pb_dist[ids_R] * datas$w.method[ids_R] * datas$cos2[ids_R]

        # mean value
        sumsin2 <- sum(array_sin2)
        sumcos2 <- sum(array_cos2)
        meansin2 <- sumsin2 / sw_R
        meancos2 <- sumcos2 / sw_R
        meanR <- sqrt(meansin2 ^ 2 + meancos2 ^ 2)

        s0 <- sqrt(-2 * log(meanR)) / (2 * deg2rad)
        meanSH_rad <- atan2(meansin2, meancos2) / 2
        meanSH <- (meanSH_rad / deg2rad) %% 180
      }

      ml <- ((k-1) * 6)
      SH[i, ml+1] <- meanSH
      SH[i, ml+2] <- s0
      SH[i, ml+3] <- R_search
      SH[i, ml+4] <- mdr
      SH[i, ml+5] <- N_in_R

      if(s0 < threshold) {
        gmt_coord[i, ] <- c(XG[i], YG[i], s0, R_search)
        gmt_azi[i, ] <- c(XG[i], YG[i], meanSH)
      } else if(R_search == min(R_range)) {
        gmt_coord[i, ] <- c(XG[i], YG[i], s0, R_search)
      }
    }
  }
}

