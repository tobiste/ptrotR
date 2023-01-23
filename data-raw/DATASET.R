library(dplyr)
library(sf)

## code to prepare `pb2002_boundaries` dataset goes here
pb2002_boundaries <- sf::st_read("data-raw/PB2002_boundaries.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_wrap_dateline() %>%
  #sf::st_make_valid() %>%
  select(-LAYER)
#plot(pb2002_boundaries)
usethis::use_data(pb2002_boundaries, overwrite = TRUE)

## code to prepare `pb2002_plates` dataset goes here
pb2002_plates <- sf::st_read("data-raw/PB2002_plates.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_wrap_dateline() %>%
  #sf::st_make_valid() %>%
  select(-LAYER)
#plot(pb2002_plates)
usethis::use_data(pb2002_plates, overwrite = TRUE)

kroner_plates2 <- sf::st_read("data-raw/pangea.shp")
plot(kroner_plates2)
usethis::use_data(kroner_plates2, overwrite = TRUE)

kroner_plates <- sf::st_read("data-raw/plates.shp")
plot(kroner_plates)
usethis::use_data(kroner_plates, overwrite = TRUE)

kroner_structures <- sf::st_read("data-raw/structures.shp")
plot(kroner_structures)
usethis::use_data(kroner_structures, overwrite = TRUE)

kroner_coastlines <- sf::st_read("data-raw/coastlines.shp")
plot(kroner_coastlines)
usethis::use_data(kroner_coastlines, overwrite = TRUE)


Torsvik_APWP <- sf::st_read("data-raw/Torsvik2012_apwp_grid_Africa.shp")
plot(Torsvik_APWP)
usethis::use_data(Torsvik_APWP, overwrite = TRUE)

