library(dplyr)
library(sf)

## code to prepare `pb2002_boundaries` dataset goes here
pb2002_boundaries <- sf::st_read("data-raw/PB2002_boundaries.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_wrap_dateline() %>%
  #sf::st_make_valid() %>%
  select(-LAYER)
plot(pb2002_boundaries)
usethis::use_data(pb2002_boundaries, overwrite = TRUE)

## code to prepare `pb2002_plates` dataset goes here
pb2002_plates <- sf::st_read("data-raw/PB2002_plates.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_wrap_dateline() %>%
  #sf::st_make_valid() %>%
  select(-LAYER)
plot(pb2002_plates)
usethis::use_data(pb2002_plates, overwrite = TRUE)
