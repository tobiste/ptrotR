library(dplyr)
library(sf)

## code to prepare `pb2002_boundaries` dataset goes here
pb2002_boundaries <- sf::st_read("data-raw/PB2002_boundaries.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_wrap_dateline() %>%
  # sf::st_make_valid() %>%
  select(-LAYER)
# plot(pb2002_boundaries)
usethis::use_data(pb2002_boundaries, overwrite = TRUE)

## code to prepare `pb2002_plates` dataset goes here
pb2002_plates <- sf::st_read("data-raw/PB2002_plates.shp") %>%
  sf::st_set_crs("EPSG:4326") %>%
  sf::st_wrap_dateline() %>%
  # sf::st_make_valid() %>%
  select(-LAYER)
# plot(pb2002_plates)
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


greiner <- readxl::read_xlsx("E:/UofC/plate_motion/Greiner.xlsx")
usethis::use_data(greiner, overwrite = TRUE)

gplates_ids <- readxl::read_xlsx("data-raw/Seton_etal_ESR2012_PlateIDs.xlsx") %>%
  mutate(Plate_ID = as.integer(Plate_ID))
usethis::use_data(gplates_ids, overwrite = TRUE)

gp_cmt <- readr::read_delim("data-raw/Muller2019-Young2019-Cao2020_CombinedRotations.rot", delim = "!", col_names = FALSE) %>%
  pull(X2)

gplates_rot <- readr::read_table("data-raw/Muller2019-Young2019-Cao2020_CombinedRotations.rot", comment = "!", col_names = c("plate.rot", "age", "lat", "lon", "angle", "plate.fix")) %>%
  mutate(plate.rot = as.integer(plate.rot), plate.fix = as.numeric(plate.fix), comment = gp_cmt) %>%
  as_tibble() %>%
  left_join(gplates_ids, by = c("plate.rot" = "Plate_ID"), suffix = c(".rot", ".fix")) %>%
  left_join(gplates_ids, by = c("plate.fix" = "Plate_ID"), suffix = c(".rot", ".fix")) %>%
  tidyr::unite("rot", Abbreviation.rot, Abbreviation.fix, sep = "-", remove = F) %>%
  select(plate.rot, age, lat, lon, angle, plate.fix, rot, Abbreviation.rot, Name_and_Description.rot, Abbreviation.fix, Name_and_Description.fix, comment)
# tidyr::separate(comment, into = c("comment", "reference"), sep = "@REF") %>%
# tidyr::separate(reference, into = c("reference", "doi"), sep = "@DOI")
usethis::use_data(gplates_rot, overwrite = TRUE)


st_cmt <- readr::read_delim("data-raw/Seton_etal_ESR2012_2012.1.rot", delim = "!", col_names = FALSE) %>%
  pull(X2)

seton_rot <- readr::read_table("data-raw/Seton_etal_ESR2012_2012.1.rot", comment = "!", col_names = c("plate.rot", "age", "lat", "lon", "angle", "plate.fix")) %>%
  mutate(plate.rot = as.integer(plate.rot), plate.fix = as.numeric(plate.fix), comment = st_cmt) %>%
  left_join(gplates_ids, by = c("plate.rot" = "Plate_ID"), suffix = c(".rot", ".fix")) %>%
  left_join(gplates_ids, by = c("plate.fix" = "Plate_ID"), suffix = c(".rot", ".fix")) %>%
  tidyr::unite("rot", Abbreviation.rot, Abbreviation.fix, sep = "-", remove = F) %>%
  select(plate.rot, age, lat, lon, angle, plate.fix, rot, Abbreviation.rot, Name_and_Description.rot, Abbreviation.fix, Name_and_Description.fix, comment) %>%
  filter(plate.rot != 999)
usethis::use_data(seton_rot, overwrite = TRUE)
