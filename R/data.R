#' @title Example rotation file
#' @description Sequence of total reconstruction rotations in GPLATES' .rot format
#'
#' @docType data
#'
#' @usage data(pangea)
#'
#' @format An object of class \code{"finite"}
#' @references <div class="csl-entry">Kroner, U., Roscher, M., &#38; Romer, R. L. (2016). Ancient plate kinematics derived from the deformation pattern of continental crust: Paleo- and Neo-Tethys opening coeval with prolonged Gondwana–Laurussia convergence. <i>Tectonophysics</i>, <i>681</i>, 220–233. https://doi.org/10.1016/j.tecto.2016.03.034</div>
#' @examples
#' data(pangea)
#' head(pangea)
"pangea"

#' @title Example rotation file
#' @description Sequence of total reconstruction rotations in GPLATES' .rot format
#'
#' @docType data
#'
#' @usage data(pannotia)
#'
#' @format An object of class \code{"finite"}
#' @references <div class="csl-entry">Kroner, U., Stephan, T., Romer, R. L., &#38; Roscher, M. (2020). Paleozoic plate kinematics during the Pannotia–Pangaea supercontinent cycle. <i>Geological Society, London, Special Publications</i>, <i>503</i>, SP503-2020–15. https://doi.org/10.1144/SP503-2020-15</div>
#' @examples
#' data(pannotia)
#' head(pannotia)
"pannotia"


#' Plate Boundaries of the Earth
#'
#' Global set of present plate boundaries on the Earth (PB2002) by Bird 2003
#'
#' @docType data
#'
#' @usage data('pb2002_boundaries')
#'
#' @format An object of class \code{"sf"}. LINESTRING
#'
#' @references  Bird, P. (2003), An updated digital model of plate boundaries,
#' *Geochem. Geophys. Geosyst.*, 4, 1027, doi: 10.1029/2001GC000252, 3.
#' @keywords datasets
#' @examples
#' data("pb2002_boundaries")
#' head("pb2002_boundaries")
"pb2002_boundaries"

#' Plates of the Earth
#'
#' Global set of present plates on the Earth (PB2002) by Bird 2003
#'
#' @docType data
#'
#' @usage data('pb2002_plates')
#'
#' @format An object of class \code{"sf"}. MULTIPOLYGON
#'
#' @references  Bird, P. (2003), An updated digital model of plate boundaries,
#' *Geochem. Geophys. Geosyst.*, 4, 1027, doi: 10.1029/2001GC000252, 3.
#' @keywords datasets
#' @examples
#' data("pb2002_plates")
#' head("pb2002_plates")
"pb2002_plates"

#' Plates of the Earth for paleogeographical reconstructions for the
#' Pannotia-Pangea cycle
#'
#' Global set of plates for Pangea-Pannotia supercontinent cycle
#'
#' @docType data
#'
#' @usage data('kroner_plates')
#'
#' @format An object of class \code{"sf"}. MULTIPOLYGON
#'
#' @references <div class="csl-entry">Kroner, U., Stephan, T., Romer, R. L., &#38; Roscher, M. (2020). Paleozoic plate kinematics during the Pannotia–Pangaea supercontinent cycle. <i>Geological Society, London, Special Publications</i>, <i>503</i>, SP503-2020–15. https://doi.org/10.1144/SP503-2020-15</div>
#' @keywords datasets
#' @examples
#' data("kroner_plates")
"kroner_plates"

#' Structures used for paleogeographical reconstructions of the
#' Pannotia-Pangea cycle
#'
#' Global set of structures for Pangea-Pannotia supercontinent cycle
#'
#' @docType data
#'
#' @usage data('kroner_structures')
#'
#' @format An object of class \code{"sf"}. MULTIPOLYGON
#'
#' @references <div class="csl-entry">Kroner, U., Stephan, T., Romer, R. L., &#38; Roscher, M. (2020). Paleozoic plate kinematics during the Pannotia–Pangaea supercontinent cycle. <i>Geological Society, London, Special Publications</i>, <i>503</i>, SP503-2020–15. https://doi.org/10.1144/SP503-2020-15</div>
#' @keywords datasets
#' @examples
#' data("kroner_structures")
"kroner_structures"

#' Coastlines used for paleogeographical reconstructions of the
#' Pannotia-Pangea cycle
#'
#' Global set of coastlines for Pangea-Pannotia supercontinent cycle
#'
#' @docType data
#'
#' @usage data('kroner_coastlines')
#'
#' @format An object of class \code{"sf"}. MULTIPOLYGON
#'
#' @references <div class="csl-entry">Kroner, U., Stephan, T., Romer, R. L., &#38; Roscher, M. (2020). Paleozoic plate kinematics during the Pannotia–Pangaea supercontinent cycle. <i>Geological Society, London, Special Publications</i>, <i>503</i>, SP503-2020–15. https://doi.org/10.1144/SP503-2020-15</div>
#' @keywords datasets
#' @examples
#' data("kroner_coastlines")
"kroner_coastlines"

#' Torsvik et al. 2012 APWP
#'
#' Torsvik et al. 2012 APWP in African coordinates
#'
#' @docType data
#'
#' @usage data('Torsvik_APWP')
#'
#' @format An object of class \code{"sf"}. MULTIPOLYGON
#'
#' @keywords datasets
#' @examples
#' data("Torsvik_APWP")
"Torsvik_APWP"

#' The rotations opening the Central and Northern Atlantic Ocean
#'
#' Total reconstruction rotations of Greiner and Neugebauer (2013)
#'
#' @docType data
#'
#' @usage data('greiner')
#'
#' @format An object of class \code{"tibble"}
#'
#' @keywords datasets
#' @examples
#' data("greiner")
"greiner"

#' ID list of plates used in GPLATES
#'
#' @docType data
#'
#' @usage data('gplates_ids')
#'
#' @format An object of class \code{"tibble"}
#'
#' @keywords datasets
#' @examples
#' data("gplates_ids")
"gplates_ids"

#' Rotation parameters of GPlates
#'
#' @docType data
#'
#' @usage data('gplates_rot')
#'
#' @format An object of class \code{"tibble"}
#'
#' @keywords datasets
#' @examples
#' data("gplates_rot")
"gplates_rot"

#' Rotation parameters from Seton et al. 2012
#'
#' @docType data
#'
#' @usage data('seton_rot')
#'
#' @format An object of class \code{"tibble"}
#'
#' @keywords datasets
#' @examples
#' data("seton_rot")
"seton_rot"
