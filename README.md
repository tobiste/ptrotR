<!-- badges: start -->
[![R-CMD-check](https://github.com/tobiste/PlateTectonicMotionR/workflows/R-CMD-check/badge.svg)](https://github.com/tobiste/PlateTectonicMotionR/actions)
<!-- badges: end -->

# ptrotR: **P**late **T**ectonic **Rot**ations in **R**
**ptrotR** is a free and open-source **R** package for analyzing and reconstructing plate motion. It provides:
- import of GPLATES *.rot files
- interpolate total reconstruction poles between time steps
- rotation of plate polygon (for quick visualization)
- extraction of relative plate motion parameters from sequence of total reconstruction poles

## Prerequisites

You must have R installed on your system (see http://r-project.org). Additionally, to install ptrotR from Github, you also need the devtools package. This can be installed by typing the following code at the R command line prompt:

```
install.packages("remotes")
```

## Installation

The most recent development version  of ptrotR is available from Github and can be installed on your system as follows:

```
remotes::install_github('tobiste/euler') % requiered package from my repository
remotes::install_github('tobiste/ptrotR')
library('ptrotR')
```

## Example

Import total reconstructions from a GPLATES *.rot file:

```
fname <- system.file("Pangea.rot", package="ptrotR")
pangea <- read.gplates(fname)
print(pangea)
```


Extract stage rotations for a specific plate:
```
stage_pols <- extract_stage_rotations(pangea, plate = 137)
```

Create a grid of plate motion vectors:
```
vectors <- plate_motion_grid(stage_pols[1,])

% plot the motion vectors:
tectonicr::axes(vectors$lon, vectors$lat, vectors$azimuth, add = FALSE)
```


## Author
Tobias Stephan

## Useful references
- <div class="csl-entry">Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions. <i>Computers and Geosciences</i>, <i>25</i>(3), 209–216. https://doi.org/10.1016/S0098-3004(98)00160-5</div>

- <div class="csl-entry">Schaeben, H., Kroner, U., &#38; Stephan, T. (2021). Euler Poles of Tectonic Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley, &#38; F. Agterberg (Eds.), <i>Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series</i> (pp. 1–7). Springer Nature Switzerland AG 2021. https://doi.org/10.1007/978-3-030-26050-7_435-1</div>

- Schaeben, H., Kroner, U., & Stephan, T. (2024). Mathematical fundamentals of spherical kinematics of plate tectonics in terms of quaternions. <i>Mathematical Methods in the Applied Sciences</i>, 47(6), 4469–4496. https://doi.org/10.1002/mma.9823
