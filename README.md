# Solar/Satellite Geolocation for Animal Tracking
<!-- badges: start -->
[![R-CMD-check](https://github.com/SWotherspoon/SGAT/workflows/R-CMD-check/badge.svg)](https://github.com/SWotherspoon/SGAT/actions)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/SWotherspoon/SGAT?branch=master&svg=true)](https://ci.appveyor.com/project/SWotherspoon/SGAT)
<!-- badges: end -->

SGAT (pronounced "tags backwards") provides facilities for estimating
broadscale animal motions from archival or satellite tag data.


## Installing

The package is easily installed from GitHub, using the devtools package. 

```R
devtools::install_github("SWotherspoon/SGAT")
```

If you don't have `devtools` installed already, install it first. 

```R
install.packages("devtools")
```

(SGAT otherwise does not need devtools for normal use.)


## TODO

- **Alternate Coordinates**.  The Metropolis sampler may be more
  efficient if locations are represented in geocentric ecliptic
  coordinates. But it is unclear whether the gains in efficiency due
  to higher acceptance rates would outweigh the costs of the
  coordinate transformation and the added code complexity. This is a
  low priority.

- **Parallelisation**.  At this point, the Metropolis samplers are
  only capable of utilizing a single core on a multicore machine.  It
  would be relatively simple to introduce coarse grain parallelism by
  having the samplers draw multiple chains in parallel, using
  something like the multicore facility in the parallel package.
  Unfortunately, at the time of writing there does not seem to be a
  good parallization solution that works equally well on all
  platforms.  This is a low priority.
