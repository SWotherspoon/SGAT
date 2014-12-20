#' Solar/Satellite Geolocation for Animal Tracking
#'
#' Provides facilities for estimating broadscale animal motions from
#' archival or satellite tag data, similar to the \pkg{GeoLight} and
#' \pkg{tripEstimation} packages.
#'
#' @name SGAT-package
#' @docType package
#' @author S. Wotherspoon, M. Sumner and S. Lisovski.
NULL



## Solar Zenith/Sunrise/Sunset calculations
##
## The functions presented here are based on code and the excel
## spreadsheet from the NOAA site
##
##       http://www.esrl.noaa.gov/gmd/grad/solcalc/
##


##' Calculate solar time, the equation of time and solar declination
##'
##' The solar time, the equation of time and the sine and cosine of
##' the solar declination are calculated for the times specified by
##' \code{tm} using the same methods as
##' \url{www.esrl.noaa.gov/gmd/grad/solcalc/}.
##' @title Solar Time and Declination
##' @param tm a vector of POSIXct times.
##' @return A list containing the following vectors.
##' \item{\code{solarTime}}{the solar time (degrees)}
##' \item{\code{eqnTime}}{the equation of time (minutes of time)}
##' \item{\code{sinSolarDec}}{sine of the solar declination}
##' \item{\code{cosSolarDec}}{cosine of the solar declination}
##' @seealso \code{\link{zenith}}
##' @examples
##' ## Current solar time
##' solar(Sys.time())
##' @export
solar <- function(tm) {

  rad <- pi/180

  ## Time as Julian day (R form)
  Jd <- as.numeric(tm)/86400.0+2440587.5

  ## Time as Julian century [G]
  Jc <- (Jd-2451545)/36525

  ## The geometric mean sun longitude (degrees) [I]
  L0 <- (280.46646+Jc*(36000.76983+0.0003032*Jc))%%360


  ## Geometric mean anomaly for the sun (degrees) [J]
  M <- 357.52911+Jc*(35999.05029-0.0001537*Jc)

  ## The eccentricity of earth's orbit [K]
  e <- 0.016708634-Jc*(0.000042037+0.0000001267*Jc)

  ## Equation of centre for the sun (degrees) [L]
  eqctr <- sin(rad*M)*(1.914602-Jc*(0.004817+0.000014*Jc))+
    sin(rad*2*M)*(0.019993-0.000101*Jc)+
      sin(rad*3*M)*0.000289

  ## The true longitude of the sun (degrees) [M]
  lambda0 <- L0 + eqctr

  ## The apparent longitude of the sun (degrees) [P]
  omega <- 125.04-1934.136*Jc
  lambda <- lambda0-0.00569-0.00478*sin(rad*omega)


  ## The mean obliquity of the ecliptic (degrees) [Q]
  seconds <- 21.448-Jc*(46.815+Jc*(0.00059-Jc*(0.001813)))
  obliq0 <- 23+(26+(seconds/60))/60

  ## The corrected obliquity of the ecliptic (degrees) [R]
  omega <- 125.04-1934.136*Jc
  obliq <- obliq0 + 0.00256*cos(rad*omega)

  ## The equation of time (minutes of time) [U,V]
  y <- tan(rad*obliq/2)^2
  eqnTime <- 4/rad*(y*sin(rad*2*L0) -
                      2*e*sin(rad*M) +
                      4*e*y*sin(rad*M)*cos(rad*2*L0) -
                      0.5*y^2*sin(rad*4*L0) -
                      1.25*e^2*sin(rad*2*M))

  ## The sun's declination (radians) [T]
  solarDec <- asin(sin(rad*obliq)*sin(rad*lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)

  ## Solar time unadjusted for longitude (degrees) [AB!!]
  ## Am missing a mod 360 here, but is only used within cosine.
  solarTime <- ((Jd-0.5)%%1*1440+eqnTime)/4
  #solarTime <- ((Jd-2440587.5)*1440+eqnTime)/4

  ## Return solar constants
  list(solarTime=solarTime,
       eqnTime=eqnTime,
       sinSolarDec=sinSolarDec,
       cosSolarDec=cosSolarDec)
}



##' Calculate the solar zenith angle for given times and locations
##'
##' \code{zenith} uses the solar time and declination calculated by
##' \code{solar} to compute the solar zenith angle for given times and
##' locations, using the same methods as
##' \url{www.esrl.noaa.gov/gmd/grad/solcalc/}.  This function does not
##' adjust for atmospheric refraction see \code{\link{refracted}}.
##' @title Solar Zenith Angle
##' @param sun list of solar time and declination computed by
##' \code{solar}.
##' @param lon vector of longitudes.
##' @param lat vector latitudes.
##' @return A vector of solar zenith angles (degrees) for the given
##' locations and times.
##' @seealso \code{\link{solar}}
##' @examples
##' ## Approx location of Sydney Harbour Bridge
##' lon <- 151.211
##' lat <- -33.852
##' ## Solar zenith angle for noon on the first of May 2000
##' ## at the Sydney Harbour Bridge
##' s <- solar(as.POSIXct("2000-05-01 12:00:00","EST"))
##' zenith(s,lon,lat)
##' @export
zenith <- function(sun,lon,lat) {

  rad <- pi/180

  ## Suns hour angle (degrees) [AC!!]
  hourAngle <- sun$solarTime+lon-180
  #hourAngle <- sun$solarTime%%360+lon-180

  ## Cosine of sun's zenith [AD]
  cosZenith <- (sin(rad*lat)*sun$sinSolarDec+
                cos(rad*lat)*sun$cosSolarDec*cos(rad*hourAngle))

  ## Limit to [-1,1] [!!]
  cosZenith[cosZenith > 1] <- 1
  cosZenith[cosZenith < -1] <- -1

  ## Ignore refraction correction
  acos(cosZenith)/rad
}



##' Adjust the solar zenith angle for atmospheric refraction.
##'
##' Given a vector of solar zeniths computed by \code{\link{zenith}},
##' \code{refracted} calculates the solar zeniths adjusted for the
##' effect of atmospheric refraction.
##'
##' \code{unrefracted} is the inverse of \code{refracted}. Given a
##' (single) solar zenith adjusted for the effect of atmospheric
##' refraction, \code{unrefracted} calculates the solar zenith as
##' computed by \code{\link{zenith}}.
##'
##' @title Atmospheric Refraction
##' @param zenith zenith angle (degrees) to adjust.
##' @return vector of zenith angles (degrees) adjusted for atmospheric
##' refraction.
##' @examples
##' ## Refraction causes the sun to appears higher on the horizon
##' refracted(85:92)
##' ## unrefracted gives unadjusted zenith
##' unrefracted(refracted(90))
##' @export
refracted <- function(zenith) {
  rad <- pi/180
  elev <- 90-zenith
  te <- tan((rad)*elev)
  ## Atmospheric Refraction [AF]
  r <- ifelse(elev>85,0,
              ifelse(elev>5,58.1/te-0.07/te^3+0.000086/te^5,
                     ifelse(elev>-0.575,
                            1735+elev*(-518.2+elev*(103.4+elev*(-12.79+elev*0.711))),-20.772/te)))
  ## Corrected Zenith [90-AG]
  zenith-r/3600
}




##' @rdname refracted
##' @export
unrefracted <- function(zenith)
  uniroot(function(x) refracted(x)-zenith,c(zenith,zenith+2))



##' Estimate time of sunrise or sunset for a given location given the
##' approximate solar time of twilight
##'
##' Solar declination and equation of time vary slowly over the day,
##' and so the values of the Solar declination and equation of time at
##' sunrise/sunset ca be caclulated approximately is an approximate
##' time of sunrise/sunset is known. The sun's hour angle and hence
##' sunrise/sunset for the required zenith can then be calculated from
##' these approximations.
##'
##' Note this function returns the time of twilight in solar time.
##' @title Solar Time of Sunrise and Sunset
##' @param solar output of \code{solar} for approximate times of
##' twilight.
##' @param lon vector of longitudes.
##' @param lat vector of latitudes.
##' @param rise logical vector indicating whether to compute rise or
##' set.
##' @param zenith the solar zenith angle that defines twilight.
##' @return a vector of twilight times in solar time (degrees)
##' @seealso \code{\link{twilight}}
##' @export
twilight.solartime <- function(solar,lon,lat,rise,zenith=96) {
  rad <- pi/180
  cosz <- cos(rad*zenith)
  cosHA <- (cosz-sin(rad*lat)*solar$sinSolarDec)/(cos(rad*lat)*solar$cosSolarDec)
  ## Compute the sun's hour angle from its declination for this location
  hourAngle <- ifelse(rise,360,0)+ifelse(rise,-1,1)*suppressWarnings(acos(cosHA)/rad)
  ## Solar time of sunrise at this zenith angle, lon and lat
  #(hourAngle+180-lon)%%360
  #360*(solar$solarTime%/%360)+solarTime
  solarTime <- (hourAngle+180-lon)%%360
  (solarTime-solar$solarTime+180)%%360-180+solar$solarTime
}


##' Estimate time of sunrsie or sunset for a given day and location
##'
##' \code{twilight} uses an iterative algorithm to estimate times of
##' sunrise and sunset.
##'
##' Note that these functions return the twilight that occurs on the
##' same date GMT as \code{tm}, and so sunset may occur before
##' sunrise, depending upon latitude.
##'
##' Solar declination and equation of time vary slowly over the day,
##' and so the values of the Solar declination and equation of time at
##' sunrise/sunset are well approximated by their values at 6AM/6PM
##' local time. The sun's hour angle and hence sunrise/sunset for the
##' required zenith can then be caclulates from these approximations.
##' The calculation is then repeated using the approximate
##' sunrise/sunset times to derive more accurate values of the Solar
##' declination and equation of time and hence better approximations
##' of sunrise/sunset.  The process is repreated and is accurate to
##' less than 2 seconds within 2 or 3 iterations.
##'
##' \code{sunrise} and \code{sunset} are simple wrappers for
##' \code{twilight}.
##' @title Times of Sunrise and Sunset
##' @param tm vector of approximate times of twilight.
##' @param lon vector of longitudes.
##' @param lat vector of latitudes.
##' @param rise logical vector indicating whether to compute rise or
##' set.
##' @param zenith the solar zenith angle that defines twilight.
##' @param iters number of iteratve refinements made to the initial
##' approximation.
##' @return a vector of twilight times.
##' @examples
##' ## Approx location of Santa Barbara
##' lon <- -119.7022
##' lat <- 34.4191
##' ## Sunrise and sunset for 8th April 2013 at Santa Barbara
##' day <- as.POSIXct("2013-04-08","GMT")
##' sunrise(day,lon,lat)
##' sunset(day,lon,lat)
##' @export
twilight <- function(tm,lon,lat,rise,zenith=96,iters=3) {

  ## Compute date
  date <- as.POSIXlt(tm)
  date$hour <- date$min <- date$sec <- 0
  date <- as.POSIXct(date,"GMT")

  lon <- (lon+180)%%360-180
  ## GMT equivalent of 6am or 6pm local time
  twl <- date+240*(ifelse(rise,90,270)-lon)
  ## Iteratively improve estimate
  for(k in seq_len(iters)) {
    s <- solar(twl)
    s$solarTime <- s$solarTime%%360
    solarTime <- 4*twilight.solartime(s,lon,lat,rise,zenith)-s$eqnTime
    twl <- date+60*solarTime
  }
  twl
}

##' @rdname twilight
##' @export
sunrise <- function(tm,lon,lat,zenith=96,iters=3)
  twilight(tm,lon,lat,rise=TRUE,zenith=zenith,iters=iters)

##' @rdname twilight
##' @export
sunset <- function(tm,lon,lat,zenith=96,iters=3)
  twilight(tm,lon,lat,rise=FALSE,zenith=zenith,iters=iters)


##' Midpoints of a path
##'
##' Compute the midpoints of a sequence of locations along a path.
##' @title Path Midpoints
##' @param p a two column matrix of (lon,lat) locations along the
##' path.
##' @param fold should the longitudes be folded into [-180,180).
##' @return a two column matrix of (lon,lat) midpoints.
##' @export
trackMidpts <- function(p,fold=FALSE) {
  n <- nrow(p)
  rad <- pi/180
  p <- rad*p
  dlon <- diff(p[,1L])
  lon1 <- p[-n,1L]
  lat1 <- p[-n,2L]
  lat2 <- p[-1L,2L]
  bx <- cos(lat2)*cos(dlon)
  by <- cos(lat2)*sin(dlon)
  lat <- atan2(sin(lat1)+sin(lat2),sqrt((cos(lat1)+bx)^2+by^2))/rad
  lon <- (lon1+atan2(by,cos(lat1)+bx))/rad
  if(fold) lon <- (lon+180)%%360-180
  cbind(lon,lat)
}


##' Distances along a path
##'
##' The \code{trackDist} computes the great circle distances (in km)
##' between successive locations along path. The \code{trackDist2}
##' accepts a second sequence of intermediate points, and computes the
##' great circle distances along the dog leg paths from \code{x[i,]}
##' to \code{z[i,]} to \code{x[i+1,]}.
##'
##' @title Distance along a path
##' @param x a two column matrix of (lon,lat) locations along the path.
##' @param z a two column matrix of (lon,lat) intermediate locations
##' along the path.
##' @return vector of interpoint distances (km)
##' @export
trackDist <- function(x) {
  n <- nrow(x)
  rad <- pi/180
  cosx2 <- cos(rad*x[,2L])
  sinx2 <- sin(rad*x[,2L])

  6378.137*acos(pmin.int(cosx2[-n]*cosx2[-1L]*cos(rad*(x[-1L,1L]-x[-n,1L]))+sinx2[-n]*sinx2[-1L],1))
}


##' @rdname trackDist
##'@export
trackDist2 <- function(x,z) {
    n <- nrow(x)
    rad <- pi/180
    cosx2 <- cos(rad*x[,2L])
    sinx2 <- sin(rad*x[,2L])
    cosz2 <- cos(rad*z[,2L])
    sinz2 <- sin(rad*z[,2L])

    6378.137*(acos(pmin.int(cosx2[-n] *cosz2*cos(rad*(z[,1L]-x[-n, 1L]))+sinx2[-n] *sinz2,1))+
              acos(pmin.int(cosx2[-1L]*cosz2*cos(rad*(z[,1L]-x[-1L,1L]))+sinx2[-1L]*sinz2,1)))
  }




##' Bearing changes along a track
##'
##' The \code{trackBearingChange} computes the change in bearing between
##' successive locations along path. The \code{trackBearingChange2}
##' accepts a second sequence of intermediate points, and computes the
##' change in bearing along the dog leg paths from \code{x[i,]}
##' to \code{z[i,]} to \code{x[i+1,]}.
##'
##' @title Distance along a path
##' @param x a two column matrix of (lon,lat) locations along the path.
##' @param z a two column matrix of (lon,lat) intermediate locations
##' along the path.
##' @return vector of changes in bearing (degrees)
##' @export
trackBearingChange <- function(x) {
  n <- nrow(x)
  rad <- pi/180
  cosx2 <- cos(rad*x[,2L])
  sinx2 <- sin(rad*x[,2L])


  ## Bearing from one x to the next
  bs <- atan2(sin(rad*(x[-1L,1L]-x[-n,1L]))*cosx2[-1L],
              cosx2[-n]*sinx2[-1L]-sinx2[-n]*cosx2[-1L]*cos(rad*(x[-1L,1L]-x[-n,1L])))/rad
  ## Difference bs and fold difference into [-180,180)
  (bs[1-n]-bs[-1]+180)%%360-180
}


##' @rdname trackBearingChange
##' @export
trackBearingChange2 <- function(x,z) {
  n <- nrow(x)
  rad <- pi/180
  cosx2 <- cos(rad*x[,2L])
  sinx2 <- sin(rad*x[,2L])
  cosz2 <- cos(rad*z[,2L])
  sinz2 <- sin(rad*z[,2L])
  dLon1 <- rad*(x[-n,1L]-z[,1L])
  dLon2 <- rad*(x[-1,1L]-z[,1L])

  ## Bearing from z to previous x
  b1 <- atan2(sin(dLon1)*cosx2[-n],
              cosz2*sinx2[-n]-sinz2*cosx2[-n]*cos(dLon1))/rad

  ## Bearing from z to next
  b2 <- atan2(sin(dLon2)*cosx2[-1L],
              cosz2*sinx2[-1L]-sinz2*cosx2[-1L]*cos(dLon2))/rad
  ## Reverse b1 and fold difference into [-180,180)
  (b2-b1)%%360-180
}






##' Convert streams of twilights to sunrise/sunset pairs
##'
##' This function converts the twilight, rise format used by Stella
##' and Estelle into successive sunrise and sunset pairs.
##' @title Extract sunrise/sunset pairs
##' @param twilight the observed times of twilight as POSIXct.
##' @param rise logical vector indicating which twilights are sunrise.
##' @return A dataframe with columns
##' \item{\code{Twilight1}}{times of earlier twilight as POSIXct objects}
##' \item{\code{Twilight2}}{times of later twilight as POSIXct objects}
##' \item{\code{Day}}{logical vector indicating whether the twilights span a day.}
##' \item{\code{Mid}}{the midpont of the two twilights.}
##' @export
twilight.pairs <- function(twilight,rise) {
  n <- length(twilight)
  t1 <- twilight[-n]
  r1 <- rise[-n]
  t2 <- twilight[-1L]
  ## Must have one rise and set per day, and must be less than 24
  ## hours apart.
  keep <- (r1!=rise[-1L]) & (as.numeric(t2)-as.numeric(t1) < 86400)
  mid <- .POSIXct(as.numeric(t1)+(as.numeric(t2)-as.numeric(t1))/2,"GMT")
  data.frame(Twilight1=t1[keep],
             Twilight2=t2[keep],
             Day=r1[keep],
             Mid=mid[keep])
}



##' Estimate location from consecutive twilights
##'
##' These functions estimate the location of a stationary observer
##' given the times at which the observer sees two successive
##' twilights. \code{threshold.estimate} estimates locations given
##' pairs of times of sunrise and sunset. \code{threshold.location}
##' is a wrapper for \code{threshold.estimate} that estimates
##' locations given a sequence twilight times and rise indicators,
##' while \code{threshold.path} interpolates the estimates generated
##' by \code{threshold.location} to give locations at a sequence of
##' arbitrary times.
##'
##' Longitude is estimated by computing apparent time of local noon
##' from sunrise and sunset, and determining the longitude for which
##' this is noon.  Latitude is estimated from the required zenith and
##' the sun's hour angle for both sunrise and sunset, and averaged.
##'
##' When the solar declination is near zero (at the equinoxes)
##' latitude estimates are extremely sensitive to errors.  Where the
##' sine of the solar declination is less than \code{tol}, the
##' latitude estimates are returned as \code{NA}.
##'
##' The \code{threshold.path} function interpolates the estimates
##' generated by \code{threshold.location} to produce estimates at the
##' arbitrary set of times specified the by the \code{time} argument.
##' If \code{time} is \code{NULL}, \code{threshold.path} returns the
##' estimates generated by \code{threshold.location}. If
##' \code{unfold=TRUE}, \code{threshold.path} attempts to construct a
##' continuous path that does not wrap longitudes into
##' (-180,180]. However, this process can fail if the observer crosses
##' the dateline near an equinox, and it may be necessary to make
##' manual adjustments.
##'
##' These functions provides the same basic functionality of the
##' \code{coord} function from \pkg{GeoLight}, but are based on
##' different astronomical approximations.
##' @title Simple Threshold Geolocation Estimates
##' @param trise vector of sunrise times.
##' @param tset vector of sunset times.
##' @param twilight the observed times of twilight as POSIXct.
##' @param rise logical vector indicating which twilights are sunrise.
##' @param time times for which locations are required.
##' @param zenith the solar zenith angle that defines twilight.
##' @param tol tolerance on the sine of the solar declination.
##' @param unfold if \code{TRUE}, unfold longitudes across the dateline.
##' @return \code{threshold.estimate} returns estimated locations as a
##' two column (lon,lat) matrix.  \code{threshold.location} and
##' \code{threshold.path} return a list with components
##' \item{\code{time}}{the time as POSIXct.}
##' \item{\code{x}}{a two column matrix of (lon,lat) locations.}
##' @seealso \code{\link{zenith}}
##' @export
threshold.estimate <- function(trise,tset,zenith=96,tol=0) {
  rad <- pi/180
  sr <- solar(trise)
  ss <- solar(tset)
  cosz <- cos(rad*zenith)
  lon <- -(sr$solarTime+ss$solarTime+ifelse(sr$solarTime<ss$solarTime,360,0))/2
  lon <- (lon+180)%%360-180

  ## Compute latitude from sunrise
  hourAngle <- sr$solarTime+lon-180
  a <- sr$sinSolarDec
  b <- sr$cosSolarDec*cos(rad*hourAngle)
  x <- (a*cosz-sign(a)*b*suppressWarnings(sqrt(a^2+b^2-cosz^2)))/(a^2+b^2)
  lat1 <- ifelse(abs(a)>tol,asin(x)/rad,NA)

  ## Compute latitude from sunset
  hourAngle <- ss$solarTime+lon-180
  a <- ss$sinSolarDec
  b <- ss$cosSolarDec*cos(rad*hourAngle)
  x <- (a*cosz-sign(a)*b*suppressWarnings(sqrt(a^2+b^2-cosz^2)))/(a^2+b^2)
  lat2 <- ifelse(abs(a)>tol,asin(x)/rad,NA)

  ## Average latitudes
  cbind(lon=lon,lat=rowMeans(cbind(lat1,lat2),na.rm=TRUE))
}



##' @rdname threshold.estimate
##' @export
threshold.location <- function(twilight,rise,zenith=96,tol=0.08) {
  ## Convert to sunrise/sunset pairs
  pr <- twilight.pairs(twilight,rise)
  ## Estimate locations
  ps <- threshold.estimate(ifelse(pr$Day,pr$Twilight1,pr$Twilight2),
                           ifelse(pr$Day,pr$Twilight2,pr$Twilight1),
                           zenith=zenith,tol=tol)
  list(time=pr$Mid,x=ps)
}


##' @rdname threshold.estimate
##' @export
threshold.path <- function(twilight,rise,time=twilight,zenith=96,tol=0.08,unfold=TRUE) {
  ## Estimate locations
  ls <- threshold.location(twilight,rise,zenith=zenith,tol=tol)
  if(!is.null(time)) {
    ## Interpolate the non-missing longitudes
    keep <- !is.na(ls$x[,1L])
    ts <- ls$time[keep]
    lon <- ls$x[keep,1L]
    if(unfold) lon <- cumsum(c(lon[1L],(diff(lon)+180)%%360-180))
    lon <- approx(x=ts,y=lon,xout=time,rule=2)$y
    ## Interpolate the non-missing latitudes
    keep <- !is.na(ls$x[,2L])
    ts <- ls$time[keep]
    lat <- ls$x[keep,2L]
    lat <- approx(x=ts,y=lat,xout=time,rule=2)$y
    ls <- list(time=time,x=cbind(lon,lat))
  }
  ls
}



##' Estimate locations by the threshold method assuming a
##' nonstationary observer and errors in estimated twilights
##'
##' Given the times of a single sunrise and sunset pair,
##' \code{threshold.sensitivity} estimates the location of the tagged
##' animal at sunrise and at sunset assuming that during this time the
##' animal moves no further than a given maximum range, and that the
##' observed times of sunrise and sunset contain an additive log
##' Normally distributed error with known mean and variance. These
##' errors are directed so that observed sunrise occurs earlier than
##' true sunrise, and the observed sunset occurs later than true
##' sunrise.
##'
##' \code{threshold.sensitivity} implements a Metropolis sampler to
##' draw samples from the posterior distribution for the sunrise and
##' sunset.
##' @title Threshold Geolocation Sensitivity
##' @param rise observed time of sunrise as POSIXct.
##' @param set observed time of sunset as POSIXct.
##' @param zenith the solar zenith angle that defines twilight.
##' @param range maximum range of travel between twilights (km).
##' @param sr.mulog log mean parameter for the Log Normal distribution
##' of sunrise errors.
##' @param sr.sdlog log standard deviation parameter for the Log
##' Normal distribution of sunrise errors.
##' @param ss.mulog log mean parameter for the Log Normal distribution
##' of sunset errors.
##' @param ss.sdlog log standard deviation parameter for the Log
##' Normal distribution of sunset errors.
##' @param sr.proposal function for drawing from the proposal
##' distribution for sunrise location.
##' @param ss.proposal function for drawing from the proposal
##' distribution for sunrise location.
##' @param n.thin rate at which to thin samples.
##' @param n.iters total number of samples to draw.
##' @return a list with three components
##' \item{\code{p0}}{the threshold estimate}
##' \item{\code{rise}}{the sampled sunrise locations as a two column
##' matrix}
##' \item{\code{set}}{the sampled sunset locations as a two column
##' matrix}
##' @export
threshold.sensitivity <- function(rise,set,zenith=96,range=100,
                                  sr.mulog,sr.sdlog,ss.mulog,ss.sdlog,
                                  sr.proposal,ss.proposal,
                                  n.thin=10,n.iters=1000) {

  ## Great circle distance (km)
  gcdist <- function(a,b) {
    rad <- pi/180
    6378.137*acos(pmin.int(
      cos(rad*a[2L])*cos(rad*b[2L])*
      cos(rad*(b[1L]-a[1L]))+sin(rad*a[2L])*
      sin(rad*b[2L]),
      1))
  }

  ## Solar properties at ss/sr
  sr <- solar(rise)
  ss <- solar(set)

  ## Initialize chain from best estimate of location if stationary
  p0 <- as.vector(threshold.estimate(rise,set,zenith))
  sr.p <- p0
  ss.p <- p0

  ## Initialize cached solar times
  sr.solar <- twilight.solartime(sr,sr.p[1L],sr.p[2L],TRUE,zenith)
  ss.solar <- twilight.solartime(ss,ss.p[1L],ss.p[2L],FALSE,zenith)

  ## Ensure approximation consistency
  sr$solarTime <- sr.solar
  ss$solarTime <- ss.solar

  ## Initialize cached log posterior
  sr.logp <- dlnorm(0,sr.mulog,sr.sdlog,log=TRUE)
  ss.logp <- dlnorm(0,ss.mulog,ss.sdlog,log=TRUE)

  P.sr <- matrix(0,2L,n.iters)
  P.ss <- matrix(0,2L,n.iters)
  for(k1 in 1:n.iters) {
    for(k2 in 1:n.thin) {

      ## Propose new sunrise location
      sr.p1 <- sr.proposal(sr.p)
      sr.solar1 <- twilight.solartime(sr,sr.p1[1L],sr.p1[2L],TRUE,zenith)
      if(is.finite(sr.solar1) && gcdist(sr.p1,ss.p) < range) {
        ## When proposal in range compute time error
        sr.delta <- 4*(sr$solarTime-sr.solar1)
        if(sr.delta>0) {
          ## Metropolis rule
          sr.logp1 <- dlnorm(sr.delta,sr.mulog,sr.sdlog,log=TRUE)
          if(sr.logp1-sr.logp > log(runif(1))) {
            ## Accept proposal
            sr.p <- sr.p1
            sr.solar <- sr.solar1
            sr.logp <- sr.logp1
          }
        }
      }

      ## Propose new sunset location
      ss.p1 <- ss.proposal(ss.p)
      ss.solar1 <- twilight.solartime(ss,ss.p1[1L],ss.p1[2L],FALSE,zenith)
      if(is.finite(ss.solar1) && gcdist(sr.p,ss.p1) < range) {
        ## When proposal in range compute time error
        ss.delta <- 4*(ss.solar1-ss$solarTime)
        if(ss.delta>0) {
          ## Metropolis rule
          ss.logp1 <- dlnorm(ss.delta,ss.mulog,ss.sdlog,log=TRUE)
          if(ss.logp1-ss.logp > log(runif(1))) {
            ## Accept proposal
            ss.p <- ss.p1
            ss.solar <- ss.solar1
            ss.logp <- ss.logp1
          }
        }
      }

    }
    ## Record locations at sr/ss
    P.sr[,k1] <- sr.p
    P.ss[,k1] <- ss.p
  }
  list(p0=p0,rise=t(P.sr),set=t(P.ss))
}


##' Simulate zenith angles, times and locations of twilight along a
##' specified track.
##'
##' Given times, longitudes and latitudes that specify a template
##' track, \code{zenith.simulate} interpolates the template onto the
##' new times specified by \code{tm.out} and computes the solar zenith
##' angle at each point along the new track. Given a dataframe
##' generated by \code{zenith.simulate}, \code{twilight.simulate}
##' computes times and locations of sunrise and sunset based on the
##' simulated zenith angles. The \code{twilight.perturb} adds a given
##' vector of errors (in minutes) to the twilights in a dataframe
##' generated by \code{twilight.simulate}, in such a way that a
##' positive error causes sunrise to occur later and sunset earlier.
##' @title Solar Zenith and Twilight Simulation
##' @param tm vector of times that specify the template track.
##' @param lon vector of longitude that specify the template track.
##' @param lat vector of latitude that specify the template track.
##' @param tm.out vector of times to which the template is resampled.
##' @param dfz a dataframe generated with \code{zenith.simulate}.
##' @param zenith the solar zenith angle that defines twilight.
##' @param dft a dataframe generated with \code{twilight.simulate}.
##' @param err a vector of adjustments (in minutes) to the twilight
##' times.
##' @return \code{zenith.simulate} returns a data frame with
##' components
##' \item{\code{Date}}{times along the simulated track}
##' \item{\code{Lon}}{longitudes along the simulated track}
##' \item{\code{Lat}}{latitudes along the simulated track}
##' \item{\code{Zenith}}{zenith angles along the simulated track}
##' \code{twilight.simulate} returns a data frame of twilights with
##' components
##' \item{\code{Twilight}}{times of twilight}
##' \item{\code{Rise}}{is this a sunrise}
##' \item{\code{Lon}}{longitude at twilight}
##' \item{\code{Lat}}{latitude at twilight}
##' @export
zenith.simulate <- function(tm,lon,lat,tm.out) {
  ## unwrap longitudes
  lon <- cumsum(c(lon[1L],(diff(lon)+180)%%360-180))
  ## Interpolate track
  keep <- !is.na(lon)
  lon.out <- approx(tm[keep],lon[keep],tm.out,rule=2)$y
  keep <- !is.na(lat)
  lat.out <- approx(tm[keep],lat[keep],tm.out,rule=2)$y
  ## Compute zenith angles
  z <- zenith(solar(tm.out),lon.out,lat.out)
  data.frame(Date=tm.out,
             Lon=lon.out,
             Lat=lat.out,
             Zenith=z)
}



##' @rdname zenith.simulate
##' @export
twilight.simulate <- function(dfz,zenith=96) {

  n <- nrow(dfz)

  ## Compute indexes for sunrise and sunset
  sr.k <- which(dfz$Zenith[-n] >= zenith & dfz$Zenith[-1L] < zenith)
  ss.k <- which(dfz$Zenith[-n] < 96 & dfz$Zenith[-1L] >= 96)
  ## Interleave sunrise and sunset
  ord <- order(c(sr.k,ss.k))
  k <- c(sr.k,ss.k)[ord]
  rise <- rep(c(T,F),c(length(sr.k),length(ss.k)))[ord]
  ## Interpolation weights
  w <- (zenith-dfz$Zenith[k])/(dfz$Zenith[k+1L]-dfz$Zenith[k])

  ## Interplated times and locations of twilight
  data.frame(Twilight=dfz$Date[k] + w*(as.vector(dfz$Date[k+1L])-as.vector(dfz$Date[k])),
             Rise=rise,
             Lon=dfz$Lon[k] + w*(dfz$Lon[k+1L]-dfz$Lon[k]),
             Lat=dfz$Lat[k] + w*(dfz$Lat[k+1L]-dfz$Lat[k]))
}


##' @rdname zenith.simulate
##' @export
twilight.perturb <- function(dft,err) {
  dft$Twilight <- dft$Twilight + ifelse(dft$Rise,60*err,-60*err)
  dft
}


##' Twilight residuals for known locations
##'
##' Given known locations \code{p}, this function calculates the
##' difference between the predicted and observed twilights for those
##' locations, where the sign is chosen so that obscuration of the
##' sensor leads to a positive residual.
##'
##' @title Twilight Residuals
##' @param twilight the observed times of twilight as POSIXct.
##' @param rise logical vector indicating which twilights are sunrise.
##' @param p a two column matrix of (lon,lat) locations
##' @param zenith the solar zenith angle that defines twilight.
##' @return a vector twilight residuals (in minutes).
##' @export
twilight.residuals <- function(twilight,rise,p,zenith=96) {
  sgn <- ifelse(rise,1,-1)
  s <- solar(twilight)
  4*sgn*(s$solarTime-twilight.solartime(s,p[,1L],p[,2L],rise,zenith))
}



##' Estimate location from two consecutive twilights by a threshold
##' method
##'
##' This package and the \pkg{GeoLight} package provide some common
##' functionality, but are based on different astronomical
##' approximations.  This function provides a drop-in replacement for
##' the \code{coord} function from the \pkg{GeoLight} package to allow
##' easy comparison of the two systems.
##' @title Geolocation Estimation by the Threshold Method
##' @param tFirst factor or character vector representing the times of
##' the first twilight in the format "Y-m-d H:M:S".
##' @param tSecond factor or character vector representing the times
##' of the second twilight in the format "Y-m-d H:M:S".
##' @param type vector with elements 1 or 2, defining the elements of
##' \code{tFirst} as sunrise or sunset respectively.
##' @param degElevation sun elevation angle (90-zenith).
##' @return location estimates stored as a two column (lon,lat) matrix.
coord <- function(tFirst,tSecond,type,degElevation=-6) {
  tFirst <- as.POSIXct(tFirst,"GMT")
  tSecond <- as.POSIXct(tSecond,"GMT")
  rise <- ifelse(type==1,tFirst,tSecond)
  set <- ifelse(type==1,tSecond,tFirst)
  threshold.estimate(rise,set,zenith=90-degElevation)
}




##' Movement model that assumes speeds are Gamma distributed
##'
##' This function implements a movement model that assumes the speed
##' of travel between locations is Gamma distributed, and for Estelle
##' models, the change in bearing along the dog-leg path segments can
##' be assumed Normally distributed with mean zero.
##'
##' For Stella models, average speeds is calculated along great circle
##' paths between primary locations (x).  For Estelle, average speed
##' is calculated along dog leg paths through the intermediate points
##' (z), and the change in bearing at each intermediate point calculated.
##'
##' If \code{beta} is a vector, then \code{beta[1]} and \code{beta[2]}
##' specify the shape and rate of the Gamma distribution of speeds.
##' If \code{beta} has three elements, then \code{beta[3]} specifies
##' the standard deviation of the change in bearing (in degrees) along
##' dog leg paths.
##'
##' Alternately, these parameters can be specified individually for
##' each track segment by passing \code{beta} as a matrix with one row
##' for each segment.
##'
##' @title Gamma Behavioural Model
##' @param beta parameters of the behavioural model.
##' @param dt time intervals for speed calculation in hours.
##' @return Functions to evaluate the contributions to the log
##' posterior for the Estelle and Stella models.
##' @seealso \code{\link{satellite.model}}, \code{\link{threshold.model}},
##' \code{\link{grouped.threshold.model}}, \code{\link{curve.model}}.
##' @export
speed.gamma.model <- function(beta,dt) {

  ## Sanity check
  if(any(dt <= 0)) stop("Data not ordered in time")

  ## Ensure beta is always a matrix
  if(!is.matrix(beta)) beta <- t(beta)

  if(ncol(beta)==2) {
    ## Contribution to log posterior from the movement
    estelle.logpb <- function(x,z) {
      spd <- pmax.int(trackDist2(x,z), 1e-06)/dt
      dgamma(spd,beta[,1L],beta[,2L],log=TRUE)
    }
  }

  if(ncol(beta)==3) {
    ## Contribution to log posterior from the movement
    estelle.logpb <- function(x,z) {
      spd <- pmax.int(trackDist2(x,z), 1e-06)/dt
      angle <- trackBearingChange2(x,z)
      dgamma(spd,beta[,1L],beta[,2L],log=TRUE)+dnorm(angle,0,beta[,3L],log=TRUE)
    }
  }

  stella.logpb <- function(x) {
    spd <- pmax.int(trackDist(x), 1e-06)/dt
    dgamma(spd,beta[,1L],beta[,2L],log=TRUE)
  }


  ## Behavioural contribution to the log posterior
  list(
    estelle.logpb=estelle.logpb,
    stella.logpb=stella.logpb,
    beta=beta,
    dt=dt)
}


##' Satellite Model Structure for Stella and Estelle
##'
##' Stella and Estelle require a model structure that describes the
##' model being fitted. This function generates a basic model
##' structure for satellite data that should provide a suitable
##' starting point for most analyses.
##'
##' The \code{satellite.model} function constructs a model structure
##' assuming that assumes the locations \code{X} determined by the
##' satellite are independently distributed about the corresponding
##' true locations.  The \code{location.model} parameter selects
##' whether \code{(X-x)/sd} is
##'
##' \describe{
##' \item{'Normal'}{Normally distributed with zero mean and
##' unit variance, or}
##' \item{'T'}{t distributed with degrees of freedom df.}
##' }
##'
##' Both Estelle and Stella variants of the model assume that the
##' average speed of travel between successive locations is Gamma
##' distributed, and for Estelle models, the change in bearing
##' (degrees) along the dog-leg path segments can be assumed Normally
##' distributed with mean zero.  By default, the speed of travel is
##' calculated based on the time intervals between the twilights (in
##' hours), but the intervals of time actually available for travel
##' can be specified directly with the \code{dt} argument.
##'
##' If \code{beta} is a vector, then \code{beta[1]} and \code{beta[2]}
##' specify the shape and rate of the Gamma distribution of speeds.
##' If \code{beta} has three elements, then \code{beta[3]} specifies
##' the standard deviation of the change in bearing (in degrees) along
##' dog leg paths.
##'
##' The \code{satellite.model0} function constructs the non-movement
##' elements of the model, and the movement elements of the model are
##' constructed by the \code{speed.gamma.model} function.
##'
##' @title Satellite Model Structures
##' @param time the times of the satellite determined locations.
##' @param X the satellite determined locations.
##' @param location.model the model for the errors in satellite
##' locations.
##' @param sd a vector or two column matrix of dispersions for the
##' location model.
##' @param df a vector or two column matrix of degrees of freedom for
##' the t location model.
##' @param beta parameters of the behavioural model.
##' @param logp.x function to evaluate any additional contribution to
##' the log posterior from the satellite estimated locations.
##' @param logp.z function to evaluate any additional contribution to
##' the log posterior from the intermediate locations.
##' @param x0 suggested starting points for the satellite locations.
##' @param z0 suggested starting points for intermediate locations.
##' @param fixedx logical vector indicating which satellite locations
##' to hold fixed.
##' @param dt time intervals for speed calculation in hours.
##' @return a list with components
##' \item{\code{logpx}}{function to evaluate the contributions to the
##' log posterior from the twilight model}
##' \item{\code{logpz}}{function to evaluate the contributions to the
##' log posterior from the prior for the z locations}
##' \item{\code{estelle.logpb}}{function to evaluate contribution to
##' the log posterior from the behavioural model for estelle.}
##' \item{\code{stella.logpb}}{function to evaluate contribution to
##' the log posterior from the behavioural model for stella.}
##' \item{\code{fixedx}}{a logical vector indicating which locations
##' should remain fixed.}
##' \item{\code{x0}}{an array of initial twilight locations.}
##' \item{\code{z0}}{an array of initial intermediate locations.}
##' \item{\code{time}}{the times of the satellite determined
##' locations.}
##' \item{\code{X}}{the satellite determined locations.}
##' @export
satellite.model <- function(time,X,
                            location.model=c("Normal","T"),
                            sd,df=NULL,beta,
                            logp.x=function(x) rep.int(0L,nrow(x)),
                            logp.z=function(z) rep.int(0L,nrow(z)),
                            x0,z0=NULL,fixedx=FALSE,dt=NULL) {

  ## Times (hours) between observations
  if(is.null(dt)) dt <- diff(as.numeric(time)/3600)

  ## Contribution to log posterior from the x and z locations
  location.model <- satellite.model0(time,X,location.model,sd,df,logp.x,logp.z,fixedx)

  ## Contribution to log posterior from the movement model
  behavioural.model <- speed.gamma.model(beta,dt)

  c(location.model,
    behavioural.model,
    list(
      ## Suggested starting points
      x0=x0,z0=z0,
      ## Data
      time=time,
      X=X))
}


##' @rdname satellite.model
##' @export
satellite.model0 <- function(time,X,
                             location.model=c("Normal","T"),
                             sd,df=NULL,
                             logp.x=function(x) rep.int(0L,nrow(x)),
                             logp.z=function(z) rep.int(0L,nrow(z)),
                             fixedx=FALSE) {

  ## Fixed x locations
  fixedx <- rep_len(fixedx,length.out=length(time))

  ## Function to compute residuals
  residuals <- function(x) X-x

  ## Contribution to log posterior from each x location
  location.model <- match.arg(location.model)
  logpx <-
    switch(location.model,
           Normal=
           function(x) {
             r <- X-x
             logp <- rowSums(dnorm(r,0,sd,log=TRUE)) + logp.x(x)
             logp[fixedx] <- 0
             logp
           },
           T=
           function(x) {
             r <- X-x
             logp <- rowSums(dt(r/sd,df,log=TRUE)) + logp.x(x)
             logp[fixedx] <- 0
             logp
           })

  ## Contribution to log posterior from each z location
  logpz <- logp.z

  list(
    ## Positional contribution to the log posterior
    logpx=logpx,
    logpz=logpz,
    ## Residuals
    residuals=residuals,
    ## Locations to be held fixed
    fixedx=fixedx)
}


##' Log density of twilight errors
##'
##' Construct a function to evalute the log density of the twilight
##' errors in a threshold model.
##'
##' One of several models models may be selected for the errors in
##' twilight times.  The errors in twilight time are defined as the
##' difference in the observed and true times of twilight, with sign
##' selected so that a positive error always corresponds to a sunrise
##' observed after the true time of sunrise, and sunset observed
##' before the true time of sunset. That is, a positive error
##' corresponds to the observed light level being lower than expected.
##'
##' The properties of the twilight model are determined by
##' \code{alpha}, which must be either a vector of parameters that are
##' to be applied to each twilight, or a matrix of parameters with one
##' row for each twilight.
##'
##' The \code{twilight.model} argument selects the distribution of the
##' twilight errors
##' \describe{
##' \item{'Normal'}{Normally distributed with mean \code{alpha[,1]} and
##' standard deviation \code{alpha[,2]},}
##' \item{'LogNormal'}{Log Normally distributed so the log errors have
##' mean \code{alpha[,1]} and standard deviation \code{alpha[,2]}, or}
##' \item{'Gamma'}{Gamma distributed with shape \code{alpha[,1]} and
##' rate \code{alpha[2]}.}
##' }
##' The 'LogNormal' and 'Gamma' models forbid negative errors, that
##' is, the observed light cannot be brighter than expected.  There
##' are modified variants of these models for which negative errors
##' are extremely unlikely, but not forbidden, and can be used to
##' generate suitable initialization locations for their unmodified
##' counterparts.
##'
##' @title Twilight Error Models
##' @param twilight.model the model for the errors in twilight times.
##' @param alpha parameters of the twilight model.
##' @return a function to evaluate the log density of the twilight
##' residuals in a threshold model.
##' @export
make.twilight.model <- function(twilight.model=c("Gamma","LogNormal","Normal","ModifiedGamma","ModifiedLogNormal"),
                                alpha) {

  twilight.model <- match.arg(twilight.model)
  logp <- switch(twilight.model,
                 Gamma={
                   shape <- alpha[,1L]
                   rate <- alpha[,2L]
                   function(r) dgamma(r,shape,rate,log=TRUE)
                 },
                 LogNormal={
                   meanlog <- alpha[,1L]
                   sdlog <- alpha[,2L]
                   function(r) dlnorm(r,meanlog,sdlog,log=TRUE)
                 },
                 Normal={
                   mean <- alpha[,1L]
                   sd <- alpha[,2L]
                   function(r) dnorm(r,mean,sd,log=TRUE)
                 },
                 ModifiedGamma={
                   shape <- alpha[,1L]
                   rate <- alpha[,2L]
                   function(r)
                     ifelse(is.finite(r) & r < 0,
                            60*r-1.0E8+dgamma(shape/rate,shape,rate,log=TRUE),
                            dgamma(r,shape,rate,log=TRUE))
                 },
                 ModifiedLogNormal={
                   meanlog <- alpha[,1L]
                   sdlog <- alpha[,2L]
                   function(r)
                     ifelse(is.finite(r) & r < 0,
                            60*r-1.0E8+dlnorm(exp(meanlog+sdlog^2/2),meanlog,sdlog,log=TRUE),
                            dlnorm(r,meanlog,sdlog,log=TRUE))
                 })
  logp
}


##' Threshold Model Structures for Stella and Estelle
##'
##' Stella and Estelle require a model structure that describes the
##' model being fitted. These function generate basic model structures
##' for threshold twilight data that should provide a suitable
##' starting point for most analyses.
##'
##' The \code{threshold.model} function constructs a model structure
##' assuming that each twilight time is associated with a single
##' location, while the \code{grouped.threshold.model} function allows
##' multiple twilight times to be associated with a single location.
##'
##' One of several models models may be selected for the errors in
##' twilight times.  The errors in twilight time are defined as the
##' difference in the observed and true times of twilight, with sign
##' selected so that a positive error always corresponds to a sunrise
##' observed after the true time of sunrise, and sunset observed
##' before the true time of sunset. That is, a positive error
##' corresponds to the observed light level being lower than expected.
##'
##' The properties of the twilight model are determined by
##' \code{alpha}, which must be either a vector of parameters that are
##' to be applied to each twilight, or a matrix of parameters with one
##' row for each twilight.
##'
##' The \code{twilight.model} argument selects the distribution of the
##' twilight errors
##' \describe{
##' \item{'Normal'}{Normally distributed with mean \code{alpha[,1]} and
##' standard deviation \code{alpha[,2]},}
##' \item{'LogNormal'}{Log Normally distributed so the log errors have
##' mean \code{alpha[,1]} and standard deviation \code{alpha[,2]}, or}
##' \item{'Gamma'}{Gamma distributed with shape \code{alpha[,1]} and
##' rate \code{alpha[2]}.}
##' }
##' The 'LogNormal' and 'Gamma' models forbid negative errors, that
##' is, the observed light cannot be brighter than expected.  There
##' are modified variants of these models for which negative errors
##' are extremely unlikely, but not forbidden, and can be used to
##' generate suitable initialization locations for their unmodified
##' counterparts.
##'
##' The initialization locations \code{x0} and \code{z0} must be
##' consistent with the chosen twilight model.  That is, if
##' 'LogNormal' or 'Gamma' models are selected, the \code{x0} cannot
##' yield negative twilight errors.
##'
##' Both Estelle and Stella variants of the model assume that the
##' average speed of travel between successive locations is Gamma
##' distributed, and for Estelle models, the change in bearing
##' (degrees) along the dog-leg path segments can be assumed Normally
##' distributed with mean zero.  By default, the speed of travel is
##' calculated based on the time intervals between the twilights (in
##' hours), but the intervals of time actually available for travel
##' can be specified directly with the \code{dt} argument.
##'
##' If \code{beta} is a vector, then \code{beta[1]} and \code{beta[2]}
##' specify the shape and rate of the Gamma distribution of speeds.
##' If \code{beta} has three elements, then \code{beta[3]} specifies
##' the standard deviation of the change in bearing (in degrees) along
##' dog leg paths.
##'
##' Twilights can be missing because either the light record was too
##' noisy at that time to estimate twilight reliably, or because the
##' tag was at very high latitude and no twilight was observed.
##' Missing twilights should be replaced with an approximate time of
##' twilight, and the vector \code{missing} used to indicate which
##' twilights are approaximate and which are true.  This should be a
##' vector of integers, one for each twilight where the integer codes
##' signify
##' \describe{
##' \item{0:}{The twilight is not missing.}
##' \item{1:}{The twilight is missing, but a twilight did occur.}
##' \item{2:}{The twilight is missing because twilight did not occur.}
##' \item{3:}{The twilight is missing and it is not known if a twilight occurred.}
##' }
##'
##' The \code{threshold.model0} and \code{grouped.threshold.model0}
##' functions construct the non-movement elements of the model, and
##' the movement elements of the model are constructed by the
##' \code{speed.gamma.model} function.
##'
##' @title Threshold Model Structures
##' @param twilight the observed times of twilight as POSIXct.
##' @param rise logical vector indicating which twilights are sunrise.
##' @param group integer vector that defines the twilight groups.  If
##' code{group[k]==j} then the k-th twilight occurs at location j.
##' @param twilight.model the model for the errors in twilight times.
##' @param alpha parameters of the twilight model.
##' @param beta parameters of the behavioural model.
##' @param logp.x function to evaluate any additional contribution to
##' the log posterior from the twilight locations.
##' @param logp.z function to evaluate any additional contribution to
##' the log posterior from the intermediate locations.
##' @param x0 suggested starting points for twilight locations.
##' @param z0 suggested starting points for intermediate locations.
##' @param fixedx logical vector indicating which twilight locations
##' to hold fixed.
##' @param missing integer vector indicating which twilights were
##' unobserved and why.
##' @param dt time intervals for speed calculation in hours.
##' @param zenith the solar zenith angle that defines twilight.
##' @return a list with components
##' \item{\code{logpx}}{function to evaluate the contributions to the
##' log posterior from the twilight model}
##' \item{\code{logpz}}{function to evaluate the contributions to the
##' log posterior from the prior for the z locations}
##' \item{\code{estelle.logpb}}{function to evaluate contribution to
##' the log posterior from the behavioural model for estelle.}
##' \item{\code{stella.logpb}}{function to evaluate contribution to
##' the log posterior from the behavioural model for stella.}
##' \item{\code{residuals}}{function to evaluate the twilight model
##' residuals.}
##' \item{\code{fixedx}}{a logical vector indicating which locations
##' should remain fixed.}
##' \item{\code{x0}}{an array of initial twilight locations.}
##' \item{\code{z0}}{an array of initial intermediate locations.}
##' \item{\code{time}}{the twilight times.}
##' \item{\code{rise}}{the sunrise indicators.}
##' \item{\code{group}}{the grouping vector.}
##' @export
threshold.model <- function(twilight,rise,
                            twilight.model=c("Gamma","LogNormal","Normal","ModifiedGamma","ModifiedLogNormal"),
                            alpha,beta,
                            logp.x=function(x) rep.int(0L,nrow(x)),
                            logp.z=function(z) rep.int(0L,nrow(z)),
                            x0,z0=NULL,fixedx=FALSE,missing=0,dt=NULL,zenith=96) {

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(twilight)/3600)

  ## Contribution to log posterior from the x and z locations
  location.model <- threshold.model0(twilight,rise,twilight.model,alpha,logp.x,logp.z,fixedx,missing,zenith)

  ## Contribution to log posterior from the movement model
  behavioural.model <- speed.gamma.model(beta,dt)

  c(location.model,
    behavioural.model,
    list(
      ## Suggested starting points
      x0=x0,z0=z0,
      ## Data
      time=twilight,
      rise=rise))
}

##' @rdname threshold.model
##' @export
threshold.model0 <- function(twilight,rise,
                             twilight.model=c("Gamma","LogNormal","Normal","ModifiedGamma","ModifiedLogNormal"),
                             alpha,
                             logp.x=function(x) rep.int(0L,nrow(x)),
                             logp.z=function(z) rep.int(0L,nrow(z)),
                             fixedx=FALSE,missing=0,zenith=96) {

  ## Convert twilights to solar time.
  s <- solar(twilight)
  ## Sign for residuals
  sgn <- ifelse(rise,1,-1)
  ## Fixed x locations
  fixedx <- rep_len(fixedx,length.out=length(twilight))

  ## Ensure alpha is alway a matrix
  if(!is.matrix(alpha)) alpha <- t(alpha)

  ## Discrepancy in expected and observed times of twilight, with sign
  ## selected so that a positive value corresponds to the observed
  ## sunrise occurring after the expected time of sunrise, and the
  ## observed sunset occurring before the expected time of sunset
  residuals <- function(x) {
    4*sgn*(s$solarTime-twilight.solartime(s,x[,1L],x[,2L],rise,zenith))
  }

  ## Select the density of the twilight residuals
  twilight.model <- match.arg(twilight.model)
  logp.residual <- make.twilight.model(twilight.model,alpha)

  ## Log density at forbidden locations
  forbid <- if(twilight.model %in% c("ModifiedGamma","ModifiedLogNormal")) -1.0E8 else -Inf

  ## Contribution to log posterior from each x location
  if(any(missing > 0)) {
    logpx <- function(x) {
      r <- residuals(x)
      logp <- logp.residual(r)
      ## Fix missing values
      logp[missing <= 1 & !is.finite(r)] <- forbid
      logp[missing == 1 & is.finite(r)] <- 0
      logp[missing == 2 & is.finite(r)] <- forbid
      logp[missing == 2 & !is.finite(r)] <- 0
      logp[missing == 3] <- 0
      logp <- logp + logp.x(x)
      logp[fixedx] <- 0
      logp
    }
  } else {
    logpx <- function(x) {
      r <- residuals(x)
      logp <- logp.residual(r)
      ## Fix missing values
      logp[!is.finite(r)] <- forbid
      logp <- logp + logp.x(x)
      logp[fixedx] <- 0
      logp
    }
  }

  ## Contribution to log posterior from each z location
  logpz <- logp.z

  list(
    ## Positional contribution to the log posterior
    logpx=logpx,
    logpz=logpz,
    ## Residuals
    residuals=residuals,
    ## Locations to be held fixed
    fixedx=fixedx)
}





##' @rdname threshold.model
##' @export
grouped.threshold.model <- function(twilight,rise,group,
                                    twilight.model=c("Gamma","LogNormal","Normal","ModifiedGamma","ModifiedLogNormal"),
                                    alpha,beta,
                                    logp.x=function(x) rep.int(0L,nrow(x)),
                                    logp.z=function(z) rep.int(0L,nrow(z)),
                                    x0,z0=NULL,fixedx=FALSE,missing=0,dt=NULL,zenith=96) {

  ## Compute median times
  time <- .POSIXct(tapply(twilight,group,median),"GMT")

  if(is.null(dt)) {
    ## Times (hours) between twilight groups
    tmin <- tapply(as.numeric(twilight)/3600,group,min)
    tmax <- tapply(as.numeric(twilight)/3600,group,max)
    dt <- tmin[-1L]-tmax[-max(group)]
  }

  ## Contribution to log posterior from the x and z locations
  location.model <- grouped.threshold.model0(twilight,rise,group,twilight.model,alpha,logp.x,logp.z,fixedx,missing,zenith)

  ## Contribution to log posterior from the movement model
  behavioural.model <- speed.gamma.model(beta,dt)

  c(location.model,
    behavioural.model,
    list(
      ## Suggested starting points
      x0=x0,z0=z0,
      ## Data
      twilight=twilight,
      rise=rise,
      group=group,
      time=time))
}


##' @rdname threshold.model
##' @export
grouped.threshold.model0 <- function(twilight,rise,group,
                                     twilight.model=c("Gamma","LogNormal","Normal","ModifiedGamma","ModifiedLogNormal"),
                                     alpha,
                                     logp.x=function(x) rep.int(0L,nrow(x)),
                                     logp.z=function(z) rep.int(0L,nrow(z)),
                                     fixedx=FALSE,missing=0,zenith=96) {

  ## Convert twilights to solar time.
  s <- solar(twilight)
  ## Sign for residuals
  sgn <- ifelse(rise,1,-1)
  ## Fixed x locations
  fixedx <- rep_len(fixedx,length.out=max(group))

  ## Ensure alpha is always a matrix
  if(!is.matrix(alpha)) alpha <- t(alpha)

  ## Discrepancy in expected and observed times of twilight, with sign
  ## selected so that a positive value corresponds to the observed
  ## sunrise occurring after the expected time of sunrise, and the
  ## observed sunset occurring before the expected time of sunset
  residuals <- function(x) {
    4*sgn*(s$solarTime-twilight.solartime(s,x[group,1L],x[group,2L],rise,zenith))
  }

  ## Select the density of the twilight residuals
  twilight.model <- match.arg(twilight.model)
  logp.residual <- make.twilight.model(twilight.model,alpha)

  ## Log density at forbidden locations
  forbid <- if(twilight.model %in% c("ModifiedGamma","ModifiedLogNormal")) -1.0E8 else -Inf

  ## Contribution to log posterior from each x location
  if(any(missing > 0)) {
    logpx <- function(x) {
      r <- residuals(x)
      logp <- logp.residual(r)
      ## Fix missing values
      logp[missing <= 1 & !is.finite(r)] <- forbid
      logp[missing == 1 & is.finite(r)] <- 0
      logp[missing == 2 & is.finite(r)] <- forbid
      logp[missing == 2 & !is.finite(r)] <- 0
      logp[missing == 3] <- 0
      logp <- tapply(logp,group,sum)+logp.x(x)
      logp[fixedx] <- 0
      logp
    }
  } else {
    logpx <- function(x) {
      r <- residuals(x)
      logp <- logp.residual(r)
      ## Fix missing values
      logp[!is.finite(r)] <- forbid
      logp <- tapply(logp,group,sum)+logp.x(x)
      logp[fixedx] <- 0
      logp
    }
  }

  ## Contribution to log posterior from each z location
  logpz <- logp.z

  list(## Positional contribution to the log posterior
    logpx=logpx,
    logpz=logpz,
    ## Residuals
    residuals=residuals,
    ## Locations to be held fixed
    fixedx=fixedx)
}




##' Light Curve Model Structures for Stella and Estelle
##'
##' Stella and Estelle require a model structure that describes the
##' model being fitted. These function generate basic model structures
##' for curve fitting methods that should provide a suitable
##' starting point for most analyses.
##'
##' The \code{curve.model} function constructs a model structure
##' assuming that each twilight profile is associated with a single
##' location.  The errors in observed log light level are assumed to
##' be Normally distributed about their expected value.
##'
##' The initialization locations \code{x0} and \code{z0} must be
##' consistent with any other constraints imposed by the data.
##'
##' Both Estelle and Stella variants of the model assume that the
##' average speed of travel between successive locations is Gamma
##' distributed, and for Estelle models, the change in bearing
##' (degrees) along the dog-leg path segments can be assumed Normally
##' distributed with mean zero.  By default, the speed of travel is
##' calculated based on the time intervals between the twilights (in
##' hours), but the intervals of time actually available for travel
##' can be specified directly with the \code{dt} argument.
##'
##' If \code{beta} is a vector, then \code{beta[1]} and \code{beta[2]}
##' specify the shape and rate of the Gamma distribution of speeds.
##' If \code{beta} has three elements, then \code{beta[3]} specifies
##' the standard deviation of the change in bearing (in degrees) along
##' dog leg paths.
##'
##' The \code{curve.model0} function constructs only the non-movement
##' elements of the model, and can be used as a basis for those
##' wishing to experiment with alternative movement models.
##'
##' @title Curve Model Structure
##' @param time vector of sample times as POSIXct.
##' @param light vector of observed (log) light levels.
##' @param segment vector of integers that assign observations to
##' twilight segments.
##' @param calibration function that maps zenith angles to expected
##' light levels.
##' @param alpha parameters of the twilight model.
##' @param beta parameters of the behavioural model.
##' @param logp.x function to evaluate any additional contribution to
##' the log posterior from the twilight locations.
##' @param logp.z function to evaluate any additional contribution to
##' the log posterior from the intermediate locations.
##' @param x0 suggested starting points for twilight locations.
##' @param z0 suggested starting points for intermediate locations.
##' @param fixedx logical vector indicating which twilight locations
##' to hold fixed.
##' @param dt time intervals for speed calculation in hours.
##' @return a list with components
##' \item{\code{logpx}}{function to evaluate the contributions to the
##' log posterior from the twilight model}
##' \item{\code{logpz}}{function to evaluate the contributions to the
##' log posterior from the prior for the z locations}
##' \item{\code{estelle.logpb}}{function to evaluate contribution to
##' the log posterior from the behavioural model for estelle.}
##' \item{\code{stella.logpb}}{function to evaluate contribution to
##' the log posterior from the behavioural model for stella.}
##' \item{\code{fitted}}{function to evaluate the fitted values for a
##' given set of locations.}
##' \item{\code{fixedx}}{a logical vector indicating which locations
##' should remain fixed.}
##' \item{\code{x0}}{an array of initial twilight locations.}
##' \item{\code{z0}}{an array of initial intermediate locations.}
##' \item{\code{time}}{the sample times.}
##' \item{\code{light}}{the recorded light levels.}
##' \item{\code{segment}}{vector of integers that assign observations
##' to twilight segments.}
##' @export
curve.model <- function(time,light,segment,
                        calibration,alpha,beta,
                        logp.x=function(x) rep.int(0L,nrow(x)),
                        logp.z=function(z) rep.int(0L,nrow(z)),
                        x0=NULL,z0=NULL,fixedx=FALSE,dt=NULL) {

  ## Median time in each segment
  tm <- .POSIXct(sapply(split(time,segment),median),"GMT")

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(tm)/3600)

  ## Contribution to log posterior from the x and z locations
  location.model <- curve.model0(time,light,segment,calibration,alpha,logp.x,logp.z,fixedx)

  ## Contribution to log posterior from the movement model
  behavioural.model <- speed.gamma.model(beta,dt)

  c(location.model,
    behavioural.model,
    list(
      ## Suggested starting points
      x0=x0,z0=z0,
      ## Median time of twilights
      time=tm))
}



##' @rdname curve.model
##' @export
curve.model0 <- function(time,light,segment,
                         calibration,alpha,
                         logp.x=function(x) rep.int(0L,nrow(x)),
                         logp.z=function(z) rep.int(0L,nrow(z)),
                         fixedx=FALSE) {

  ## Convert to solar time.
  sun <- solar(time)
  ## Fixed x locations
  fixedx <- rep_len(fixedx,length.out=nlevels(factor(segment)))

  ## Ensure alpha is alway a matrix
  if(!is.matrix(alpha)) alpha <- t(alpha)

  ## Return a dataframe of the fitted zenith and light observations
  ## for a set of estimated locations.
  fitted <- function(x) {
    xs <- x[segment,]
    zenith <- zenith(sun,xs[,1L],xs[,2L])
    fitted <- calibration(zenith)+xs[,3L]
    ## Contributions to log posterior
    logp <- dnorm(light,fitted,alpha[,1L],log=TRUE)
    data.frame(Time=time,
               Segment=segment,
               Zenith=zenith,
               Fitted=fitted,
               Light=light,
               LogL=logp)
  }

  ## Contribution to log posterior from each x location
  logpx <- function(x) {
    xs <- x[segment,]
    zenith <- zenith(sun,xs[,1L],xs[,2L])
    fitted <- calibration(zenith)+xs[,3L]
    ## Contributions to log posterior
    logp <- dnorm(light,fitted,alpha[,1L],log=TRUE)
    sapply(split(logp,segment),sum) + logp.x(x) + dnorm(x[,3L],0,alpha[,2L],log=TRUE)
  }

  ## Contribution to log posterior from each z location
  logpz <- logp.z

  list(
    ## Positional contribution to the log posterior
    logpx=logpx,
    logpz=logpz,
    ## Fitted values
    fitted=fitted,
    ## Locations to be held fixed
    fixedx=fixedx)
}



##' Metropolis samplers for Stella or Estelle.
##'
##' These functions draw samples form posterior for the Stella or
##' Estelle model by the Metropolis algorithm.
##' @title Metropolis Samplers
##' @param model a model structure as generated by
##' \code{threshold.model}.
##' @param proposal.x function for drawing proposals for x.
##' @param proposal.z function for drawing proposals for z.
##' @param x0 Starting values for twilight locations x.
##' @param z0 Starting values for intermediate locations z.
##' @param iters number of samples to draw.
##' @param thin rate at which to thin samples.
##' @param chains number of chains to sample.
##' @param verbose report progress at prompt?
##' @return If there are r samples drawn for each of q chains of p
##' parameters at n locations, Stella will return a list containing
##' \item{\code{model}}{the model structure}
##' \item{\code{x}}{a list of n x p x r arrays of twilight locations
##' from the q chains}
##' While in addition Estelle will return
##' \item{\code{z}}{a list of (n-1) x p x r arrays of intermediate
##' locations from the q chains}.
##' @seealso \code{\link{threshold.model}}
##' @export
estelle.metropolis <- function(model,
                               proposal.x,proposal.z,
                               x0=NULL,z0=NULL,
                               iters=1000L,thin=10L,chains=1L,
                               verbose=interactive()) {

  ## Initialize x,z
  if(is.null(x0)) x0 <- model$x0
  if(is.null(z0)) z0 <- model$z0
  ## Expand starting values for multiple chains
  x0 <- rep(if(is.list(x0)) x0 else list(x0),length.out=chains)
  z0 <- rep(if(is.list(z0)) z0 else list(z0),length.out=chains)

  ## Number of locations
  n <- nrow(x0[[1]])
  ## Number of parameters
  m <- ncol(x0[[1]])

  ## Extract model components
  logpx <- model$logpx
  logpz <- model$logpz
  logpb <- model$estelle.logpb
  fixedx <- model$fixedx

  ## Lists of chains
  ch.xs <- vector(mode="list",chains)
  ch.zs <- vector(mode="list",chains)

  ## PARALLEL - parallelise this loop
  for(k1 in 1:chains) {
    ## Allocate chains
    ch.x <- array(0,c(n,m,iters))
    ch.z <- array(0,c(n-1L,2L,iters))

    ## Initialize
    x1 <- x0[[k1]]
    z1 <- z0[[k1]]
    ## Drop dimnames for speed
    dimnames(x1) <- NULL
    dimnames(z1) <- NULL

    ## Contribution to logp from the initial x,z
    logp.x1 <- logpx(x1)
    logp.z1 <- logpz(z1)
    logp.b1 <- logpb(x1,z1)

    k2 <- 0
    if(verbose) {
      cat("iter ",sprintf("%6d",k2))
      flush.console()
    }

    for(k2 in 1:iters) {

      if(verbose && k2%%10==0) {
        cat("\b\b\b\b\b\b");
        cat(sprintf("%6d",k2));
        flush.console()
      }

      for(k3 in 1:thin) {

        ## Propose all x at once, and calculate contribution to the log
        ## posterior
        x2 <- proposal.x(x1)
        x2[fixedx,] <- x1[fixedx,]
        logp.x2 <- logpx(x2)

        x <- x1
        x[c(1L,n),] <- x2[c(1L,n),]
        logp.b2 <- logpb(x,z1)


        ## Update x
        ## In each case we compute full contribution (positional +
        ## behavourial) to the log posterior for current and proposed
        ## points, and apply the MH rule. If the proposal is accepted,
        ## we update both x and the cached contributions to the log
        ## posterior.


        ## Accept/reject first x
        if(!fixedx[1L]) {
          logp1 <- logp.x1[1L]+logp.b1[1L]
          logp2 <- logp.x2[1L]+logp.b2[1L]
          if(logp2-logp1 > log(runif(1))) {
            x1[1L,] <- x2[1L,]
            logp.x1[1L] <- logp.x2[1L]
            logp.b1[1L] <- logp.b2[1L]
          }
        }


        ## Accept/reject last x
        if(!fixedx[n]) {
          logp1 <- logp.x1[n]+logp.b1[n-1L]
          logp2 <- logp.x2[n]+logp.b2[n-1L]
          if(logp2-logp1 > log(runif(1))) {
            x1[n,] <- x2[n,]
            logp.x1[n] <- logp.x2[n]
            logp.b1[n-1L] <- logp.b2[n-1L]
          }
        }


        ## Red/Black update for interior x
        for(rb in 2:3) {
          is <- seq.int(rb,n-1L,by=2L)
          x <- x1
          x[is,] <- x2[is,]
          logp.b2 <- logpb(x,z1)

          logp1 <- logp.x1[is]+logp.b1[is-1L]+logp.b1[is]
          logp2 <- logp.x2[is]+logp.b2[is-1L]+logp.b2[is]
          ## MH rule - compute indices of the accepted points.
          accept <- is[logp2-logp1 > log(runif(length(is)))]
          x1[accept,] <- x[accept,]
          logp.x1[accept] <- logp.x2[accept]
          logp.b1[accept] <- logp.b2[accept]
          logp.b1[accept-1L] <- logp.b2[accept-1L]
        }

        ## Update z
        ## Here we need only consider the behavioural contributions to
        ## the log posterior (the position contributions are constant
        ## and would cancel), and so we can update all the z in parallel.
        z2 <- proposal.z(z1)
        logp.z2 <- logpz(z2)
        logp.b2 <- logpb(x1,z2)
        logp1 <- logp.z1+logp.b1
        logp2 <- logp.z2+logp.b2
        ## MH rule - compute indices of the accepted points.
        accept <- logp2-logp1 > log(runif(n-1L))
        z1[accept,] <- z2[accept,]
        logp.z1[accept] <- logp.z2[accept]
        logp.b1[accept] <- logp.b2[accept]
      }
      ## Store the current state
      ch.x[,,k2] <- x1
      ch.z[,,k2] <- z1
    }
    ch.xs[[k1]] <- ch.x
    ch.zs[[k1]] <- ch.z
    if(verbose) cat("\n")
  }
  list(model=model,x=ch.xs,z=ch.zs)
}



##' @rdname estelle.metropolis
##' @export
stella.metropolis <- function(model,
                              proposal.x,
                              x0=NULL,
                              iters=1000L,thin=10L,chains=1L,
                              verbose=interactive()) {

  ## Initialize x
  if(is.null(x0)) x0 <- model$x0
  ## Expand starting values for multiple chains
  x0 <- rep(if(is.list(x0)) x0 else list(x0),length.out=chains)

  ## Number of locations
  n <- nrow(x0[[1]])
  ## Number of parameters
  m <- ncol(x0[[1]])

  ## Extract model components
  logpx <- model$logpx
  logpb <- model$stella.logpb
  fixedx <- model$fixedx

  ## Lists of chains
  ch.xs <- vector(mode="list",chains)

  ## PARALLEL - parallelise this loop
  for(k1 in 1:chains) {
    ## Allocate chain
    ch.x <- array(0,c(n,m,iters))

    ## Initialize
    x1 <- x0[[k1]]
    ## Drop dimnames for speed
    dimnames(x1) <- NULL

    ## Contribution to logp from the initial x
    logp.x1 <- logpx(x1)
    logp.b1 <- logpb(x1)

    k2 <- 0
    if(verbose) {
      cat("iter ",sprintf("%6d",k2))
      flush.console()
    }

    for(k2 in 1:iters) {

      if(verbose && k2%%10==0) {
        cat("\b\b\b\b\b\b");
        cat(sprintf("%6d",k2));
        flush.console()
      }

      for(k3 in 1:thin) {

        ## Propose all x at once, and calculate contribution to the log
        ## posterior
        x2 <- proposal.x(x1)
        x2[fixedx,] <- x1[fixedx,]
        logp.x2 <- logpx(x2)

        x <- x1
        x[c(1L,n),] <- x2[c(1L,n),]
        logp.b2 <- logpb(x)


        ## Update x
        ## In each case we compute full contribution (positional +
        ## behavourial) to the log posterior for current and proposed
        ## points, and apply the MH rule. If the proposal is accepted,
        ## we update both x and the cached contributions to the log
        ## posterior.


        ## Accept/reject first x
        if(!fixedx[1L]) {
          logp1 <- logp.x1[1L]+logp.b1[1L]
          logp2 <- logp.x2[1L]+logp.b2[1L]
          if(logp2-logp1 > log(runif(1))) {
            x1[1L,] <- x2[1L,]
            logp.x1[1L] <- logp.x2[1L]
            logp.b1[1L] <- logp.b2[1L]
          }
        }


        ## Accept/reject last x
        if(!fixedx[n]) {
          logp1 <- logp.x1[n]+logp.b1[n-1L]
          logp2 <- logp.x2[n]+logp.b2[n-1L]
          if(logp2-logp1 > log(runif(1))) {
            x1[n,] <- x2[n,]
            logp.x1[n] <- logp.x2[n]
            logp.b1[n-1L] <- logp.b2[n-1L]
          }
        }


        ## Red/Black update for interior x
        for(rb in 2:3) {
          is <- seq.int(rb,n-1L,by=2L)
          x <- x1
          x[is,] <- x2[is,]
          logp.b2 <- logpb(x)

          logp1 <- logp.x1[is]+logp.b1[is-1L]+logp.b1[is]
          logp2 <- logp.x2[is]+logp.b2[is-1L]+logp.b2[is]
          ## MH rule - compute indices of the accepted points.
          accept <- is[logp2-logp1 > log(runif(length(is)))]
          x1[accept,] <- x[accept,]
          logp.x1[accept] <- logp.x2[accept]
          logp.b1[accept] <- logp.b2[accept]
          logp.b1[accept-1L] <- logp.b2[accept-1L]

        }
      }
      ## Store the current state
      ch.x[,,k2] <- x1
    }
    ch.xs[[k1]] <- ch.x
    if(verbose) cat("\n")
  }
  list(model=model,x=ch.xs)
}



##' Number of locations
##'
##' A convience function to determine the number of locations a chain,
##' or set of initial locations or a location summary
##' describe. Assumes \code{s} is either an array or a list of arrays
##' in which the first dimension corresponds to location, and returns
##' the length of the first dimension.
##'
##' @title Number of locations
##' @param s an array or a list of arrays.
##' @return size of the first dimension of the array.
##' @export
nlocation <- function(s) {
  dim(if(is.list(s)) s[[1]] else s)[1]
}


##' Summarize a set of location samples
##'
##' These functions compute various summaries of a sample or list of
##' samples generated by \code{estelle.metropolis} or
##' \code{stella.metropolis}.
##'
##' These functions accept either a sample from a single mcmc run, or
##' a list of samples from parallel mcmc runs.  When \code{collapse}
##' is true, multiple samples are merged and single result is
##' returned, otherwise a result is returned for each sample.
##'
##' @rdname location.summary
##' @title Summaries of Location Samples
##' @param s a single chain or a list of parallel chains generated by
##' \code{estelle.metropolis} or \code{stella.metropolis}.
##' @param time the times corresponding to the x locations.
##' @param discard number of initial samples to discard.
##' @param alpha coverage of the credible intervals calculated by
##' \code{location.summary}.
##' @param collapse whether to collapse parallel chains to a single chain
##' @param chains the set of chains to retain, or \code{NULL}.
##' @return
##' \item{\code{location.summary}}{returns a dataframe or a list of
##' dataframes of summary quantities for each location.}
##' \item{\code{location.mean}}{returns an array or a list of arrays
##' of the means of the samples for each location.}
##' @export
location.summary <- function(s,time=NULL,discard=0,alpha=0.95,collapse=TRUE,chains=NULL) {
  summary <- function(s) {
     stat <- function(x) c(mean=mean(x),sd=sd(x),quantile(x,prob=c(0.5,(1-alpha)/2,1-(1-alpha)/2)))
    lon <- t(apply(s[,1L,],1L,stat))
    colnames(lon) <- paste("Lon",colnames(lon),sep=".")
    lat <- t(apply(s[,2L,],1L,stat))
    colnames(lat) <- paste("Lat",colnames(lat),sep=".")
    d <- as.data.frame(cbind(lon,lat))
    if(!is.null(time)) {
      ## Add timing information
      n <- nrow(d)
      if(length(time)==n)
        d <- cbind(Time=time,d)
      else
        d <- cbind(Time1=time[1:n],Time2=time[2:(n+1L)],d)
    }
    d
   }

  s <- chain.collapse(s,collapse=collapse,discard=discard,chains=chains)
  if(is.list(s)) lapply(s,summary) else summary(s)
}

##' @rdname location.summary
##' @export
location.mean <- function(s,discard=0,collapse=TRUE,chains=NULL) {
  locmean <- function(s) apply(s[,1:2,],1:2,mean)

  s <- chain.collapse(s,collapse=collapse,discard=discard,chains=chains)
  if(is.list(s)) lapply(s,locmean) else locmean(s)
}



##' Bin locations to form a 2D density image
##'
##' Bins the samples for a sequence of locations to produce 2D array
##' suitable for plotting with \code{image}.  Optionally, a vector of
##' weights can be provided to differentially weight samples by
##' location.
##'
##' This function accepts either a sample from a single mcmc run, or
##' a list of samples from parallel mcmc runs.  If
##' \code{collapse} is true, multiple samples are merged and single
##' image is returned, otherwise an image is returned for each sample.
##'
##' @title Location Density Image
##' @param s a single chain or a list of parallel chains generated by
##' \code{estelle.metropolis} or \code{stella.metropolis}.
##' @param xlim range of the first coordinate.
##' @param ylim range of the second coordinate.
##' @param nx number of cells in the first coordinate.
##' @param ny number of cells in the second coordinate.
##' @param weight weights for each location.
##' @param discard number of initial samples to discard.
##' @param collapse whether to collapse parallel chains to a single chain
##' @param chains the set of chains to retain, or \code{NULL}.
##' @return A list with elesments
##' \item{\code{x}}{the x-ordinates that bound the bins}
##' \item{\code{y}}{the y-ordinates that bound the bins}
##' \item{\code{W}}{the weighted image.}
##' @export
location.image <- function(s,xlim,ylim,nx,ny,weight=rep_len(1,dim(s)[1L]),
                           discard=0,collapse=TRUE,chains=NULL) {
  nx <- round(nx)
  ny <- round(ny)
  xbin <- seq.int(xlim[1L],xlim[2L],length.out=nx+1L)
  ybin <- seq.int(ylim[1L],ylim[2L],length.out=ny+1L)

  bin <- function(s) {
    W <- 0
    for(k in 1:dim(s)[1L]) {
      W <- W+weight[k]*table(
        factor(.bincode(s[k,1L,],xbin),levels=1:nx),
        factor(.bincode(s[k,2L,],ybin),levels=1:ny))
    }
    W[W==0] <- NA
    list(x=xbin,y=ybin,W=W)
  }

  s <- chain.collapse(s,collapse=collapse,discard=discard,chains=chains)
  if(is.list(s)) lapply(s,bin) else bin(s)
}




##' Manipulate the samples generated by the Metropolis samplers.
##'
##' These functions provide some basic operations on the samples
##' generated by the Metropolis samplers for Estelle and Stella.
##'
##' @rdname chain.summary
##' @title Manipulate MCMC samples
##' @param s a single chain or a list of parallel chains generated by
##' \code{estelle.metropolis} or \code{stella.metropolis}.
##' @param discard number of initial samples to discard.
##' @param thin rate at which to thin the sample.
##' @param collapse whether to collapse parallel chains to a single chain
##' @param chains the set of chains to retain, or \code{NULL}.
##' @return
##' \code{chain.summary} returns a summary of the sample
##' \code{chain.tail} discards the initial samples from each chain
##' \code{chain.last} returns the last sample for each location in each chain
##' \code{chain.collapse} collapses multiple chains into a single sample
##' \code{chain.cov} returns the covariance of the parameters location by location as a pxpxn array.
##' \code{chain.bcov} returns the joint covariance of the parameters as an (np)x(np) array.
##' \code{chain.acceptance} returns the acceptance rate in the (thinned) chain
##' @export
chain.summary <- function(s) {
  dm <- dim(s[[1]])
  cat("Sample of",dm[3L],"from",length(s),"chains of",dm[2L],"parameters for",dm[1L],"locations\n")
}

##' @rdname chain.summary
##' @export
chain.tail <- function(s,discard=0,thin=1) {
  tail <- function(s) s[,,seq.int(from=1+max(discard,0),to=dim(s)[3L],by=thin)]
  if(!is.list(s)) tail(s) else lapply(s,tail)
}

##' @rdname chain.summary
##' @export
chain.last <- function(s) {
  last <- function(s) s[,,dim(s)[3L]]
  if(!is.list(s)) last(s) else lapply(s,last)
}


##' @rdname chain.summary
##' @export
chain.collapse <- function(s,collapse=TRUE,discard=0,thin=1,chains=NULL) {
  subset <- function(s) s[,,seq.int(from=1+max(discard,0),to=dim(s)[3L],by=thin)]
  if(!is.list(s)) {
    if(thin>1 || discard>0)
      s <- subset(s)
  } else {
    if(!is.null(chains)) s <- s[chains]
    if(thin>1 || discard>0) s <- lapply(s,subset)
    if(collapse) {
      dm <- dim(s[[1]])
      s <- array(unlist(s),c(dm[1:2],length(s)*dm[3]))
    }
  }
  s
}



##' @rdname chain.summary
##' @export
chain.cov <- function(s,discard=0,chains=NULL) {
  s <- chain.collapse(s,collapse=FALSE,discard=discard,chains=chains)

  if(!is.list(s)) {
    V <- apply(s,1L,function(x) var(t(x)))
  } else {
    dm <- dim(s[[1]])
    V <- apply(array(unlist(lapply(s, function(s) apply(s,1L,function(y) var(t(y))))),
                     c(dm[c(2L,2L,1L)],length(s))),
               1:3,mean)
  }
  V
}



##' @rdname chain.summary
##' @export
chain.bcov <- function(s,discard=0,chains=NULL) {
  bcov <- function(s) {
    dm <- dim(s)
    dim(s) <- c(prod(dm[1:2]),dm[3])
    var(t(s))
  }

  s <- chain.collapse(s,collapse=FALSE,discard=discard,chains=chains)
  if(is.list(s))
    apply(simplify2array(lapply(s,bcov)),1:2,mean)
  else
    bcov(s)
}

##' @rdname chain.summary
##' @export
chain.acceptance <- function(s,collapse=FALSE,chains=NULL) {
  rate <- function(s) mean(apply(s,1,function(x) mean(rowMeans(x[,-1L]-x[,-ncol(x)]!=0))))

  s <- chain.collapse(s,collapse=FALSE,chains=chains)
  r <- if(is.list(s)) lapply(s,rate) else rate(s)
  if(collapse & is.list(r)) do.call(mean,r) else r
}


##' Apply a function to chains of samples
##'
##' The function \code{chain.apply} applies a function of a single
##' argument to the samples of a chain or a list of parallel chains
##' generated by \code{estelle.metropolis} or
##' \code{stella.metropolis}.
##'
##' The function \code{chain.apply2} applies a function of two
##' arguments to samples of two chain or two lists of parallel chains
##' generated by \code{estelle.metropolis} or
##' \code{stella.metropolis}.
##'
##' @title Apply a function to chains
##' @param f the function to apply to samples
##' @param s1 a single chain or a list of parallel chains generated by
##' \code{estelle.metropolis} or \code{stella.metropolis}.
##' @param s2 a single chain or a list of parallel chains generated by
##' \code{estelle.metropolis} or \code{stella.metropolis}.
##' @param collapse whether to collapse parallel chains to a single chain.
##' @param thin rate at which to thin the sample.
##' @param chains the set of chains to retain, or \code{NULL}.
##' @export
chain.apply <- function(s1,f,collapse=TRUE,thin=1,chains=NULL) {
  s1 <- chain.collapse(s1,collapse=collapse,thin=thin,chains=chains)
  applyf <- function(s1) sapply(seq_len(dim(s1)[3]),function(k) f(s1[,,k]))
  if(!collapse) lapply(s1,applyf) else applyf(s1)
}


##' @rdname chain.apply
##' @export
chain.apply2 <- function(s1,s2,f,collapse=TRUE,thin=1,chains=NULL) {
  s1 <- chain.collapse(s1,collapse=collapse,thin=thin,chains=chains)
  s2 <- chain.collapse(s2,collapse=collapse,thin=thin,chains=chains)
  applyf <- function(s1,s2) sapply(seq_len(dim(s1)[3]),function(k) f(s1[,,k],s2[,,k]))
  if(!collapse) mapply(applyf,s1,s2,SIMPLIFY=FALSE) else applyf(s1,s2)
}



##' Convert to Coda objects
##'
##' Convert samples generated by \code{estelle.metropolis} or
##' \code{stella.metropolis} to a \pkg{coda} object.
##'
##' @title Export to Coda
##' @param s a list of chains generated by \code{estelle.metropolis}
##' or \code{stella.metropolis}.
##' @return a \pkg{coda} object.
##' @importFrom coda mcmc mcmc.list
##' @export
chain.coda <- function(s) {
  coda <- function(s) {
      dm <- dim(s)
      dim(s) <- c(prod(dm[1:2]),dm[3])
      nms <- c("Lon","Lat")
      if(dm[2]>2)
          nms <- c("Lon","Lat",paste0("P",seq.int(length.out=dm[2]-2)))
      nms <- as.vector(t(outer(nms,1:dm[1],paste,sep=".")))
      rownames(s) <- nms
      mcmc(t(s))
  }
  if(is.list(s)) {
      if(length(s)==1) coda(s[[1]]) else do.call(mcmc.list,lapply(s,coda))
  } else {
      coda(s)
  }
}






##' Calculate shape and rate parameters of the Gamma distribution
##'
##' The Gamma distribution is usually parameterized in terms of the
##' shape and rate parameters.  The \code{parameters.gamma} function
##' deterimines the shape and rate parameters that will yield a Gamma
##' distribution with a desired mean and standard deviation.
##' @title Alternate Gamma Parametrization
##' @param mean vector of means.
##' @param sd vector of standard deviations.
##' @return Returns a 2 column array of the shape and rate parameters
##' that will yield the required mean and standard deviation.
##' @examples
##' ## Shape and rate that give a mean of 3 and sd of 2
##' parameters.gamma(3,2)
##' @export
parameters.gamma <- function(mean,sd) {
  drop(cbind(shape=(mean/sd)^2,rate=mean/sd^2))
}

##' Permute and collapse the dimensions of an array
##'
##' Like \code{aperm} this function permutes the dimensions of an
##' array, but can also collapse adjacent dimensions.  The permutation
##' of dimensions is specified by a list of vectors.  The index of
##' each dimension must appear exactly once, but dimensions that occur
##' together in a vector are reduced to a single dimension.
##' @title Generalized Array Permutation
##' @param x an array.
##' @param perm the subscript permutation list.
##' @return The redimensioned array.
##' @examples
##' x <- array(1:120,2:5)
##' dim(gperm(x,list(4,3,2,1)))
##' dim(gperm(x,list(c(4,3),2,1)))
##' dim(gperm(x,list(1,c(3,2),4)))
##'
##' @export
gperm <- function(x,perm) {
  dm <- dim(x)
  x <- aperm(x,unlist(perm))
  dim(x) <- sapply(perm,function(k) prod(dm[k]))
  x
}




##' Construct a sampler to draw multivariate Normal deviates.
##'
##' Construct a sampler that draws Multivariate Normal deviates with
##' covariances determined by \code{S} and mean determined by its
##' first argument.
##'
##' The \code{mvnorm} function constructs a function that generates
##' \code{n} independent Multivariate Normal deviates, where the
##' covariance of each deviate is specified by \code{S} which must be
##' either a \code{m}x\code{m} covariance matrix or an
##' \code{n}x\code{m}x\code{m} array of covariance matrices.
##'
##' The \code{bmvnorm} constructs a function that generates \code{n}
##' correlated Multivariate Normal deviates, where the joint
##' covariance is specified by \code{S} which must be a
##' \code{nm}x\code{nm} covariance matrix (as generated by
##' \code{chain.bcov}).
##' @title Multivariate Normal Samples
##' @param S a covariance matrix or an array of covariance matrices.
##' @param s a scale factor applied to S.
##' @param n number of deviates to draw.
##' @param m dimension of each deviate.
##' @param tol minimum allowable variance.
##' @return A function that draws bivariate Normal deviates with mean
##' given by its first argument.
##' @export
mvnorm <- function(S,s=1,n=1,tol=1.0E-6) {

  ## Fault tolerant cholesky
  fchol <- function(V) {
    diag(V) <- pmax.int(diag(V),tol)
    tryCatch(chol(V),
             error=function(e) {
               d <- diag(V)
               V <- V/4
               diag(V) <- pmax.int(d/2,tol)
               chol(V)
             })
  }

  m <- dim(S)[1L]
  if(length(dim(S))==2) {
    S <- array(s*fchol(S),c(m,m,n))
  } else {
    n <- dim(S)[3L]
    for(k in 1:n) {
      S[,,k] <- s*fchol(S[,,k])
    }
  }
  S <- aperm(S,c(1L,3L,2L))
  dim(S) <- c(m*n,m)

  function(mu) {
    z <- rnorm(m*n)*S
    dim(z) <- c(m,m*n)
    z <- colSums(z)
    dim(z) <- c(n, m)
    mu+z
  }
}

##' @rdname mvnorm
##' @export
bmvnorm <- function(S,m,s=1) {
  S <- chol(s*S)
  n <- nrow(S)/m

  function(mu) {
    z <- rnorm(m*n)%*%S
    dim(z) <- c(n,m)
    mu+z
  }
}

##' Convert GeoLight data
##'
##' This function converts from the tFirst, tSecond format used by
##' GeoLight to the twilight, rise format used by Stella and Estelle.
##' @title Convert GeoLight Format
##' @param tFirst times of first twilight.
##' @param tSecond times of second twilight.
##' @param type type of twilight.
##' @return A data frame with columns \item{\code{twilight}}{times of
##' twilight as POSIXct objects.}  \item{\code{rise}}{logical vector
##' indicating which twilights are sunrise.}
##' @export
geolight.convert <- function(tFirst,tSecond,type) {
  tm <- .POSIXct(c(as.POSIXct(tFirst,"GMT"),
                   as.POSIXct(tSecond,"GMT")),"GMT")
  keep <- !duplicated(tm)
  tm <- tm[keep]
  rise <- c(type==1,type!=1)[keep]
  ord <- order(tm)
  data.frame(Twilight=tm[ord],Rise=rise[ord])
}


##' Twilight times from a flight of a Short-tailed Shearwater.
##'
##' Twilight times derived from light intensity measurements from an
##' archival tag on a Short-tailed Shearwater (\emph{Puffinus
##' tenuirostris}).  The bird was tagged at its burrow on Wedge
##' Island, Tasmania Australia (147.670E 43.133S).  The bird makes two
##' foraging trips, returning to its burrow after each trip.
##'
##' When the bird is in its burrow, no light is recorded by the tag
##' and the twilights reported in the data set are the calculated
##' twilights Wedge Island when solar zenith is 95 degrees.
##'
##' These data supplied courtesy of Delphi Ward and Mark Hindell,
##' Institute of Marine and Antarctic Studies, University of
##' Tasmania.
##' @name Shearwater
##' @docType data
##' @title Short-tailed Shearwater Track
##' @format A data frame with 3 columns and 72 rows.  The columns
##' represent
##' \tabular{rl}{
##' \code{twilight} \tab times of twilight. \cr
##' \code{rise} \tab logical indicating which twilights are sunrise. \cr
##' \code{burrow} \tab logical indicating when the bird in its burrow.
##' }
##' @keywords data
NULL

##' Light data from the foraging trips of the Southern elephant seal.
##'
##' Light intensity measurements over time from an archival tag on a
##' Southern elephant seal (\emph{Mirounga leonina}).  The seal was
##' tagged at the isthmus on Macquarie Island (158.950E, 54.5S), with
##' data from a time-depth-recorder (Mk9 TDR; Wildlife Computers,
##' Seattle, WA, USA).  These tags provides regular time series of
##' measurements of depth, water temperature, and ambient light
##' level. The original data for \code{ElephantSeal1} were processed
##' to remove values at depths greater than 15m and to classify
##' periods of twilight. The data for \code{ElephantSeal2} have also
##' been processed for twilight periods. The seals makes one single
##' foraging trip, returning to the isthmus where they were tagged.
##' Data recorded while the seal is at the isthmus are used for
##' calibration (\code{\link{ElephantSeal1calib}}).
##'
##' These data supplied courtesy of Mark Hindell, Institute of Marine
##' and Antarctic Studies, University of
##' Tasmania. \code{ElephantSeal1} is B362_99 and \code{ElephantSeal2}
##' is C699_02
##' @name ElephantSeal1
##' @aliases ElephantSeal1calib ElephantSeal2 ElephantSeal2calib
##' @docType data
##' @title Southern Elephant seal tag data
##' @format \code{ElephantSeal1} A data frame with 3 columns.  The
##' columns represent
##' \tabular{rl}{
##' \code{time} \tab times of measurement \cr
##' \code{light} \tab  (log) light values \cr
##' \code{segment} \tab integer indicating sequence of twilight periods \cr
##' }
##' \code{ElephantSeal2} This tag is similar to \code{ElephantSeal1}
##' but also has columns
##' \tabular{rl}{
##' \code{depth} \tab depth in the water column in metres (positive) \cr
##' \code{temp} \tab temperature in the water column in degrees celcius \cr
##' }
##' \code{ElephantSeal1calib} and \code{ElephantSealcalib2} A data
##' frame with 2 columns.  The columns represent
##' \tabular{rl}{
##' \code{zenith} \tab zenith values \cr
##' \code{light} \tab  (log) light values \cr
##' }
##' @keywords data
NULL

