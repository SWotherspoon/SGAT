##' Calculate great circle distancces between points.
##'
##' Compute the great circle distances from x1 to x2.
##' @title Great Circle Distance
##' @param x1 a two column matrix of (lon,lat) locations.
##' @param x2 a two column matrix of (lon,lat) locations.
##' @return vector of interpoint distances (km)
##' @export
gcDist <- function(x1,x2) {
  rad <- pi/180
  s1 <- sin(rad*x1[,2L])
  s2 <- sin(rad*x2[,2L])
  6378.137*acos(pmin.int(sqrt(1-s1*s1)*sqrt(1-s2*s2)*cos(rad*(x2[,1L]-x1[,1L]))+s1*s2,1))
}

## Alternate code
gcDist0 <- function(x1,x2) {
  rad <- pi/180
  6378.137*acos(pmin.int(
    cos(rad*x1[,2L])*cos(rad*x2[,2L])*cos(rad*(x2[,1L]-x1[,1L]))+
      sin(rad*x1[,2L])*sin(rad*x2[,2L]),1))
}


##' Calculate great circle distancces between all pairs of points.
##'
##' Compute the great circle distances from every point in x1 to every
##' point is x2, and return the result as a matrix.
##' @title Great Circle Distance
##' @param x1 a two column matrix of (lon,lat) locations.
##' @param x2 a two column matrix of (lon,lat) locations.
##' @return vector of interpoint distances (km)
##' @export
gcOuterDist <-function(x1,x2) {
  rad <- pi/180
  s1 <- sin(rad*x1[,2L])
  s2 <- sin(rad*x2[,2L])
  matrix(6378.137*acos(pmin.int(
    (sqrt(1-s1*s1)%o%sqrt(1-s2*s2))*cos(outer(rad*x1[,1L],rad*x2[,1L],"-"))+s1%o%s2,1)),
         length(x1),length(x2))
}







##' Threshold Model Structures for Essie
##'
##' Essie requires a model structure that describes the model being
##' fitted. This function generates basic model structures for
##' threshold twilight data that should provide a suitable starting
##' point for most analyses.
##'
##' The \code{essieThresholdModel} function constructs a model structure
##' assuming that each twilight time is associated with a single
##' location.
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
##' \code{alpha}, which must be a vector of parameters that are
##' to be applied to each twilight.
##'
##' The \code{twilight.model} argument selects the distribution of the
##' twilight errors
##' \describe{
##' \item{'Normal'}{Normally distributed with mean \code{alpha[,1]} and
##' standard deviation \code{alpha[,2]},}
##' \item{'LogNormal'}{Log Normally distributed so the log errors have
##' mean \code{alpha[,1]} and standard deviation \code{alpha[,2]}, or}
##' \item{'Gamma'}{Gamma distributed with shape \code{alpha[,1]} and
##' rate \code{alpha[,2]}.}
##' }
##'
##' The initialization locations \code{x0} are only required when
##' specifying fixed locations.
##'
##' Essie assumes that the average speed of travel between successive
##' locations is Gamma distributed. By default, the speed of travel is
##' calculated based on the time intervals between the twilights (in
##' hours), but the intervals of time actually available for travel
##' can be specified directly with the \code{dt} argument. The
##' parameters \code{beta[1]} and \code{beta[2]} specify the shape and
##' rate of the Gamma distribution of speeds.
##'
##' Twilights can be missing because either the light record was too
##' noisy at that time to estimate twilight reliably, or because the
##' tag was at very high latitude and no twilight was observed.
##' Missing twilights should be replaced with an approximate time of
##' twilight, and the vector \code{missing} used to indicate which
##' twilights are approximate and which are true.  This should be a
##' vector of integers, one for each twilight where the integer codes
##' signify
##' \describe{
##' \item{0:}{The twilight is not missing.}
##' \item{1:}{The twilight is missing, but a twilight did occur.}
##' \item{2:}{The twilight is missing because twilight did not occur.}
##' \item{3:}{The twilight is missing and it is not known if a twilight occurred.}
##' }
##'
##' @title Threshold Model Structures (Essie)
##' @param twilight the observed times of twilight as POSIXct.
##' @param rise logical vector indicating which twilights are sunrise.
##' @param twilight.model the model for the errors in twilight times.
##' @param alpha parameters of the twilight model.
##' @param beta parameters of the behavioural model.
##' @param logp0 function to evaluate any additional contribution to
##' the log posterior from the twilight locations.
##' @param x0 suggested starting points for twilight locations.
##' @param fixed logical vector indicating which twilight locations
##' to hold fixed.
##' @param missing integer vector indicating which twilights are missing.
##' @param dt time intervals for speed calculation in hours.
##' @param zenith the solar zenith angle that defines twilight.
##' @return a list with components
##' \item{\code{logpk}}{function to evaluate the contributions to the
##' log posterior from the k-th twilight}
##' \item{\code{logpbk}}{function to evaluate contribution to
##' the log posterior from the behavioural model for the k-th track segment.}
##' \item{\code{fixed}}{a logical vector indicating which locations
##' should remain fixed.}
##' \item{\code{x0}}{an array of initial twilight locations.}
##' \item{\code{time}}{the twilight times.}
##' \item{\code{rise}}{the sunrise indicators.}
##' \item{\code{alpha}}{the twilight model parameters.}
##' \item{\code{beta}}{the behavioural model parameters.}
##' @importFrom stats dgamma dnorm dlnorm
##' @export
essieThresholdModel <- function(twilight,rise,
                                  twilight.model=c("LogNormal","Gamma","Normal"),
                                  alpha,beta,
                                  logp0=function(k,x) 0,
                                  x0,fixed=FALSE,missing=0,dt=NULL,zenith=96) {

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(twilight)/3600)

  ## Fixed locations
  fixed <- rep_len(fixed,length.out=length(twilight))

  ## Extend missing
  missing <- rep_len(missing,length.out=length(twilight))

  ## Twilight residuals
  residuals <- function(k,x) {
    sgn <- ifelse(rise[k],1,-1)
    s <- solar(twilight[k])
    4*sgn*(s$solarTime-twilightSolartime(s,x[,1L],x[,2L],rise[k],zenith))
  }


  ## Select the density of the twilight residuals
  twilight.model <- match.arg(twilight.model)
  logp.residual <- switch(twilight.model,
                          Gamma=function(r) dgamma(r,alpha[1L],alpha[2L],log=TRUE),
                          LogNormal=function(r) dlnorm(r,alpha[1L],alpha[2L],log=TRUE),
                          Normal=function(r) dnorm(r,alpha[1L],alpha[2L],log=TRUE))

  ## Contribution to log posterior from each x location
  logpk <- function(k,x) {
    if(missing[k]<3) {
      logp <- logp.residual(residuals(k,x))
      if(missing[k]==0) {
        logp[!is.finite(logp)] <- -Inf
      } else {
        if(missing[k]==1)
          logp <- ifelse(is.finite(logp),0,-Inf)
        if(missing[k]==2)
          logp <- ifelse(is.finite(logp),-Inf,0)
      }
      logp+logp0(k,x)
    } else {
      logp0(k,x)
    }
  }

  ## Behavioural contribution to the log posterior
  logbk <- function(k,x1,x2) {
    spd <- pmax.int(gcDist(x1,x2), 1e-06)/dt[k]
    dgamma(spd,beta[1L],beta[2L],log=TRUE)
  }

  list(
    logpk=logpk,
    logbk=logbk,
    fixed=fixed,
    x0=x0,
    time=twilight,
    rise=rise,
    alpha=alpha,
    beta=beta)
}



##' Curve Model Structures for Essie
##'
##' Essie requires a model structure that describes the model being
##' fitted. This function generates basic model structures for curve
##' fitting metthods that should provide a suitable starting point for
##' most analyses.
##'
##' The \code{essieCurveModel} function constructs a model structure
##' assuming that each twilight profile is associated with a single
##' location.  The errors in observed log light level are assumed to
##' be Normally distributed about their expected value.
##'
##' The properties of the twilight model are determined by
##' \code{alpha}, which must be a vector of parameters that are
##' to be applied to each twilight.
##'
##' The initialization locations \code{x0} are only required when
##' specifying fixed locations.
##'
##' Essie assumes that the average speed of travel between successive
##' locations is Gamma distributed. By default, the speed of travel is
##' calculated based on the time intervals between the twilights (in
##' hours), but the intervals of time actually available for travel
##' can be specified directly with the \code{dt} argument. The
##' parameters \code{beta[1]} and \code{beta[2]} specify the shape and
##' rate of the Gamma distribution of speeds.
##'
##' @title Curve Model Structures (Essie)
##' @param time vector of sample times as POSIXct.
##' @param light vector of observed (log) light levels.
##' @param segment vector of integers that assign observations to
##' twilight segments.
##' @param calibration function that maps zenith angles to expected
##' light levels.
##' @param alpha parameters of the twilight model.
##' @param beta parameters of the behavioural model.
##' @param logp0 function to evaluate any additional contribution to
##' the log posterior from the twilight locations.
##' @param x0 suggested starting points for twilight locations.
##' @param fixed logical vector indicating which twilight locations
##' to hold fixed.
##' @param dt time intervals for speed calculation in hours.
##' @return a list with components
##' \item{\code{logpk}}{function to evaluate the contributions to the
##' log posterior from the k-th twilight}
##' \item{\code{logpbk}}{function to evaluate contribution to
##' the log posterior from the behavioural model for the k-th track segment.}
##' \item{\code{fixed}}{a logical vector indicating which locations
##' should remain fixed.}
##' \item{\code{x0}}{an array of initial twilight locations.}
##' \item{\code{time}}{the twilight times.}
##' \item{\code{rise}}{the sunrise indicators.}
##' \item{\code{alpha}}{the twilight model parameters.}
##' \item{\code{beta}}{the behavioural model parameters.}
##' @importFrom stats dgamma dnorm median
##' @export
essieCurveModel <- function(time,light,segment,
                            calibration,alpha,beta,
                            logp0=function(k,x) 0,
                            x0,fixed=FALSE,dt=NULL) {

  ## Median time in each segment
  tm <- .POSIXct(sapply(split(time,segment),median),"GMT")

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(tm)/3600)

  ## Fixed locations
  fixed <- rep_len(fixed,length.out=length(tm))

  ## Convert to solar time
  sun <- solar(time)

  ## Contribution to log posterior from each x location
  logpk <- function(k,x) {
    ## Extract the (solar) times and light for this segment
    seg <- as.numeric(segment)==k
    sun <- solar(time[seg])
    ls <- light[seg]
    ## Contributions to log posterior from the light
    logp <- sapply(seq_len(nrow(x)),
                   function(i) {
                     fs <- calibration(zenith(sun,x[i,1L],x[i,2L]))
                     off <- mean(ls)-mean(fs)
                     sum(dnorm(ls,fs+off,alpha[1L],log=TRUE))+dnorm(off,0,alpha[2L],log=TRUE)
                   })
    ## Add contribution from prior
    logp <- logp + logp0(k,x)
    logp
  }

  ## Behavioural contribution to the log posterior
  logbk <- function(k,x1,x2) {
    spd <- pmax.int(gcDist(x1,x2), 1e-06)/dt[k]
    dgamma(spd,beta[1L],beta[2L],log=TRUE)
  }

  list(
    logpk=logpk,
    logbk=logbk,
    fixed=fixed,
    x0=x0,
    time=tm,
    alpha=alpha,
    beta=beta)
}



##' Fit an Essie model
##'
##' Essie is a discrete analog of Stella that only considers locations
##' on a lattice of grid points.
##'
##' The posterior probability that the tag is at a particular location
##' is determined by a two pass recursive algorithm.  The forward
##' sweep propogates location information forward in time, the
##' backward sweep propogates location information backward in time,
##' and the full posterior is a compromise of these two.
##'
##' The method is only approximate and is controlled by the two
##' tolerances \code{epsilon1} and \code{epsilon2}.  The first limits
##' the precision of the likelihood - locations for which the
##' likelihood fall below \code{epsilon1} of its maximum are ignored.
##' The second limits the precision of the transition probabilities -
##' locations for which the probability below \code{epsilon2} of the
##' maximum are ignored in the movement calculation.  The smaller
##' these parameters are set, the slower but more accurate the fit.
##'
##' @title Essie
##' @param model a model structure generated by \code{essieThresholdModel}.
##' @param grid a raster object defining the grid
##' @param epsilon1 likelihood tolerance
##' @param epsilon2 transition probability tolerance.
##' @param verbose report progress at prompt?
##' @return a list with elements
##' \item{\code{grid}}{a raster object that defines the grid.}
##' \item{\code{times}}{the times corresponding to the location estimates.}
##' \item{\code{lattice}}{a sparse grid representation of the posterior location probabilities.}
##' @importFrom utils flush.console
##' @importFrom raster ncell values
##' @export
essie <- function(model,grid,epsilon1=1.0E-3,epsilon2=1.0E-16,verbose=interactive()) {

  normalize <- function(x) {
    s <- sum(x)
    if(s==0) x else x/sum(x)
  }

  n <- length(model$time)
  pts <- lonlatFromCell(grid,seq_len(ncell(grid)))
  fixed <- integer(n)
  fixed[model$fixed] <- cellFromLonLat(grid,model$x0[model$fixed,])

  ## Compute likelihood
  cs <- seq_len(ncell(grid))[values(grid)!=0]
  lattice <- lapply(seq_len(n),function(k) {
    if(fixed[k]!=0) {
      list(cs=fixed[k],ps=1)
    } else {
      logps <- model$logpk(k,pts[cs,,drop=FALSE])
      keep <- logps > log(epsilon1)+max(logps)
      list(cs=cs[keep],ps=normalize(exp(logps[keep])))
    }
  })

  ## Forward iteration
  if(verbose) {
    cat("Fwd ",sprintf("%6d",1))
    flush.console()
  }
  xs <- pts[lattice[[1]]$cs,,drop=FALSE]
  as <- lattice[[1]]$as <- lattice[[1]]$ps
  for(k in 2:n) {
    if(verbose) {
      cat("\b\b\b\b\b\b")
      cat(sprintf("%6d",k))
      flush.console()
    }
    xs0 <- xs
    xs <- pts[lattice[[k]]$cs,,drop=FALSE]
    as0 <- as
    as <- 0
    for(i in which(as0 > epsilon2*max(as0)))
      as <- as + as0[i]*exp(model$logbk(k-1,xs0[i,,drop=FALSE],xs))
    as <- normalize(as*lattice[[k]]$ps)
    lattice[[k]]$as <- as
  }


  ## Backward iteration
  if(verbose) {
    cat("\nBwd ",sprintf("%6d",1))
    flush.console()
  }
  xs <- pts[lattice[[n]]$cs,,drop=FALSE]
  bs <- lattice[[n]]$bs <- lattice[[n]]$ps
  for(k in (n-1):1) {
    if(verbose) {
      cat("\b\b\b\b\b\b")
      cat(sprintf("%6d",k))
      flush.console()
    }
    xs0 <- xs
    xs <- pts[lattice[[k]]$cs,,drop=FALSE]
    bs0 <- bs
    bs <- 0
    for(i in which(bs0 > epsilon2*max(bs0)))
      bs <- bs + bs0[i]*exp(model$logbk(k,xs,xs0[i,,drop=FALSE]))
    bs <- normalize(bs*lattice[[k]]$ps)
    lattice[[k]]$bs <- bs
  }

  list(grid=grid,time=model$time,lattice=lattice)
}



##' Return the probabilities from an Essie model as a raster object.
##'
##' This function returns the likelihood or the posterior probability
##' of a single location as a raster object.
##'
##' @title Essie probability rasters
##' @param obj a fitted object generated by \code{essie}
##' @param k the index of the location
##' @param type whether to return the full estimate, the estimate from
##' the forward or backward pass of the algorithm or the likelihood.
##' @return a raster of probabilities/likelihoods
##' @importFrom raster raster
##' @export
essieRaster <- function(obj,k,type=c("full","forward","backward","likelihood")) {
  g <- raster(obj$grid)
  names(g) <- obj$time[k]
  l <- obj$lattice[[k]]
  type <- match.arg(type)
  g[l$cs] <- switch(type,
                    full=l$as*l$bs/l$ps,
                    forward=l$as,
                    backward=l$bs,
                    likelihood=l$ps)
  g
}





##' Compute posterior mean or mode locations from an Essie fit
##'
##' These functions compute the posterior mean or posterior mode
##' locations for each time point from and Essie fit. The
##' \code{essieCircularMean} function assumes the computational grid
##' is in (lon,lat) projection and computes a circular mean for
##' longitude.
##'
##' @title Essie locations
##' @param obj a fitted object generated by \code{essie}
##' @param type whether to return the full estimate or just the
##' estimate from the forward or backward pass of the algorithm.
##' @return a list with elements
##' \item{\code{times}}{the times corresponding to the location estimates.}
##' \item{\code{x}}{a two column matrix of positions.}
##' @export
essieMean <- function(obj,type=c("full","forward","backward")) {
  type <- match.arg(type)
  pts <- lonlatFromCell(obj$grid,seq_len(ncell(obj$grid)))
  x <- do.call(rbind,lapply(obj$lattice,function(l) {
      ps <- switch(type,
                   full=l$as*l$bs/l$ps,
                   forward=l$as,
                   backward=l$bs)
      colSums(ps*pts[l$cs,,drop=FALSE])/sum(ps)
    }))
  colnames(x) <- c("lon","lat")
  list(time=obj$time,x=x)
}


##' @rdname essieMean
##' @export
essieCircularMean <- function(obj,type=c("full","forward","backward")) {
  type <- match.arg(type)
  pts <- lonlatFromCell(obj$grid,seq_len(ncell(obj$grid)))
  ## Convert longitude
  pts <- cbind(sin(pts[,1]*pi/180),cos(pts[,1]*pi/180),pts[,2])
  x <- do.call(rbind,lapply(obj$lattice,function(l) {
      ps <- switch(type,
                   full=l$as*l$bs/l$ps,
                   forward=l$as,
                   backward=l$bs)
      mn <- colSums(ps*pts[l$cs,,drop=FALSE])/sum(ps)
      c(180/pi*atan2(mn[1],mn[2]),mn[3])
    }))
  colnames(x) <- c("lon","lat")
  list(time=obj$time,x=x)
}


##' @rdname essieMean
##' @export
essieMode <- function(obj,type=c("full","forward","backward")) {
  type <- match.arg(type)
  cs <- unlist(lapply(obj$lattice,function(l) {
    ps <- switch(type,
                 full=l$as*l$bs/l$ps,
                 forward=l$as,
                 backward=l$bs)
    l$cs[which.max(ps)]
  }))
  x <- lonlatFromCell(obj$grid,cs)
  colnames(x) <- c("lon","lat")
  list(time=obj$time,x=x)
}
