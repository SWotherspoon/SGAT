cellFromLonLat <- function(raster,ps) {
  if(!(is.na(projection(raster)) || isLonLat(raster))) {
    ps <- SpatialPoints(ps,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
    ps <- coordinates(spTransform(ps,projection(raster,FALSE)))
  }
  cellFromXY(raster,ps)
}



gcDist <- function(x1,x2) {
  rad <- pi/180
  6378.137*acos(pmin.int(
    cos(rad*x1[,2L])*cos(rad*x2[,2L])*cos(rad*(x2[,1L]-x1[,1L]))+
      sin(rad*x1[,2L])*sin(rad*x2[,2L]),1))
}



essie.threshold.model <- function(twilight,rise,
                                  twilight.model=c("LogNormal","Gamma","Normal"),
                                  alpha,beta,
                                  logp0=function(k,x) 0,
                                  x0,fixed=FALSE,dt=NULL,zenith=96) {

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(twilight)/3600)

  ## Fixed locations
  fixed <- rep_len(fixed,length.out=length(twilight))

  ## Twilight residuals
  residuals <- function(k,x) {
    sgn <- ifelse(rise[k],1,-1)
    s <- solar(twilight[k])
    4*sgn*(s$solarTime-twilight.solartime(s,x[,1L],x[,2L],rise[k],zenith))
  }


  ## Select the density of the twilight residuals
  twilight.model <- match.arg(twilight.model)
  logp.residual <- switch(twilight.model,
                          Gamma=function(r) dgamma(r,alpha[1L],alpha[2L],log=TRUE),
                          LogNormal=function(r) dlnorm(r,alpha[1L],alpha[2L],log=TRUE),
                          Normal=function(r) dnorm(r,alpha[1L],alpha[2L],log=TRUE))

  ## Contribution to log posterior from each x location
  logpk <- function(k,x) {
    logp <- logp.residual(residuals(k,x))
    logp[!is.finite(logp)] <- -Inf
    logp <- logp+logp0(k,x)
    logp
  }

  ## Behavioural contribution to the log posterior
  logbk <- function(k,x1,x2) {
    spd <- pmax.int(gcDist(x1,x2), 1e-06)/dt[k]
    dgamma(spd,beta[1L],beta[2L],log=TRUE)
  }

  list(
    time=twilight,
    rise=rise,
    n=length(twilight),
    alpha=alpha,
    beta=beta,
    logpk=logpk,
    logbk=logbk,
    x0=x0,
    fixed=fixed)
}



essie <- function(model,grid,epsilon=1.0E-2,verbose=interactive()) {


  normalize <- function(x) {
    s <- sum(x)
    if(s < 1.0E-12) rep(0,length(x)) else x/sum(x)
  }

  n <- model$n
  pts <- lonlatFromCell(grid,1:ncell(grid))
  fixed <- integer(n)
  fixed[model$fixed] <- cellFromLonLat(grid,model$x0[model$fixed,])

  ## Compute likelihood
  cs <- (1:ncell(grid))[values(grid)>0]
  lattice <- lapply(1:n,function(k) {
    if(fixed[k]!=0) {
      list(cs=fixed[k],ps=1)
    } else {
      logps <- model$logpk(k,pts[cs,,drop=FALSE])
      keep <- logps > log(epsilon)+max(logps)
      list(cs=cs[keep],ps=normalize(exp(logps[keep])))
    }
  })

  ## Forward iteration
  if(verbose) {
    cat("Fwd ",sprintf("%6d",1))
    flush.console()
  }
  xs <- pts[lattice[[1]]$cs,,drop=F]
  as <- lattice[[1]]$as <- lattice[[1]]$ps
  for(k in 2:n) {
    if(verbose) {
      cat("\b\b\b\b\b\b");
      cat(sprintf("%6d",k));
      flush.console()
    }

    xs0 <- xs
    xs <- pts[lattice[[k]]$cs,,drop=F]
    as0 <- as
    as <- 0
    for(i in which(as0>0))
      as <- as + as0[i]*exp(model$logbk(k-1,xs0[i,,drop=F],xs))
    as <- normalize(as*lattice[[k]]$ps)
    lattice[[k]]$as <- as
  }


  ## Backward iteration
  if(verbose) {
    cat("Bwd ",sprintf("%6d",1))
    flush.console()
  }
  xs <- pts[lattice[[n]]$cs,,drop=F]
  bs <- lattice[[n]]$bs <- lattice[[n]]$ps
  for(k in (n-1):1) {
    if(verbose) {
      cat("\b\b\b\b\b\b");
      cat(sprintf("%6d",k));
      flush.console()
    }

    xs0 <- xs
    xs <- pts[lattice[[k]]$cs,,drop=F]
    bs0 <- bs
    bs <- 0
    for(i in which(bs0>0))
      bs <- bs + bs0[i]*exp(model$logbk(k,xs,xs0[i,,drop=F]))
    bs <- normalize(bs*lattice[[k]]$ps)
    lattice[[k]]$bs <- bs
  }

  list(grid=grid,time=model$time,lattice=lattice)
}


lattice.raster <- function(obj,k,type=c("posterior","forward","backward","likelihood")) {
  g <- raster(obj$grid)
  l <- obj$lattice[[k]]
  type <- match.arg(type)
  g[l$cs] <- switch(type,
                    posterior=l$as*l$bs/l$ps,
                    forward=l$as,
                    backward=l$bs,
                    likelihood=l$ps)
  g
}


lattice.mean <- function(obj,type=c("posterior","forward","backward")) {
  type <- match.arg(type)
  cs <- unlist(lapply(obj$lattice,function(l) {
    ps <- switch(type,
                 posterior=l$as*l$bs/l$ps,
                 forward=l$as,
                 backward=l$bs)
    l$cs[which.max(ps)]
  }))
  lonlatFromCell(obj$grid,cs)
}


lattice.maxp <- function(obj,type=c("posterior","forward","backward")) {
  type <- match.arg(type)
  pts <- lonlatFromCell(obj$grid,1:ncell(obj$grid))
  do.call(rbind,lapply(obj$lattice,function(l) {
    ps <- switch(type,
                 posterior=l$as*l$bs/l$ps,
                 forward=l$as,
                 backward=l$bs)
    colSums(ps*pts[l$cs,,drop=F])/sum(ps)
  }))
}

