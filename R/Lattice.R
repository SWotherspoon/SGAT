library(SGAT)
library(raster)
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
                                  alpha,beta,
                                  logpk0=function(k,x) 0,
                                  x0,fixed=FALSE,dt=NULL,zenith=96) {

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(twilight)/3600)

  fixed <- rep_len(fixed,length.out=length(twilight))

  logpk <- function(k,x) {
    sgn <- ifelse(rise[k],1,-1)
    s <- solar(twilight[k])
    r <- 4*sgn*(s$solarTime-twilight.solartime(s,x[,1L],x[,2L],rise[k],zenith))
    ifelse(!is.finite(r) | r < 0, -Inf,dlnorm(r,alpha[1L],alpha[2L],log=TRUE)+logpk0(k,x))
  }

  logbk <- function(k,x1,x2) {
    spd <- pmax.int(gcDist(x1,x2), 1e-06)/dt[k]
    dgamma(spd,beta[1L],beta[2L],log=TRUE)
  }


  list(
    n=length(twilight),
    logpk=logpk,
    logbk=logbk,
    x0=x0,
    fixed=fixed)
}

zenith <- 96.2
threshold <- 3
twl <- read.csv(paste0(19745,"twl.csv"),header=T)
twl$Twilight <- as.POSIXct(twl$Twilight,"GMT")
x0 <- matrix(c(147.67,-43.133),nrow(twl),2,byrow=T)
fixed <- twl$Marker > 0


alpha <- c(2.2,1.0)
beta <- c(6, 0.12)
model <- essie.threshold.model(twl$Twilight,twl$Rise,alpha=alpha,beta=beta,x0=x0,fixed=fixed,zenith=zenith)
#grid <- raster(ncols=1*(227-89),nrows=1*(66+75),xmn=89,xmx=227,ymn=-75,ymx=66)

grid <- crop(raster("C:/Reynolds/lsmask.nc"),extent(89,227,-75,66))
plot(grid)

essie.forback <- function(model,grid,epsilon=1.0E-10) {
  n <- model$n
  pts <- lonlatFromCell(grid,1:ncell(grid))
  fixed <- integer(n)
  fixed[model$fixed] <- cellFromLonLat(grid,model$x0[model$fixed,])

  normalize <- function(x) x/sum(x)

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
  xs <- pts[lattice[[1]]$cs,,drop=F]
  as <- lattice[[1]]$as <- lattice[[1]]$ps
  for(k in 2:n) {
    print(k)
    xs0 <- xs
    as0 <- as
    xs <- pts[lattice[[k]]$cs,,drop=F]
    as <- 0
    for(i in which(as0 > epsilon*max(as0)))
      as <- as + as0[i]*normalize(exp(model$logbk(k-1,xs0[i,,drop=F],xs)))
    lattice[[k]]$as <- normalize(as*lattice[[k]]$ps)
  }


  ## Backward iteration
  xs <- pts[lattice[[n]]$cs,,drop=F]
  bs <- lattice[[n]]$bs <- lattice[[n]]$ps
  for(k in (n-1):1) {
    print(k)
    xs0 <- xs
    bs0 <- bs
    xs <- pts[lattice[[k]]$cs,,drop=F]
    bs <- 0
    for(i in which(bs0 > epsilon*max(bs0)))
      bs <- bs + bs0[i]*normalize(exp(model$logbk(k-1,xs,xs0[i,,drop=F])))
    lattice[[k]]$bs <- normalize(bs*lattice[[k]]$ps)
  }


}

plot1 <- function(lattice,grid,k=NULL) {
  opar <- par(mfrow=c(2,2))
  ks <- unlist(lapply(lattice,function(l) l$cs[which.max(l$as)]))
  plot(grid)
  lines(pts[ks,])
  if(!is.null(k)) points(pts[ks[k],,drop=F],pch=16,col="red")
  ks <- unlist(lapply(lattice,function(l) l$cs[which.max(l$bs)]))
  plot(grid)
  lines(pts[ks,])
  if(!is.null(k)) points(pts[ks[k],,drop=F],pch=16,col="red")
  ks <- unlist(lapply(lattice,function(l) l$cs[which.max(l$as*l$bs/l$ps)]))
  plot(grid)
  lines(pts[ks,])
  if(!is.null(k)) points(pts[ks[k],,drop=F],pch=16,col="red")
  par(opar)
}

plot2 <- function(lattice,grid,k) {
  opar <- par(mfrow=c(2,2))
  g <- raster(grid)
  g[lattice[[k]]$cs] <- with(lattice[[k]],normalize(ps))
  plot(g)
  g <- raster(grid)
  g[lattice[[k]]$cs] <- with(lattice[[k]],normalize(as))
  plot(g)
  g <- raster(grid)
  g[lattice[[k]]$cs] <- with(lattice[[k]],normalize(bs))
  plot(g)
  g <- raster(grid)
  g[lattice[[k]]$cs] <- with(lattice[[k]],normalize(as*bs/ps))
  plot(g)
  par(opar)
}
