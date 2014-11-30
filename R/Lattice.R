cellFromLonLat <- function(raster,ps) {
  if(!(is.na(projection(raster)) || isLonLat(raster))) {
    ps <- SpatialPoints(ps,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
    ps <- coordinates(spTransform(ps,projection(raster,FALSE)))
  }
  cellFromXY(raster,ps)
}


essie.forback <- function(model,raster,epsilon=1.0E-4) {
  n <- model$n
  cs <- 1:ncell(raster)
  pts <- longlatFromCell(raster,cs)
  fixed.cells <- cellFromLonLat(raster,model$fixed.locations)

  scale <- function(x) x/sum(x)

  lattice <- lapply(1:n,function(k) {
    if(fixed[k]!=0) {
      list(cell=fixed.cells[fixed[k]],fs=1)
    } else {
      ps <- exp(model$logpk(k,pts))
      keep <- ps > epsilon*max(ps)
      list(cs=cs[keep],ps=ps[keep])
    }
  })

  ## Forward iteration
  as0 <- scale(lattice[[1]]$ps)
  lattice[[1]]$as <- as0
  xs <- pts[lattice[[1]]$cs,,drop=F]
  for(k in 2:n) {
    xs0 <- xs
    xs <- pts[lattice[[k]]$cs,,drop=F]
    as <- 0
    for(i in 1:nrow(xs0))
      as <- as + as0[i]*exp(logbk(k-1,xs0[i,],xs))
    as0 <- scale(as*lattice[[k]]$ps)
    lattice[[k]]$as <- as0
  }

  ## Backward iteration
  bs1 <- lattice[[n]]$ps
  xs1 <- pts[lattice[[n]]$cs,,drop=F]

  ## Store trimmed grid
  cs <- lattice[[n]]$cs
  ps <- lattice[[n]]$as
  keep <- ps > epsilon*max(ps)
  lattice[[n]] <- list(cell=cs[keep],p=scale(ps[keep]))

  for(k in (n-1):1) {
    xs <- pts[lattice[[k]]$cs,,drop=F]
    bs <- 0
    for(i in 1:nrow(xs1))
      bs <- bs + bs1[i]*exp(logbk(k,xs,xs1[i,]))
    bs <- scale(bs)
    bs1 <- bs*lattice[[k]]$ps
    xs1 <- xs

    ## Store trimmed grid
    cs <- lattice[[k]]$cs
    ps <- bs*lattice[[k]]$as
    keep <- ps > epsilon*max(ps)
    lattice[[k]] <- list(cell=cs[keep],p=scale(ps[keep]))
  }
  list(raster=raster,lattice=lattice)
}

