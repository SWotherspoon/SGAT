  ## Calculate dog-leg distances along an x-z track
    trkdist.xz <- function(x,z) {
    n <- nrow(x)
    cosx2 <- cos(pi/180*x[,2])
    sinx2 <- sin(pi/180*x[,2])
    cosz2 <- cos(pi/180*z[,2])
    sinz2 <- sin(pi/180*z[,2])

    6378.137*(acos(pmin.int(cosx2[-n]*cosz2*cos(pi/180*(z[,1]-x[-n,1]))+sinx2[-n]*sinz2,1))+
              acos(pmin.int(cosx2[-1]*cosz2*cos(pi/180*(z[,1]-x[-1,1]))+sinx2[-1]*sinz2,1)))
  }

  ## Calculate distances along an x track
  trkdist.x <- function(x) {
    n <- nrow(x)
    cosx2 <- cos(pi/180*x[,2])
    sinx2 <- sin(pi/180*x[,2])

    6378.137*acos(pmin.int(cosx2[-n]*cosx2[-1]*cos(pi/180*(x[-1,1]-x[-n,1]))+sinx2[-n]*sinx2[-1],1))
  }



curve.model <- function(datetime, light, segments, calibration,
                        twilight.model = c("LogNormal"),
                        alpha, beta,
                        b.model = c("Gamma"),

                        logp.x = function(x) rep.int(0L, nrow(x)),
                        logp.z = function(z) rep.int(0L, nrow(z)),
                        x0, z0 = NULL,
                        fixedx = FALSE) {
    sun <- solar(datetime)
    tm <- tapply(datetime, segments, median)
    dt <- diff(as.numeric(tm)/3600)
    logpx <-
    switch(twilight.model,
           LogNormal= function(x) {
               elev <- zenith(sun, x[segments, 1], x[segments, 2])
                att <- calibration(elev) ##+ x[, 3]
                logp <- dnorm(light, att, alpha[1], log = TRUE)


               sapply(split(logp, segments), sum) ##+ dnorm(x[,3], 0, alpha[2], log = TRUE)
           }
       )
    ## Contribution to log posterior from each z location
    logpz <- function(z) {
        logp.z(z)
    }
    ## Contribution to log posterior from the movement
    estelle.logpb <- function(x,z) {
        spd <- pmax.int(trkdist.xz(x,z), 1e-06)/dt
        dgamma(spd,beta[1],beta[2],log=TRUE)
    }
    stella.logpb <- function(x) {
        spd <- pmax.int(trkdist.x(x), 1e-06)/dt
        dgamma(spd,beta[1],beta[2],log=TRUE)
    }

    list(
        ## Positional contribution to the log posterior
         logpx=logpx,
         logpz=logpz,
        ## Behavioural contribution to the log posterior
        estelle.logpb = estelle.logpb,
        stella.logpb=stella.logpb,
         ## Locations to be held fixed
       fixedx=fixedx,
       ## Suggested starting points
       x0=x0,
       z0=z0)

}


    ## alpha <- b.mean^2/b.sd^2
    ## beta <- b.mean/b.sd^2
    ## log.sigma <- sqrt(log(1 + bmean^2/b.sd^2))
    ## log.mu <- log(b.mean) - log.sigma^2/2

    ## logpb <-
    ##     switch(behav.model,
    ##            Gamma = function(k, x1, z, x2) {
    ##                spd <- pmax.int(dist(x1, z) + dist(z, x2), 1e-06)/dt[k]
    ##                dgamma(spd, alpha, beta, log = TRUE)
    ##     }
    ## }
    ## else {

    ##     logp.behavioural <- function(k, x1, z, x2) {
    ##         spd <- (dist(x1, z) + dist(z, x2))/dt[k]
    ##         dnorm(spd, log.mu, log.sigma, log = TRUE)
    ##     }
    ## }
