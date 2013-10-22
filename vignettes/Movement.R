### R code from vignette source 'Movement.Rnw'

###################################################
### code chunk number 1: Movement.Rnw:37-38
###################################################
set.seed(32)


###################################################
### code chunk number 2: Movement.Rnw:72-100
###################################################
## Construct a single bezier spline from (t1,p1) to (t2,p2) with
## control points c1 and c2.
bezier <- function(t1,p1,c1,t2,p2,c2) {
  function(t) {
    t <- (unclass(t)-unclass(t1))/(unclass(t2)-unclass(t1))
    outer((1-t)^3,p1)+
      outer(3*(1-t)^2*t,c1)+
        outer(3*(1-t)*t^2,c2)+
          outer(t^3,p2)
  }
}
## Endpoints
p0 <- c(130,-30)
p1 <- c(140,-50)
p2 <- c(130,-40)
## Control points
c0 <- p0+c(6,-4)
c1 <- p1+c(-4,8)
c2 <- p1+c(-6,2)
c3 <- p2+c(2,-8)
f1 <- bezier(30,p0,c0,60,p1,c1)
f2 <- bezier(90,p1,c2,120,p2,c3)
days <- c(1:30,31:60,61:90,91:120,121:150)
ps <- rbind(cbind(rep(p0[1],30),rep(p0[2],30)),
            f1(31:60),
            cbind(rep(p1[1],30),rep(p1[2],30)),
            f2(91:120),
            cbind(rep(p2[1],30),rep(p2[2],30)))


###################################################
### code chunk number 3: Movement.Rnw:104-109
###################################################
library(SGAT)
## Resample to two minute intervals and compute zenith angles
tms <- as.POSIXct("2013-04-01","GMT")+(24*60*60)*days
tms.out <- seq(min(tms),max(tms),by=2*60)
d <- zenith.simulate(tms,ps[,1],ps[,2],tms.out)


###################################################
### code chunk number 4: Movement.Rnw:116-119
###################################################
## Compute times and locations at which twilight is observed.
twl <- twilight.simulate(d)
twl <- twilight.perturb(twl,rlnorm(nrow(twl),0.3,0.6))


###################################################
### code chunk number 5: Movement.Rnw:127-131
###################################################
## Initial x locations
x0 <- threshold.path(twl$Twilight,twl$Rise)$x
## Initial z locations
z0 <- trackMidpts(x0)


###################################################
### code chunk number 6: Movement.Rnw:134-146
###################################################
## Construct model and use residuals to adjust x0
model <- threshold.model(twl$Twilight,twl$Rise,
                         twilight.model="LogNormal",
                         alpha=c(0.3,0.6),beta=c(8,3.5),
                         x0=x0,z0=z0)
r <- model$residuals(x0)
x0 <- cbind(x0[,1]-ifelse(twl$Rise,1,-1)*pmin(0,r-0.1)/1440*360,x0[,2])
z0 <- trackMidpts(x0)
model <- threshold.model(twl$Twilight,twl$Rise,
                         twilight.model="LogNormal",
                         alpha=c(0.3,0.6),beta=c(8,3.5),
                         x0=x0,z0=z0)


###################################################
### code chunk number 7: Movement.Rnw:151-154
###################################################
## Define initial proposals
x.proposal <- mvnorm(S=diag(c(0.002,0.002)),n=nrow(twl))
z.proposal <- mvnorm(S=diag(c(0.002,0.002)),n=nrow(twl)-1)


###################################################
### code chunk number 8: Movement.Rnw:157-159
###################################################
## Short test run
fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=300,thin=20)


###################################################
### code chunk number 9: Movement.Rnw:162-174
###################################################
## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x,discard=100),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z,discard=100),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)
## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)


###################################################
### code chunk number 10: LFtrace1 (eval = FALSE)
###################################################
## opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
## matplot(scale(t(fit$x[150,,]),scale=F),type="l",lty=1,col=c(2,4),
##         xlab="",ylab=expression(x[150]))
## matplot(scale(t(fit$z[150,,]),scale=F),type="l",lty=1,col=c(2,4),
##         xlab="",ylab=expression(z[150]))
## par(opar)


###################################################
### code chunk number 11: Movement.Rnw:188-189
###################################################
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
matplot(scale(t(fit$x[150,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(x[150]))
matplot(scale(t(fit$z[150,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(z[150]))
par(opar)


###################################################
### code chunk number 12: Movement.Rnw:197-203
###################################################
## Draw final sample
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=2000,thin=50)


###################################################
### code chunk number 13: LFci (eval = FALSE)
###################################################
## ## Plot time series of lat and lon
## opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
## s <- location.summary(fit$x)
## matplot(cbind(twl$Lon,s[,c("Lon.mean","Lon.2.5%","Lon.97.5%")]),type="l",lty=1,
##         col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lon")
## matplot(cbind(twl$Lat,s[,c("Lat.mean","Lat.2.5%","Lat.97.5%")]),type="l",lty=1,
##         col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lat")
## par(opar)


###################################################
### code chunk number 14: Movement.Rnw:227-228
###################################################
## Plot time series of lat and lon
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
s <- location.summary(fit$x)
matplot(cbind(twl$Lon,s[,c("Lon.mean","Lon.2.5%","Lon.97.5%")]),type="l",lty=1,
        col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lon")
matplot(cbind(twl$Lat,s[,c("Lat.mean","Lat.2.5%","Lat.97.5%")]),type="l",lty=1,
        col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lat")
par(opar)


###################################################
### code chunk number 15: LFpath (eval = FALSE)
###################################################
## ## Plot sequence of mean x's
## plot(location.mean(fit$x),pch=16,cex=0.4,col="grey80",xlab="Lon",ylab="Lat")
## points(t(fit$x[30,,]),pch=16,cex=0.2,col="dodgerblue1")
## points(t(fit$x[80,,]),pch=16,cex=0.2,col="dodgerblue1")
## points(t(fit$x[200,,]),pch=16,cex=0.2,col="dodgerblue1")
## points(t(fit$z[150,,]),pch=16,cex=0.2,col="firebrick1")
## points(t(fit$z[100,,]),pch=16,cex=0.2,col="firebrick1")
## points(t(fit$z[220,,]),pch=16,cex=0.2,col="firebrick1")
## lines(twl$Lon,twl$Lat,col="grey50")


###################################################
### code chunk number 16: Movement.Rnw:251-252
###################################################
## Plot sequence of mean x's
plot(location.mean(fit$x),pch=16,cex=0.4,col="grey80",xlab="Lon",ylab="Lat")
points(t(fit$x[30,,]),pch=16,cex=0.2,col="dodgerblue1")
points(t(fit$x[80,,]),pch=16,cex=0.2,col="dodgerblue1")
points(t(fit$x[200,,]),pch=16,cex=0.2,col="dodgerblue1")
points(t(fit$z[150,,]),pch=16,cex=0.2,col="firebrick1")
points(t(fit$z[100,,]),pch=16,cex=0.2,col="firebrick1")
points(t(fit$z[220,,]),pch=16,cex=0.2,col="firebrick1")
lines(twl$Lon,twl$Lat,col="grey50")


###################################################
### code chunk number 17: Movement.Rnw:282-283
###################################################
ps <- cbind(ave(ps[,1],rep(1:50,each=3)),ave(ps[,2],rep(1:50,each=3)))


###################################################
### code chunk number 18: Movement.Rnw:286-294
###################################################
library(SGAT)
## Resample to two minute intervals and compute zenith angles
tms <- as.POSIXct("2013-04-01","GMT")+(24*60*60)*days
tms.out <- seq(min(tms),max(tms),by=2*60)
d <- zenith.simulate(tms,ps[,1],ps[,2],tms.out)
## Compute times and locations at which twilight is observed.
twl <- twilight.simulate(d)
twl <- twilight.perturb(twl,rlnorm(nrow(twl),0.3,0.6))


###################################################
### code chunk number 19: Movement.Rnw:301-303
###################################################
## Initial x locations
x0 <- threshold.path(twl$Twilight,twl$Rise)$x


###################################################
### code chunk number 20: Movement.Rnw:306-317
###################################################
model <- threshold.model(twl$Twilight,twl$Rise,
                         twilight.model="LogNormal",
                         alpha=c(0.3,0.6),beta=c(8,3.5),
                         x0=NULL,z0=NULL)
r <- model$residuals(x0)
x0 <- cbind(x0[,1]-ifelse(twl$Rise,1,-1)*pmin(0,r-0.1)/1440*360,x0[,2])
z0 <- trackMidpts(x0)
model <- threshold.model(twl$Twilight,twl$Rise,
                         twilight.model="LogNormal",
                         alpha=c(0.3,0.6),beta=c(8,3.5),
                         x0=x0,z0=z0)


###################################################
### code chunk number 21: Movement.Rnw:321-324
###################################################
## Define initial proposals
x.proposal <- mvnorm(S=diag(c(0.002,0.002)),n=nrow(twl))
z.proposal <- mvnorm(S=diag(c(0.002,0.002)),n=nrow(twl)-1)


###################################################
### code chunk number 22: Movement.Rnw:327-329
###################################################
## Short test run
fit <- estelle.metropolis(model,x.proposal,z.proposal,iters=300,thin=20)


###################################################
### code chunk number 23: Movement.Rnw:332-344
###################################################
## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x,discard=100),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z,discard=100),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)
## Tune proposals based on previous run
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=300,thin=20)


###################################################
### code chunk number 24: HFtrace1 (eval = FALSE)
###################################################
## opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
## matplot(scale(t(fit$x[150,,]),scale=F),type="l",lty=1,col=c(2,4),
##         xlab="",ylab=expression(x[150]))
## matplot(scale(t(fit$z[150,,]),scale=F),type="l",lty=1,col=c(2,4),
##         xlab="",ylab=expression(z[150]))
## par(opar)


###################################################
### code chunk number 25: Movement.Rnw:358-359
###################################################
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
matplot(scale(t(fit$x[150,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(x[150]))
matplot(scale(t(fit$z[150,,]),scale=F),type="l",lty=1,col=c(2,4),
        xlab="",ylab=expression(z[150]))
par(opar)


###################################################
### code chunk number 26: Movement.Rnw:367-373
###################################################
## Draw final sample
x.proposal <- mvnorm(chain.cov(fit$x),s=0.3)
z.proposal <- mvnorm(chain.cov(fit$z),s=0.3)
fit <- estelle.metropolis(model,x.proposal,z.proposal,
                          x0=chain.last(fit$x),z0=chain.last(fit$z),
                          iters=2000,thin=50)


###################################################
### code chunk number 27: HFci (eval = FALSE)
###################################################
## ## Plot time series of lat and lon
## opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
## s <- location.summary(fit$x)
## matplot(cbind(twl$Lon,s[,c("Lon.mean","Lon.2.5%","Lon.97.5%")]),type="l",lty=1,
##         col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lon")
## matplot(cbind(twl$Lat,s[,c("Lat.mean","Lat.2.5%","Lat.97.5%")]),type="l",lty=1,
##         col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lat")
## par(opar)


###################################################
### code chunk number 28: Movement.Rnw:397-398
###################################################
## Plot time series of lat and lon
opar <- par(mfrow=c(2,1),mar=c(3,5,2,1)+0.1)
s <- location.summary(fit$x)
matplot(cbind(twl$Lon,s[,c("Lon.mean","Lon.2.5%","Lon.97.5%")]),type="l",lty=1,
        col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lon")
matplot(cbind(twl$Lat,s[,c("Lat.mean","Lat.2.5%","Lat.97.5%")]),type="l",lty=1,
        col=c("grey70","firebrick1","dodgerblue1","dodgerblue1"),ylab="Lat")
par(opar)


###################################################
### code chunk number 29: HFpath (eval = FALSE)
###################################################
## ## Plot sequence of mean x's
## plot(location.mean(fit$x),pch=16,cex=0.4,col="grey80",xlab="Lon",ylab="Lat")
## points(t(fit$x[30,,]),pch=16,cex=0.2,col="dodgerblue1")
## points(t(fit$x[80,,]),pch=16,cex=0.2,col="dodgerblue1")
## points(t(fit$x[200,,]),pch=16,cex=0.2,col="dodgerblue1")
## points(t(fit$z[150,,]),pch=16,cex=0.2,col="firebrick1")
## points(t(fit$z[100,,]),pch=16,cex=0.2,col="firebrick1")
## points(t(fit$z[220,,]),pch=16,cex=0.2,col="firebrick1")
## lines(twl$Lon,twl$Lat,col="grey50")


###################################################
### code chunk number 30: Movement.Rnw:421-422
###################################################
## Plot sequence of mean x's
plot(location.mean(fit$x),pch=16,cex=0.4,col="grey80",xlab="Lon",ylab="Lat")
points(t(fit$x[30,,]),pch=16,cex=0.2,col="dodgerblue1")
points(t(fit$x[80,,]),pch=16,cex=0.2,col="dodgerblue1")
points(t(fit$x[200,,]),pch=16,cex=0.2,col="dodgerblue1")
points(t(fit$z[150,,]),pch=16,cex=0.2,col="firebrick1")
points(t(fit$z[100,,]),pch=16,cex=0.2,col="firebrick1")
points(t(fit$z[220,,]),pch=16,cex=0.2,col="firebrick1")
lines(twl$Lon,twl$Lat,col="grey50")


