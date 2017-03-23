######################
# MA smoothing for mortality
par(mar=c(1,2.2,0,0)+.5, mgp=c(1.5,.6,0))
ma5 = filter(cmort, sides=2, rep(1,5)/5)
ma53 = filter(cmort, sides=2, rep(1,53)/53)
plot(cmort, type="n", ylab="mortality")
grid(lty=1); points(cmort)
lines(ma5, col=2, lwd=2); lines(ma53, col=4, lwd=2)
dev.off()


# Regression smoothing for mortality
par(mar=c(1,2.2,0,0)+.5, mgp=c(1.5,.6,0))
wk = time(cmort) - mean(time(cmort))
wk2 = wk^2; wk3 = wk^3
cs = cos(2*pi*wk); sn = sin(2*pi*wk)
reg1 = lm(cmort~wk + wk2 + wk3, na.action=NULL)
reg2 = lm(cmort~wk + wk2 + wk3 + cs + sn, na.action=NULL)
plot(cmort, type="n", ylab="mortality")
grid(lty=1); points(cmort)
lines(fitted(reg1), col=2, lwd=2); lines(fitted(reg2), col=4, lwd=2)
dev.off()

# Kernel smoothing
par(mar=c(1,2.2,0,0)+.5, mgp=c(1.5,.6,0))
plot(cmort, type="n", ylab="mortality")
grid(lty=1); points(cmort)
lines(ksmooth(time(cmort), cmort, "normal", bandwidth=5/52), 
      lwd=2, col=2)
lines(ksmooth(time(cmort), cmort, "normal", bandwidth=2),
      lwd=2, col=4)
dev.off()

par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,.5), mgp=c(1.6,.6,0))
plot(cmort, type="p", ylab="mortality", main="nearest neighbor")
lines(supsmu(time(cmort), cmort, span=.5), lwd=2, col=2)
lines(supsmu(time(cmort), cmort, span=.01),lwd=2, col=4)
plot(cmort, type="p", ylab="mortality", main="lowess")
lines(lowess(cmort, f=.02), lwd=2, col=4)
lines(lowess(cmort, f=2/3), lwd=2, col=2)
dev.off()

# smoothing splines
par(mar=c(2.5,2.5,.5,.5), mgp=c(1.6,.6,0))
plot(cmort, type="n", ylab="mortality")
grid(lty=1); points(cmort)
lines(smooth.spline(time(cmort), cmort), lwd=2, col=2)
lines(smooth.spline(time(cmort), cmort, spar=1), lwd=2, col=4)
dev.off()

#############################################
# Regress motality and temperature
par(mfrow=c(2,1),mar=c(2.5,2.5,2.5,.5), mgp=c(1.6,.6,0))
plot(tempr, cmort, main="lowess", type='n', xlab="Temperature", 
     ylab="Mortality")
grid(lty=1); points(tempr,cmort, col=gray(.5))
lines(lowess(tempr,cmort), col=4, lwd=2)

plot(tempr, cmort, main="smoothing splines", type='n', xlab="Temperature", 
     ylab="Mortality")
grid(lty=1); points(tempr,cmort, col=gray(.5))
lines(smooth.spline(tempr, cmort), col=4, lwd=2)
dev.off()




##############################################
##############################################
library('astsa')
source('grid.r')

##########
# AR(1) with phi=+-0.9
set.seed(101010)
par(mfrow = c(2,1), mar=c(1.5,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
x<-arima.sim(list(order=c(1,0,0), ar=.9), n=100)
plot(x, ylab="x", xlab="", main=(expression(AR(1)~~~phi==+.9)), type='n')
grid(lty=1)
lines(x)
x<-arima.sim(list(order=c(1,0,0), ar=-.9), n=100)
plot(x, ylab="x",  xlab="",  main=(expression(AR(1)~~~phi==-.9)), type='n')
grid(lty=1)
lines(x)
mtext('Time', side=1, line=1)
dev.off()


#####################
# MA(1) with theta=+-0.9
par(mfrow = c(2,1), mar=c(1.5,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
set.seed(101010)
plot(x<-arima.sim(list(order=c(0,0,1), ma=.9), n=100), ylab="x", xlab="", 
     main=(expression(MA(1)~~~theta==+.9)), type='n')
grid(lty=1)
lines(x)
plot(x<-arima.sim(list(order=c(0,0,1), ma=-.9), n=100), ylab="x", xlab='', 
     main=(expression(MA(1)~~~theta==-.9)), type='n')
grid(lty=1)
lines(x)
mtext('Time', side=1, line=1)
dev.off()