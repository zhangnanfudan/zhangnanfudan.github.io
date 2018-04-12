library(astsa)

####################
# Example 2.11 Moving Average Smoother
wgts = c(.5, rep(1,11), .5)/12
soif = filter(soi, sides=2, filter=wgts)
tsplot(soi)
lines(soif, lwd=2, col=4)
par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert
nwgts = c(rep(0,20), wgts, rep(0,20))
plot(nwgts, type="l", ylim = c(-.02,.1), xaxt='n', yaxt='n', ann=FALSE)
dev.off()

####################
# Example 2.12 Kernel Smoothing
tsplot(soi)
lines(ksmooth(time(soi), soi, "normal", bandwidth=1), lwd=2, col=4)
par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert
gauss = function(x) { 1/sqrt(2*pi) * exp(-(x^2)/2) }
x = seq(from = -3, to = 3, by = 0.001)
plot(x, gauss(x), type ="l", ylim=c(-.02,.45), xaxt='n', yaxt='n', ann=FALSE)
dev.off()

####################
# Example 2.13
tsplot(soi)
lines(lowess(soi, f=.05), lwd=2, col=4) # El Nino cycle
lines(lowess(soi), lty=2, lwd=2, col=2) # trend (with default span)
dev.off()

####################
# Example 2.14
tsplot(soi)
lines(smooth.spline(time(soi), soi, spar=.5), lwd=2, col=4)
lines(smooth.spline(time(soi), soi, spar= 1), lty=2, lwd=2, col=2)
dev.off()

####################
# Example 2.15
plot(tempr, cmort, xlab="Temperature", ylab="Mortality")
lines(lowess(tempr, cmort))  
dev.off()



####################
# Chapter 3
####################
# Example 3.2
par(mfrow=c(2,1))                         
# in the expressions below, ~ is a space and == is equal
tsplot(arima.sim(list(order=c(1,0,0), ar=.9), n=100), ylab="x", main=(expression(AR(1)~~~phi==+.9))) 
tsplot(arima.sim(list(order=c(1,0,0), ar=-.9), n=100), ylab="x", main=(expression(AR(1)~~~phi==-.9))) 
dev.off()

####################
# Example 3.5
par(mfrow=c(2,1))                                   
tsplot(arima.sim(list(order=c(0,0,1), ma=.9), n=100), ylab="x", main=(expression(MA(1)~~~theta==+.9)))    
tsplot(arima.sim(list(order=c(0,0,1), ma=-.9), n=100), ylab="x", main=(expression(MA(1)~~~theta==-.9)))    
dev.off()