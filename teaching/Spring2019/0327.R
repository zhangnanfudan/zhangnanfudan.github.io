library(astsa)

####################
# Chapter 3
####################
# Example 3.2 AR(p)
par(mfrow=c(2,1))                         
# in the expressions below, ~ is a space and == is equal
tsplot(arima.sim(list(order=c(1,0,0), ar=.9), n=100), ylab="x", main=(expression(AR(1)~~~phi==+.9))) 
tsplot(arima.sim(list(order=c(1,0,0), ar=-.9), n=100), ylab="x", main=(expression(AR(1)~~~phi==-.9))) 
dev.off()

####################
# Example 3.5 MA(1)
par(mfrow=c(2,1))                                   
tsplot(arima.sim(list(order=c(0,0,1), ma=.9), n=100), ylab="x", main=(expression(MA(1)~~~theta==+.9)))    
tsplot(arima.sim(list(order=c(0,0,1), ma=-.9), n=100), ylab="x", main=(expression(MA(1)~~~theta==-.9)))    

dev.off()

####################
# Example 3.11 AR(2) with complex roots
par(mfrow=c(3,1))
ACF = ARMAacf(ar=c(1,-.25), ma=0, 50)
plot(ACF, type="h", xlab="lag", main="equal real roots")
abline(h=0)

ACF = ARMAacf(ar=c(.75,-.125), ma=0, 50)
plot(ACF, type="h", xlab="lag", main="distinct real roots")
abline(h=0)

ACF = ARMAacf(ar=c(1.5,-.75), ma=0, 50)
plot(ACF, type="h", xlab="lag", main="complex roots")
abline(h=0)

dev.off()
