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

####################
# Example 3.16 ACF and PACF of AR(p)
ar2.acf = ARMAacf(ar=c(1.5,-.75), ma=0, 24)[-1]
ar2.pacf = ARMAacf(ar=c(1.5,-.75), ma=0, 24, pacf=TRUE)
par(mfrow=c(1,2))
plot(ar2.acf, type="h", xlab="lag")
abline(h=0)
plot(ar2.pacf, type="h", xlab="lag")
abline(h=0)
dev.off()

####################
# Example 3.18 Preliminary analysis of Rec
acf2(rec, 48)     # will produce values and a graphic 
(regr = ar.ols(rec, order=2, demean=F, intercept=TRUE))  # regression
regr$asy.se.coef  # standard errors                             
dev.off()



