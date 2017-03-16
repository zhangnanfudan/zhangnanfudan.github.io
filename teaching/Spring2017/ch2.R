library('astsa')
source('grid.r')


###############
# Global temperature with linear trend
par(mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
summary(fit <- lm(gtemp~time(gtemp))) # regress gtemp on time
plot(gtemp, type="n", ylab="Global Temperature Deviation")
grid(lty=1)
lines(gtemp, type="o")
abline(fit) # add regression line to the plot
dev.off()


################
# cardiovascular mortality, temperature, polution in LA, 1970-1979
par(mfrow = c(3,1), mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
plot(cmort, main="Cardiovascular Mortality", xlab="", ylab="",   type='n')
grid(lty=1); lines(cmort)
plot(tempr, main="Temperature", xlab="", ylab="",   type='n')
grid(lty=1); lines(tempr)
plot(part, main="Particulates", xlab="", ylab="",   type='n')
grid(lty=1); lines(part)
dev.off()


###########
# Scatterplot matrix
pairs(cbind(Mortality=cmort, Temperature=tempr, Particulates=part))
dev.off()


temp2 = temp^2
trend = time(cmort) # time
num=length(trend)
fit1 = lm(cmort~ trend, na.action=NULL)
fit2 = lm(cmort~ trend + temp, na.action=NULL)
fit3 = lm(cmort~ trend + temp + temp2, na.action=NULL)
fit4 = lm(cmort~ trend + temp + temp2 + part, na.action=NULL)

sum1=summary(fit1) # regression results
sum2=summary(fit2) # regression results
sum3=summary(fit3) # regression results
sum4=summary(fit4) # regression results

# R.squared
rsquared=c(lm1=sum1$r.squared, 
           lm2=sum2$r.squared, 
           lm3=sum3$r.squared, 
           lm4=sum4$r.squared)
rsquared

# AIC
aic=c(lm1=AIC(fit1), 
      lm2=AIC(fit2), 
      lm3=AIC(fit3), 
      lm4=AIC(fit4))/num- log(2*pi)
aic

# BIC
bic=c(lm1=AIC(fit1, k=log(num)), 
      lm2=AIC(fit2, k=log(num)), 
      lm3=AIC(fit3, k=log(num)), 
      lm4=AIC(fit4, k=log(num)))/num- log(2*pi)
bic

################
# Exploratory data analysis
##########
# Data examples
##########
# Johnson & Johnson Quarterly Earnings
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(jj, ylab="Quarterly Earnings per Share",   type='n')
grid(lty=1)
lines(jj, type="o")
dev.off()

###########
# Global temperature
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(globtemp, ylab="Global Temperature Deviations",    type='n')   
grid(lty=1); lines(globtemp, type='o')
dev.off()

#################
# RW with Drift
par(mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
set.seed(154)                # so you can reproduce the results
w = rnorm(200,0,1);  x = cumsum(w)   # two commands in one line
wd = w +.2;   xd = cumsum(wd)
plot.ts(xd, ylim=c(-5,55), main="random walk", ylab='',   type='n')
grid(lty=1); lines(xd)
lines(x, col=4); abline(h=0, col=4, lty=2)
abline(a=0, b=.2, lty=2)
dev.off()
#################

##########
# Detrend global temperature
fit = lm(gtemp~time(gtemp), na.action=NULL) # regress gtemp on time
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
plot(resid(fit), type="n", main="detrended")
grid(lty=1); lines(resid(fit), type='o')
plot(diff(gtemp), type="n", main="first difference")
grid(lty=1); lines(diff(gtemp), type='o')
dev.off()

mean(diff(gtemp))
sd(diff(gtemp))/sqrt(length(diff(gtemp)))
# 95% C.I. contains 0
c(mean(diff(gtemp))-1.96*sd(diff(gtemp))/sqrt(length(diff(gtemp))), 
  mean(diff(gtemp))+1.96*sd(diff(gtemp))/sqrt(length(diff(gtemp))))


# plot ACFs
par(mfrow = c(3,1), mar=c(2,2,.75,0)+.5, mgp=c(1.6,.6,0))
acf(gtemp, 48, xlab="", main="", panel.first=grid(lty=1))
mtext("gtemp", side=3, line=.1, cex=1, font=2)
acf(resid(fit), 48, xlab="", main="", panel.first=grid(lty=1))
mtext("detrended", side=3, line=.1, cex=1, font=2)
acf(diff(gtemp), 48, xlab="", main="first difference", 
    panel.first=grid(lty=1))
mtext("first difference", side=3, line=.1, cex=1, font=2)
mtext("LAG", side=1, line=1.5, cex=.76)
dev.off()

############
# Paleoclimatic glacial varves
par(mfrow = c(2,1), mar=c(2,1.5,.75,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
plot(varve, main="varve", ylab="", xlab="",  type='n')
grid(lty=1); lines(varve)
plot(log(varve), main="log(varve)", ylab="",   type='n')
grid(lty=1); lines(log(varve))
dev.off()


##############
# Scatterplot matrix of SOI
lag1.plot(soi, 12)

acf(soi, 48, xlab=' ', main="", panel.first=grid(lty=1))
mtext(side=3, 'Southern Oscillation Index', font=2)
dev.off()

# Scatterplot matrix of Recruitment v.s. SOI
lag2.plot(soi, rec, 8)
dev.off()


#####################
# Signal in noise: periodic behavior 
par(mfrow=c(2,1), mar=c(1,2.2,0,0)+.5, mgp=c(1.5,.6,0))
set.seed(1000)  # so you can reproduce these results
x = 2*cos(2*pi*1:500/50 + .6*pi) + rnorm(500,0,5)
z1 = cos(2*pi*1:500/50)
z2 = sin(2*pi*1:500/50)
fit <- lm(x~0+z1+z2) # zero to exclude the intercept, 
plot.ts(x, xlab='', type='n')
grid(); lines(x)
plot.ts(x, ylab=expression(hat(x)), type='n') 
lines(x, col=gray(.5))
grid();
lines(fitted(fit), col=2, lwd=2)
lines(2*cos(2*pi*1:500/50 + .6*pi), col=4, lwd=2)
dev.off()


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
 