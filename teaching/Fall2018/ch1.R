library(astsa)

###################
# Mar 8, 2018
###################
################
# Johnson & Johnson Quarterly Earnings
tsplot(jj, type="o", ylab="Quarterly Earnings per Share")
dev.off()

###########
# Global Warming
tsplot(globtemp, type="o", ylab="Global Temperature Deviations")
dev.off()

###########
# Speech
tsplot(speech)  
dev.off()

########### 
# Dow Jones Industrial Average
library(quantmod)
getSymbols("^DJI", from="2006-04-20", to="2016-04-20", periodicity="daily")
djia = DJI
library(xts)                              # if don't use TTR
djiar = diff(log(djia$DJI.Close))[-1]         # approximate returns
plot(djiar, main="DJIA Returns", type="n")  
lines(djiar) 

###########
# New York Stock Exchange
# tsplot(nyse, ylab='NYSE Returns')
# dev.off()

#########
# SOI and Fish
par(mfrow = c(2,1))  # set up the graphics
tsplot(soi, ylab="", main="Southern Oscillation Index")
tsplot(rec, ylab="", main="Recruitment") 
dev.off()


#########
# fMRI
par(mfrow=c(2,1), mar=c(3,2,1,0)+.5, mgp=c(1.6,.6,0))  
ts.plot(fmri1[,2:5], col=1:4, ylab="BOLD", xlab="", main="Cortex")
ts.plot(fmri1[,6:9], col=1:4, ylab="BOLD", xlab="", main="Thalamus & Cerebellum")
mtext("Time (1 pt = 2 sec)", side=1, line=2) 
dev.off()


#####################
# Earthquakes and Explosions
par(mfrow=c(2,1))
tsplot(EQ5, main="Earthquake")
tsplot(EXP6, main="Explosion")
dev.off()


####################
# Example 1.9 white noise and moving average
w = rnorm(500,0,1)  # 500 N(0,1) variates
v = filter(w, sides=2, rep(1/3,3))  # moving average
par(mfrow=c(2,1))
tsplot(w, main="white noise")
tsplot(v, ylim=c(-3,3), main="moving average")

# now try this (not in text):  
# dev.new()  # open a new graphic device
# ts.plot(w, v, lty=2:1, col=1:2, lwd=1:2)
dev.off()


####################
# Example 1.10 Autoregression
w = rnorm(550,0,1)  # 50 extra to avoid startup problems
x = filter(w, filter=c(1,-.9), method="recursive")[-(1:50)]
tsplot(x, main="autoregression")
dev.off()


#################
# Example 1.11 RW with Drift
set.seed(154) # so you can reproduce the results
w = rnorm(200); x = cumsum(w) # two commands in one line
wd = w +.2;    xd = cumsum(wd)
tsplot(xd, ylim=c(-5,55), main="random walk", ylab='')
lines(x, col=4) 
abline(h=0, col=4, lty=2)
abline(a=0, b=.2, lty=2)
dev.off()


####################
# Example 1.12 Signal in Noise
cs = 2*cos(2*pi*(1:500)/50 + .6*pi)
w = rnorm(500,0,1)
par(mfrow=c(3,1), mar=c(3,2,2,1), cex.main=1.5)
tsplot(cs, ylab="", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)))
tsplot(cs + w, ylab="", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)+N(0,1)))
tsplot(cs + 5*w, ylab="", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)+N(0,25)))
dev.off()

###################
###################
###################
###################
# Mar 15, 2018
###################
# Autocovariance of three-point MA
ACF = c(0,0,0,1/3,2/3,1,2/3,1/3,0,0,0)/3
LAG = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
par(mar=c(2,2,0,0)+.5, mgp=c(1.6,.6,0))
plot(LAG, ACF, type='h', lwd=3,  panel.first=grid(lty=1))
abline(h=0)
points(LAG[-(4:8)],ACF[-(4:8)], pch=20)
dev.off()



###################
# CCF of y and x (leads)
par(mar=c(2.2,2,.7,0)+.5, mgp=c(1.6,.6,0))
set.seed(2)
x = rnorm(100)
y = lag(x,-5) + rnorm(100)
ccf(y,x, ylab='CCovF', xlab="LAG", type='covariance',panel.first=grid(lty=1, col=gray(.9)), main="", lwd=2)
abline(v=0, lty=2)
text(11, .9, 'x leads')
text(-9, .9, 'y leads')
title(main="y & x", cex.main=1)
dev.off()

########
# Extra example on ACF and CCF
set.seed(1492)
num=120; t=1:num
X = ts(2*cos(2*pi*t/12) + rnorm(num), freq=12)
Y = ts(2*cos(2*pi*(t+5)/12) + rnorm(num), freq=12)
Yw = resid( lm(Y~ cos(2*pi*t/12) + sin(2*pi*t/12), na.action=NULL) )
par(mfrow=c(3,2), mgp=c(1.6,.6,0), mar=c(3,3,1,1) )
plot(X, type='n'); grid(lty=1, col=gray(.9)); lines(X)
plot(Y, type='n'); grid(lty=1, col=gray(.9)); lines(Y)
acf(X,48, ylab='ACF(X)', panel.first=grid(lty=1, col=gray(.9)))
acf(Y,48, ylab='ACF(Y)', panel.first=grid(lty=1, col=gray(.9)))
ccf(X,Y,24, ylab='CCF(X,Y)', panel.first=grid(lty=1, col=gray(.9)))
ccf(X,Yw,24, ylab='CCF(X,Yw)', ylim=c(-.6,.6), panel.first=grid(lty=1, col=gray(.9)))
dev.off()


######################
# Speech data ACF
ACF = acf(speech, 250, plot = FALSE)$acf[-1]
LAG = 1:250
minA = min(ACF)
maxA = max(ACF)
num=length(speech)
U = 2/sqrt(num)
L = -U
minu = min(minA, L) - 0.01
maxu = min(maxA + 0.1, 1)
plot(LAG, ACF, type = "n", ylim = c(minu, maxu))
grid(lty = 1, col = gray(0.9))
abline(h = c(0, L, U), lty = c(1, 2, 2), col = c(1, 4, 4))
lines(LAG, ACF, type = "h")
#
acf(speech, 250, panel.first=grid(lty=1))
dev.off()

###################
# SOI Recruitment ACF and CCF
par(mfrow=c(3,1), mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
acf(soi, 48, xlab=' ', main="", panel.first=grid(lty=1))
mtext(side=3, 'Southern Oscillation Index', font=2)
acf(rec, 48, xlab='',main="", panel.first=grid(lty=1))
mtext(side=3, 'Recruitment', font=2)
ccf(soi, rec, 48, xlab='LAG', main="", ylab="CCF",panel.first=grid(lty=1))
mtext(side=3, 'SOI vs Recruitment', font=2)
dev.off()


###################
# Soil Temperature 
par(mar=c(0,1,0,0)+.5, mgp=c(1.6,.6,0))
persp(1:64, 1:36, soiltemp,  phi=25, theta=25, scale=FALSE, expand=4, ticktype="detailed", xlab="rows", ylab="cols", zlab="temperature")
dev.off()

##
# 64 row means
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot.ts(rowMeans(soiltemp), xlab="row", ylab="Average Temperature" , type='n')
grid(lty=1); lines(rowMeans(soiltemp))
dev.off()


#######################
# Spatial ACF
par(mar=c(1,1,0,0)+.5)
fs = Mod(fft(soiltemp-mean(soiltemp)))^2/(64*36)
cs = Re(fft(fs, inverse=TRUE)/sqrt(64*36)) # ACovF
rs = cs/cs[1,1] # ACF
rs2 = cbind(rs[1:41,21:2], rs[1:41,1:21])
rs3 = rbind(rs2[41:2,], rs2)
persp(-40:40, -20:20, rs3, phi=30, theta=30, expand=30, scale="FALSE",
ticktype="detailed", xlab="row lags", ylab="column lags", zlab="ACF")
dev.off()



