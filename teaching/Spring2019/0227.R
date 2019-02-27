library(astsa)

###################
# Feb 27, 2019
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
# djia = read.zoo('DJI.csv')
djia = DJI
djiar = diff(log(djia$DJI.Close))[-1]         # approximate returns
plot(djiar, main="DJIA Returns", type="n")  
lines(djiar) 
dev.off()

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
# overlapping view
par(mfrow=c(1,1))
tsplot(w, col=2, ylim=c(-3,3), main="moving average")
lines(v)

# now try this (not in text):  
# dev.new()  # open a new graphic device
# ts.plot(w, v, lty=2:1, col=1:2, lwd=1:2)
dev.off()


####################
# Example 1.10 Autoregression
w = rnorm(550,0,1)  # 50 extra to avoid startup problems
x = filter(w, filter=c(1,-.9), method="recursive")[-(1:50)]
tsplot(x, main="autoregression")
lines(w[-(1:50)], col=2)
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
