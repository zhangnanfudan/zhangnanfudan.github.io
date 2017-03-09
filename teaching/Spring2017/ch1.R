library('astsa')
source('grid.r')  # I changed the defaults so that
#  grid() gives grid(lty=1, col = gray(.9))

###################
# Mar 2, 2017
###################
################
# Johnson & Johnson Quarterly Earnings
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(jj, ylab="Quarterly Earnings per Share",   type='n')
grid(lty=1)
lines(jj, type="o")

dev.off()


###########
# Global Warming
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(globtemp, ylab="Global Temperature Deviations",    type='n')   
grid(lty=1); lines(globtemp, type='o')
dev.off()

###########
# New York Stock Exchange
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(nyse,   type='n', ylab='NYSE Returns')
grid(lty=1);   lines(nyse); 
dev.off()

###########
# Speech
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(speech,   type='n', ylab='speech')
grid(lty=1);   lines(speech); 
dev.off()


#########
# SOI and Fish
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
plot(soi, ylab="", xlab="", main="Southern Oscillation Index",   type='n')
grid(lty=1); lines(soi)
plot(rec, ylab="", main="Recruitment",  type='n')
grid(lty=1); lines(rec)
dev.off()


#########
# fMRI
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
ts.plot(fmri1[,2:5], ylab="BOLD", xlab="", main="Cortex", type='n')
grid(lty=1); par(new=TRUE)
ts.plot(fmri1[,2:5], col=1:4, ylab="BOLD", xlab="", main="Cortex")
#
ts.plot(fmri1[,6:9], ylab="BOLD", xlab="", main="Thalamus & Cerebellum", type='n')
grid(lty=1); par(new=TRUE)
ts.plot(fmri1[,6:9], col=1:4, ylab="BOLD", xlab="", main="Thalamus & Cerebellum")
mtext("Time (1 pt = 2 sec)", side=1, line=1.5)
dev.off()


#####################
# Earthquakes and Explosions
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
plot(EQ5, main="Earthquake", xlab="", type='n')
grid(lty=1); lines(EQ5)
plot(EXP6, main="Explosion", xlab="", type='n')  
grid(lty=1); lines(EXP6)
mtext("Time", side=1, line=1.5)
dev.off()


####################
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
set.seed(1)
w = rnorm(500,0,1)                    # 500 N(0,1) variates
v = filter(w, sides=2, filter=rep(1/3,3))  # moving average
plot.ts(w, main="white noise", type='n')
grid(lty=1, col=gray(.9)); lines(w)
plot.ts(v, ylim=c(-3,3), main="moving average", type='n')
grid(lty=1, col=gray(.9)); lines(v)
dev.off()


####################
# Autoregression
par(mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0))
w = rnorm(550,0,1)  # 50 extra to avoid startup problems
x = filter(w, filter=c(1,-.9), method="recursive")[-(1:50)]
plot.ts(x, main="autoregression", type='n')
grid(lty=1, col=gray(.9)); lines(x)
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


####################
# Signal in Noise
par(mfrow = c(3,1), mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
cs = 2*cos(2*pi*1:500/50 + .6*pi);  w = rnorm(500,0,1)
par(mfrow=c(3,1), mar=c(3,2,2,1), cex.main=1.05)
plot.ts(cs, ylab='',xlab='', main=expression(2*cos(2*pi*t/50+.6*pi)), type='n', cex.main=1.5)
 grid(lty=1, col=gray(.9)); lines(cs) 
plot.ts(cs+w, ylab='',xlab='',main=expression(2*cos(2*pi*t/50+.6*pi) + N(0,1)), type='n', cex.main=1.5)
 grid(lty=1, col=gray(.9)); lines(cs+w) 
plot.ts(cs+5*w, ylab='', main=expression(2*cos(2*pi*t/50+.6*pi) + N(0,5^2)), type='n', cex.main=1.5)
 grid(lty=1, col=gray(.9)); lines(cs+5*w) 
dev.off()

###################
###################
###################
###################
# Mar 2, 2017
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



