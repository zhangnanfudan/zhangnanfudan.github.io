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
plot(x<-arima.sim(list(order=c(0,0,1), ma=.9), n=100), 
     ylab="x", xlab="", main=(expression(MA(1)~~~theta==+.9)), type='n')
grid(lty=1)
lines(x)
plot(x<-arima.sim(list(order=c(0,0,1), ma=-.9), n=100), 
     ylab="x", xlab='', main=(expression(MA(1)~~~theta==-.9)), type='n')
grid(lty=1)
lines(x)
mtext('Time', side=1, line=1)
dev.off()

#####################
#####################
#####################


################
# ACF and PACF of AR(2)
ACF = ARMAacf(ar=c(1.5,-.75), ma=0, 24)[-1]
PACF = ARMAacf(ar=c(1.5,-.75), ma=0, 24, pacf=TRUE) 
par(mfrow=c(1,2), mar=c(2.5,2.5,.5,0)+.5, mgp=c(1.6,.6,0))
plot(ACF, type="h", xlab="lag", ylim=c(-.8,1), panel.first=grid(lty=1)); abline(h=0)
plot(PACF, type="h", xlab="lag", ylim=c(-.8,1), panel.first=grid(lty=1)); abline(h=0)
dev.off()


#######################
# ACF and PACF of Recruitment data
u = acf2(rec, 48); dev.off(); ACF=u[,1]; PACF=u[,2]
LAG = 1:48/frequency(rec)
num = length(rec)
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
minu = min(minA, minP, L) - 0.01
maxu = min(max(maxA + 0.1, maxP + 0.1), 1)
par(mfrow=c(2,1), mar=c(2,2,0,0)+.5, mgp=c(1.5,.6,0))
plot(LAG, ACF, type="h", xlab="LAG", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(LAG, PACF, type="h", xlab="LAG",  ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
dev.off()

#####################
#####################
#####################

############
# Recruitment series AR(2) forecasting
par(mar=c(2.5,2.5,0,0)+.5, mgp=c(1.6,.6,0))
regr = ar.ols(rec, order=2, demean=FALSE, intercept=TRUE)
fore = predict(regr, n.ahead=24)
ts.plot(rec, fore$pred, col=1:2, xlim=c(1980,1990), ylab="Recruitment", type='n')
grid(lty=1); par(new=TRUE)
ts.plot(rec, fore$pred, col=1:2, xlim=c(1980,1990), ylab="Recruitment")
 U = fore$pred+fore$se
 L = fore$pred-fore$se	
 xx = c(time(U), rev(time(U)))
 yy = c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
lines(fore$pred, type="p", col=2)
dev.off()


#####################
#####################
#####################
#####################
# Example 3.38: Analysis of GNP Data
layout(matrix(c(1, 1, 2), ncol = 1))
par(mar=c(2.75,2.5,.5,.5), mgp=c(1.6,.6,0), cex.lab=1.1) 
plot(gnp, ylab="Billions of Dollars",  type='n')
grid(lty=1, col=gray(.9)); lines(gnp)
# acf  
acf(gnp, 48, panel.first=grid(lty=1))
dev.off()


### Plots of difference of GNP itself and of logGNP (growth rate)
par(mfrow=c(2,1), mar=c(2.75,2.5,.5,.5), mgp=c(1.6,.6,0)) 
plot(diff(gnp), ylab="diff(GNP)", type='n')
grid(lty=1, col=gray(.9)); lines(diff(gnp))
abline(h=mean(diff(gnp)), col=4)

plot(diff(log(gnp)), ylab="GNP Growth Rate", type='n')
grid(lty=1, col=gray(.9)); lines(diff(log(gnp)))
abline(h=mean(diff(log(gnp))), col=4)
dev.off()


### ACF and PACF
ACF = acf(diff(log(gnp)), 24, plot=FALSE)$acf[-1]
PACF = pacf(diff(log(gnp)), 24, plot=FALSE)$acf
num = length(gnp)-1
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
LAG = 1:24/4
minu = min(minA, minP, L) - 0.01
maxu = min(maxA + 0.2, maxP + 0.2, 1)
par(mfrow=c(2,1), mar=c(2,2.5,0,0)+.5, mgp=c(1.4,.6,0))
plot(LAG, ACF, type="h", xlab="LAG", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(LAG, PACF, type="h", xlab="LAG",  ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
dev.off()

#####################
# Example 3.39: Diagnostics for GNP Growth Rate
sarima(diff(log(gnp)), 0, 0, 2) # MA(2)
dev.off()



#####################
# Example 3.39: Diagnostics for the Glacial Varve Series
par(mfrow=c(2,1), mar=c(2,2.5,.5,0)+.5, mgp=c(1.4,.6,0))
rs = resid(arima(log(varve), order=c(0,1,1)))
  pval=c()
  nlag=20
  ppq=1
  for (i in (ppq+1):nlag) {
        u <- stats::Box.test(rs, i, type = "Ljung-Box")$statistic
        pval[i] <- stats::pchisq(u, i - ppq, lower.tail = FALSE)
    }
    plot((ppq + 1):nlag, pval[(ppq + 1):nlag], xlab = "lag", cex.main=1, font.main=1,
        ylab = "p value", ylim = c(0, 1), main = "p values for Ljung-Box statistic")
    abline(h = 0.05, lty = 2, col = "blue")
#
rs = resid(arima(log(varve), order=c(1,1,1)))
  pval=c()
  nlag=20
  ppq=2
  for (i in (ppq+1):nlag) {
        u <- stats::Box.test(rs, i, type = "Ljung-Box")$statistic
        pval[i] <- stats::pchisq(u, i - ppq, lower.tail = FALSE)
    }
    plot((ppq + 1):nlag, pval[(ppq + 1):nlag], xlab = "lag", 
        ylab = "p value", ylim = c(0, 1), main = "")
    abline(h = 0.05, lty = 2, col = "blue")
dev.off()
	

#####################
# Example 3.41: A Problem with Overfitting in US Population Data
par(mar=c(2,2.5,.5,0)+.5, mgp=c(1.6,.6,0))
dat = read.table("uspop.dat")
y = dat[,2]
x = dat[,1]
b = dat[,3]
g = function(x) (
  b[1]+b[2]*(x-1955)+b[3]*(x-1955)^2+b[4]*(x-1955)^3+b[5]*(x-1955)^4+b[6]*(x-1955)^5+b[7]*(x-1955)^6+b[8]*(x-1955)^7+b[9]*(x-1955)^8)/10^8
curve(g, 1910,2002, ylab="Population", xlab="Year", main="U.S. Population by Official Census", panel.first=grid(ny=NULL,lty=1), cex.main=1, font.main=1)
abline(v=c(1910,1930,1950,1970,1990), lty=1, col='lightgray')
points(x, y/10^8, pch=16)
mtext(expression(""%*% 10^8), side=2, line=1.5, adj=.95)
axis(1, seq(1910,1990,by=10), labels=FALSE)
dev.off()

#####################
# Example 3.42: Model Choice for the US GNP Series
gnpgr=diff(log(gnp))
sarima(gnpgr, 1, 0, 0) # AR(1)
sarima(gnpgr, 0, 0, 2) # MA(2)


####################
####################
####################
####################
pdf(file="corerr1.pdf",width=7, height=4)
par(mfrow=c(2,1), mar=c(2.5,2.5,0,0)+.5, mgp=c(1.5,.6,0))
trend = time(cmort); temp = tempr - mean(tempr); temp2 = temp^2
fit <- lm(cmort~trend + temp + temp2 + part, na.action=NULL)
ACF = acf(resid(fit), 52, plot=FALSE)$acf[-1]
PACF = pacf(resid(fit), 52, plot=FALSE)$acf
num = length(cmort)
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
LAG = 1:52/52
minu = min(minA, minP, L) - 0.01
maxu = min(maxA + 0.2, maxP + 0.2, 1)
plot(LAG, ACF, type="h", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(LAG, PACF, type="h", ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
dev.off()



#######################
pdf(file="sAR1.pdf",width=7.5, height=4)
set.seed(666)
phi = c(rep(0,11),.9)
sAR = arima.sim(list(order=c(12,0,0), ar=phi), n=37)
sAR = ts(sAR, freq=12)
layout(matrix(c(1,1,2, 1,1,3), nc=2))
par(mar=c(2.5,2.5,2,1), mgp=c(1.6,.6,0))
plot(sAR, axes=FALSE, col='#808080', main='seasonal AR(1)', xlab="year", type='c')
abline(v=1:4, lty=2, col=gray(.6))
abline(h=seq(-4,2,2), col=gray(.9), lty=1)
Months = c("J","F","M","A","M","J","J","A","S","O","N","D")
points(sAR, pch=Months, cex=1.35, font=4, col=1:4) 
axis(1,1:4) 
axis(2)
box()
#
ACF = ARMAacf(ar=phi, ma=0, 100)[-1]  # [-1] removes 0 lag
PACF = ARMAacf(ar=phi, ma=0, 100, pacf=TRUE)
plot(ACF, type="h", xlab="LAG", ylim=c(-.1,1), axes=FALSE);
segments(0,0,0,1)
axis(1, seq(0,100,by=12))
axis(2)
box()
abline(h=0)
plot(PACF, type="h", xlab="LAG", ylim=c(-.1,1), axes=FALSE);
axis(1, seq(0,100,by=12))
axis(2)
box()
abline(h=0)
dev.off()


###################
pdf(file="sarmaacf.pdf",width=7.25,height=3.25)
phi = c(rep(0,11),.8)
ACF = ARMAacf(ar=phi, ma=-.5, 50)[-1]     # [-1] removes 0 lag
PACF = ARMAacf(ar=phi, ma=-.5, 50, pacf=TRUE)
par(mfrow=c(1,2), mar=c(2.5,2.5,2,1), mgp=c(1.6,.6,0))
plot(ACF,  type="h", xlab="LAG", ylim=c(-.4,.8), axes=FALSE)  
abline(h=0)
axis(1, seq(0,50,by=12))
axis(2)
box()
plot(PACF, type="h", xlab="LAG", ylim=c(-.4,.8), axes=FALSE)  
abline(h=0)
axis(1, seq(0,50,by=12))
axis(2)
box()
dev.off()


#####################
pdf(file="AirPdata.pdf", width=7.5, height=6)
par(mfrow=c(4,1), mar = c(0, 3, 0, 3), oma=c(3,0,2,0), mgp=c(1.6,.6,0), cex.lab=1.5)
x = AirPassengers
lx = log(x); dlx = diff(lx); ddlx = diff(dlx, 12)
u = ts.union(x,lx,dlx,ddlx)
plot.ts(u[,1], ylab='x', xaxt="no", type='n')
grid(lty=1, col=gray(.9)); lines(u[,1])
plot.ts(u[,2], ylab='lx', xaxt="no", type='n', yaxt='no', ylim=c(4.5,6.5))
grid(lty=1, col=gray(.9)); axis(4); lines(u[,2]) 
plot.ts(u[,3], ylab='dlx', xaxt="no", type='n')
grid(lty=1, col=gray(.9)); lines(u[,3])
plot.ts(u[,4], ylab='ddlx', yaxt='no', type='n')
grid(lty=1, col=gray(.9)); axis(4); lines(u[,4])
title(xlab="Time", outer=TRUE)
dev.off()


#################
pdf(file="AirPacf.pdf",width=7, height=3.5)
x = AirPassengers
lx = log(x); dlx = diff(lx); ddlx = diff(dlx, 12)
ACF = acf(ddlx, 50, plot=FALSE)$acf[-1]
PACF = pacf(ddlx, 50, plot=FALSE)$acf
num = length(x)-13
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
LAG = 1:50/12
minu = min(minA, minP, L) - 0.01
maxu = min(maxA + 0.2, maxP + 0.2, 1)
par(mfrow=c(2,1), mar=c(2,2.5,0,0)+.5, mgp=c(1.4,.6,0))
plot(LAG, ACF, type="h", xlab="LAG", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(LAG, PACF, type="h", xlab="LAG",  ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
dev.off()


##############
pdf(file="AirPresid.pdf",width=7.5)
sarima(log(AirPassengers), 0, 1, 1, 0, 1, 1, 12) 
dev.off()


####################
x = AirPassengers
xdata = log(x)
fore = sarima.for(xdata, 12, 0,1,1, 0,1,1,12)
 dev.off()
# -- for publication
pdf(file="AirPfore.pdf",width=7.5, height=3.5)
par(mar=c(2,2,0,0)+.5, mgp=c(1.4,.6,0))
    n = length(xdata)
    U = fore$pred + 2 * fore$se
    L = fore$pred - 2 * fore$se
	U1 = fore$pred + fore$se
    L1 = fore$pred - fore$se
    a = max(1, n - 100)
    minx = min(xdata[a:n], L)
    maxx = max(xdata[a:n], U)
	xnew = window(xdata, start=1953)
    ts.plot(xnew, fore$pred, col = 1:2, ylim = c(minx, maxx), type='n')
    grid(lty=1); par(new=TRUE)
	ts.plot(xnew, fore$pred, col = 1:2, type = "o", ylim = c(minx, maxx), ylab='log(AirPassengers)')
	xx = c(time(U), rev(time(U)))
    yy = c(L, rev(U))
    polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
    yy1 = c(L1, rev(U1))
    polygon(xx, yy1, border = 8, col = gray(0.6, alpha = 0.2))
    lines(fore$pred, col = "red", type = "o")	
dev.off()


###########################



