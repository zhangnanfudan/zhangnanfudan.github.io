library('astsa')
source('grid.r')

#######################
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