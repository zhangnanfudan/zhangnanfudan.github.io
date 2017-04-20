library('astsa')
source('grid.r')

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