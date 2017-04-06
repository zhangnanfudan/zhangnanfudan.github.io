library('astsa')
source('grid.r')

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