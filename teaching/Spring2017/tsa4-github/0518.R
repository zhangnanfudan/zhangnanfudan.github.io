library('astsa')
source('grid.r')

##############
# A Periodic Series
x1 = 2*cos(2*pi*1:100*6/100) + 3*sin(2*pi*1:100*6/100)
x2 = 4*cos(2*pi*1:100*10/100) + 5*sin(2*pi*1:100*10/100)
x3 = 6*cos(2*pi*1:100*40/100) + 7*sin(2*pi*1:100*40/100)
x = x1 + x2 + x3
par(mfrow = c(2,2), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0))
plot.ts(x1, ylim=c(-10,10), main=expression(omega==6/100~~~A^2==13),  panel.first=grid(lty=1))
plot.ts(x2, ylim=c(-10,10), main=expression(omega==10/100~~~A^2==41),  panel.first=grid(lty=1))
plot.ts(x3, ylim=c(-10,10), main=expression(omega==40/100~~~A^2==85),  panel.first=grid(lty=1))
plot.ts(x,  ylim=c(-16,16), main="sum",  panel.first=grid(lty=1), font.main=1)
dev.off()

##############
# Periodogram
par(mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
P = Mod(2*fft(x)/100)^2;  Fr = 0:99/100
plot(Fr, P, type="o", xlab="frequency", ylab="scaled periodogram", panel.first=grid(lty=1), ylim=c(0,90) )
abline(v=.5, lty=2, col=4)
abline(v=c(.1,.3,.7,.9), lty=1, col='lightgray')
dev.off()