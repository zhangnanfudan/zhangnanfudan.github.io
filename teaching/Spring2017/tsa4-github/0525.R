library('astsa')
source('grid.r')


##############
# Spectral densities
par(mfrow=c(3,1), mar=c(3,3,1.5,1), mgp=c(1.6,.6,0), cex.main=1.1)
arma.spec(log="no", main="White Noise")
arma.spec(ma=.5, log="no", main="Moving Average")
arma.spec(ar=c(1,-.9), log="no", main="Autoregression")
dev.off()