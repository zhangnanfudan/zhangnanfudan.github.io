library(geoR)

#
# computing variograms:
#
# binned variogram
vario.b <- variog(s100, max.dist=1)
# variogram cloud
vario.c <- variog(s100, max.dist=1, op="cloud")
#binned variogram and stores the cloud
vario.bc <- variog(s100, max.dist=1, bin.cloud=TRUE)
# smoothed variogram
vario.s <- variog(s100, max.dist=1, op="sm", band=0.2)
#
#
# plotting the variograms:
par(mfrow=c(2,2))
plot(vario.b, main="binned variogram") 
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")  
plot(vario.s, main="smoothed variogram") 
dev.off()

# computing a directional variogram
var4 <- variog4(s100, max.dist=1)
plot(var4)
