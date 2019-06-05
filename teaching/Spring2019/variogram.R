library(fields)
data=read.table('data1.txt',header=T)
nlon=length(unique(data[,1]))
nlat=length(unique(data[,2]))
loc=data[,1:2]
y=data[,3]
varigram=vgram(loc, y, lon.lat=TRUE,N=15)

par(mfrow=c(1,2))
##here obtained is semivarigram,then times 2 to get varigram##
##Plot the varigram cloud 
plot(varigram$d,2*varigram$vgram,xlab="distance",ylab="varigram",main="Variogram cloud")
lines(varigram$centers, varigram$stats["mean",], col=4,lwd=3)


bplot.xy(varigram$d, sqrt(2*varigram$vgram),ylab="sqrt(VG)", N=15, outlier=TRUE, main="Boxplots of variogram")
lines(varigram$centers, sqrt(2*varigram$stats["mean",]), col=4,lwd=3)

graphics.off()


###directional variogram, fix the latitude.
data1=data[101:150,]
data2=data[201:250,]
data3=data[301:350,]
data4=data[401:450,]

vgram_lat1=vgram(data1[,1:2], data1[,3], lon.lat=TRUE, N=10)
vgram_lat2=vgram(data2[,1:2], data2[,3], lon.lat=TRUE, N=10)
vgram_lat3=vgram(data3[,1:2], data3[,3], lon.lat=TRUE, N=10)
vgram_lat4=vgram(data4[,1:2], data4[,3], lon.lat=TRUE, N=10)

# pdf('lat.pdf',width=14,height=14)
par(mfrow=c(2,2))
bplot.xy(vgram_lat1$d, sqrt(2*vgram_lat1$vgram),ylab="sqrt(VG)", N=10, outlier=TRUE, main="latitude -49.72")
lines(vgram_lat1$centers, sqrt(2*vgram_lat1$stats["mean",]), col=4, lwd=3)

bplot.xy(vgram_lat2$d, sqrt(2*vgram_lat2$vgram),ylab="sqrt(VG)", N=10, outlier=TRUE, main="latitude -9.94")
lines(vgram_lat2$centers, sqrt(2*vgram_lat2$stats["mean",]), col=4, lwd=3)

bplot.xy(vgram_lat3$d, sqrt(2*vgram_lat3$vgram),ylab="sqrt(VG)", N=10, outlier=TRUE, main="latitude 29.83")
lines(vgram_lat3$centers, sqrt(2*vgram_lat3$stats["mean",]), col=4, lwd=3)

bplot.xy(vgram_lat4$d, sqrt(2*vgram_lat4$vgram),ylab="sqrt(VG)", N=10, outlier=TRUE, main="latitude 69.61")
lines(vgram_lat4$centers, sqrt(2*vgram_lat4$stats["mean",]), col=4, lwd=3)

graphics.off()

# 
# ###use variog4 to plot variogram in 4 directions
# library(geoR)
# data2=as.list(data)
# var4=variog4(coords =loc,data =y)
# 
# plot(var4,lwd=3)
# graphics.off()