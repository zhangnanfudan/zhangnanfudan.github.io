library(fields)

### ozone data set, day 16
## Quilt plot
data(ozone2)
# plot 16 day of ozone data set
quilt.plot(ozone2$lon.lat, ozone2$y[16,])
US(add=TRUE, col="grey", lwd=2)

# compute variogram for the midwest ozone field 

good<- !is.na(ozone2$y[16,])
x<- ozone2$lon.lat[good,] 
y<- ozone2$y[16,good]

## variogram plots
par(mfrow=c(1,2))
look<-vgram(x,y, N=15, lon.lat=TRUE)
plot(look, pch=19)

## or some boxplot bin summaries
brk<- seq(0, 250,, (25 + 1) ) # will give 25 bins.
boxplotVGram(look, breaks=brk, plot.args=list(type="o"))
plot(look, add=TRUE, breaks=brk, col=4)
dev.off()


#############################
# Kriging
# a 2-d example 
# fitting a surface to ozone  
# measurements. Exponential covariance, range parameter is 20 (in miles) 

fit <- Krig(ChicagoO3$x, ChicagoO3$y, theta=20)  

summary(fit) # summary of fit 
par(mfrow=c(2,2))
plot(fit) # four diagnostic plots of fit  
dev.off()

surface(fit, type="C") # look at the surface 

# predict at data
predict(fit)

# predict using 7.5 effective degrees of freedom:
predict(fit, df=7.5)

# predict on a grid (grid chosen here by defaults)
out<- predictSurface(fit)
surface(out, type="C") # option "C" our favorite

# predict at arbitrary points (10,-10) and (20, 15)
xnew<- rbind(c(10, -10), c(20, 15))
predict(fit, xnew)

# standard errors of prediction based on covariance model.  
predictSE(fit, xnew)

# surface of standard errors on a default grid
predictSurfaceSE(fit)-> out.p # this takes some time!
surface(out.p, type="C")
points(fit$x)

## Using another stationary covariance. 
# smoothness is the shape parameter for the Matern
fit.matern <- Krig(ChicagoO3$x, ChicagoO3$y, 
                   Covariance="Matern", theta=10, smoothness=1.0)  
summary( fit.matern)

set.panel(2,2)
plot(fit.matern)

set.panel(1,2)
surface(fit, type="C") # look at the surface 
surface(fit.matern, type="C") # look at the surface 


