library(fields)
#
# compute variogram for the midwest ozone field  day 16
# (BTW this looks a bit strange!)
#
data( ozone2)
good<- !is.na(ozone2$y[16,])
x<- ozone2$lon.lat[good,] 
y<- ozone2$y[16,good]

look<-vgram( x,y, N=15, lon.lat=TRUE) # locations are in lon/lat so use right
#distance
# take a look:
plot(look, pch=19)
#lines(look$centers, look$stats["mean",], col=4)

brk<- seq( 0, 250,, (25 + 1) ) # will give 25 bins.

## or some boxplot bin summaries

boxplotVGram(look, breaks=brk, plot.args=list(type="o"))
plot(look, add=TRUE, breaks=brk, col=4)

#############################
#############################
#############################
# Kriging
# a 2-d example 
# fitting a surface to ozone  
# measurements. Exponential covariance, range parameter is 20 (in miles) 

fit <- Krig(ChicagoO3$x, ChicagoO3$y, theta=20)  

summary( fit) # summary of fit 
set.panel( 2,2) 
plot(fit) # four diagnostic plots of fit  
set.panel()
surface( fit, type="C") # look at the surface 

# predict at data
predict( fit)

# predict using 7.5 effective degrees of freedom:
predict( fit, df=7.5)


# predict on a grid ( grid chosen here by defaults)
out<- predictSurface( fit)
surface( out, type="C") # option "C" our favorite

# predict at arbitrary points (10,-10) and (20, 15)
xnew<- rbind( c( 10, -10), c( 20, 15))
predict( fit, xnew)

# standard errors of prediction based on covariance model.  
predictSE( fit, xnew)

# surface of standard errors on a default grid
predictSurfaceSE( fit)-> out.p # this takes some time!
surface( out.p, type="C")
points( fit$x)

## Not run: 
# Using another stationary covariance. 
# smoothness is the shape parameter for the Matern. 

fit <- Krig(ChicagoO3$x, ChicagoO3$y, 
            Covariance="Matern", theta=10, smoothness=1.0)  
summary( fit)

#
# Roll your own: creating very simple user defined Gaussian covariance 
#

test.cov <- function(x1,x2,theta,marginal=FALSE,C=NA){
  # return marginal variance
  if( marginal) { return(rep( 1, nrow( x1)))}
  
  # find cross covariance matrix     
  temp<- exp(-(rdist(x1,x2)/theta)**2)
  if( is.na(C[1])){
    return( temp)}
  else{
    return( temp%*%C)}
} 
#
# use this and put in quadratic polynomial fixed function 


fit.flame<- Krig(flame$x, flame$y, cov.function="test.cov", m=3, theta=.5)

#
# note how range parameter is passed to Krig.   
# BTW:  GCV indicates an interpolating model (nugget variance is zero) 
# This is the content of the warning message.

# take a look ...
surface(fit.flame, type="I") 

## End(Not run)

# 
# Thin plate spline fit to ozone data using the radial 
# basis function as a generalized covariance function 
#
# p=2 is the power in the radial basis function (with a log term added for 
# even dimensions)
# If m is the degree of derivative in penalty then p=2m-d 
# where d is the dimension of x. p must be greater than 0. 
#  In the example below p = 2*2 - 2 = 2  
#

out<- Krig( ChicagoO3$x, ChicagoO3$y,cov.function="Rad.cov", 
            m=2,p=2,scale.type="range") 

# See also the Fields function Tps
# out  should be identical to  Tps( ChicagoO3$x, ChicagoO3$y)
# 

# A Knot example

data(ozone2)
y16<- ozone2$y[16,] 

# there are some missing values -- remove them 
good<- !is.na( y16)
y<- y16[good] 
x<- ozone2$lon.lat[ good,]

#
# the knots can be arbitrary but just for fun find them with a space 
# filling design. Here we select  50 from the full set of 147 points
#
xknots<- cover.design( x, 50, num.nn= 75)$design  # select 50 knot points

out<- Krig( x, y, knots=xknots,  cov.function="Exp.cov", theta=300)  
summary( out)
# note that that trA found by GCV is around 17 so 50>17  knots may be a 
# reasonable approximation to the full estimator. 
#
## Not run: 
# the plot 
surface( out, type="C")
US( add=TRUE)
points( x, col=2)
points( xknots, cex=2, pch="O")

## End(Not run)
## Not run: 
## A quick way to deal with too much data if you intend to smooth it anyway
##  Discretize the locations to a grid, then apply Krig 
##  to the discretized locations:
## 
RM.approx<- as.image(RMprecip$y, x=RMprecip$x, nx=20, ny=20)

# take a look:
image.plot( RM.approx)
# discretized data (observations averaged if in the same grid box)
# 336 locations -- down form the  full 806

# convert the image format to locations, obs and weight vectors
yd<- RM.approx$z[RM.approx$ind]
weights<- RM.approx$weights[RM.approx$ind] # takes into account averaging
xd<- RM.approx$xd

obj<- Krig( xd, yd, weights=weights, theta=4)

# compare to the full fit:
# Krig( RMprecip$x, RMprecip$y, theta=4) 

## End(Not run)

