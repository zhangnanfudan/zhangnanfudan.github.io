### Covariance Functions
library(fields)

#########
# Several Materns of different smoothness with a similar correlation 
# range

# find ranges for nu = .5, 1.0 and 2.0 
# where the correlation drops to .1 at a distance of 10 units.

r1<- Matern.cor.to.range( 10, nu=.5, cor.target=.1)
r2<- Matern.cor.to.range( 10, nu=1.0, cor.target=.1)
r3<- Matern.cor.to.range( 10, nu=2.0, cor.target=.1)
r4<- Matern.cor.to.range( 10, nu=100, cor.target=.1)

# note that these equivalent ranges
# with respect to this correlation length are quite different
# due the different smoothness parameters. 

d<- seq( 0, 15,,200)
y<- cbind(  Matern( d, range=r1, nu=.5),
            Matern( d, range=r2, nu=1.0),
            Matern( d, range=r3, nu=2.0),
            Matern( d, range=r4, nu=100))
matplot( d, y, type="l", lty=1, lwd=2)

#########
#Simulate a Gaussian random field with an exponential covariance function,  
#range parameter = 2.0 and the domain is  [0,5]X [0,5] evaluating the 
#field at a 100X100 grid.  
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100)) 
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
# Now simulate another ... 
look2<- sim.rf( obj)

# Suppose one requires an exponential, range = 2
# but marginal variance = 10 ( rho in fields notation)
look3<- sqrt( 10)*  sim.rf( obj)

# take a look at first two  
set.panel(2,1)
image.plot( grid$x, grid$y, look) 
title("simulated gaussian fields")
image.plot( grid$x, grid$y, look2) 
title("another realization ...")



# multiply 2-d isotropic exponential with theta=4 by a random vector 

junk<- matrix(rnorm(100*100), 100,100)

cov.obj<- stationary.image.cov( setup=TRUE, 
                                grid=list(x=1:100,y=1:100),theta=8) 
result<-  stationary.image.cov(Y=junk,cov.obj=cov.obj)

image( matrix( result, 100,100)) # NOTE that is also a smoother!

# to do it again, no setup is needed 
#  e.g. 
#  junk2<- matrix(rnorm(100**2, 100,100))
#  result2<-  stationary.image.cov(Y=junk2, cov.obj=cov.obj)

# generate a grid and set of indices based on discretizing the locations
# in the precip dataset

out<-as.image( RMprecip$y, x= RMprecip$x)
ind1<- out$ind
grid<- list( x= out$x, y=out$y)

#
# discretized x locations  to use for comparison
xd<- cbind( out$x[ out$ind[,1]], out$y[ out$ind[,2]] )

# setup to create cov.obj for exponential covariance with range= 1.25

cov.obj<- stationary.image.cov( setup=TRUE, grid=grid, theta=1.25) 

# multiply covariance matrix by an arbitrary vector
junk<-  rnorm(nrow( ind1))
result<- stationary.image.cov( ind1, ind1, Y= junk,cov.obj=cov.obj)

# The brute force way would be  
#   result<- stationary.cov( xd, xd, theta=1.25, C=junk)
# or 
#   result<- stationary.cov( xd, xd, theta=1.25) %*% junk
# both of these take much longer 


# evaluate the covariance between all grid points and the center grid point
Y<- matrix(0,cov.obj$m, cov.obj$n)
Y[32,32]<- 1
result<- stationary.image.cov( Y= Y,cov.obj=cov.obj)
# covariance surface with respect to the grid point at (32,32)
# 
# reshape "vector" as an image
temp<-  matrix( result, cov.obj$m,cov.obj$n)
image.plot(cov.obj$grid$x,cov.obj$grid$y, temp)
# or persp( cov.obj$grid$x,cov.obj$grid$y, temp) 

# check out the Matern
grid<- list( x= seq(-105,-99,,64), y= seq( 40,45,,64)) 
cov.obj<- matern.image.cov( 
    setup=TRUE, grid=grid, theta=.55, smoothness=1.0)
Y<- matrix(0,64,64)
Y[16,16]<- 1

result<- matern.image.cov( Y= Y,cov.obj=cov.obj)
temp<-  matrix( result, cov.obj$m,cov.obj$n)
image.plot( cov.obj$grid$x,cov.obj$grid$y, temp)

# Note we have centered at the location (grid$x[16],grid$y[16]) for this case
#  using sim.rf to simulate an Matern field
look<- sim.rf( cov.obj)
image.plot( grid$x, grid$y, look)
