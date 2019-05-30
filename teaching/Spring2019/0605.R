library(MASS)
library(ggplot2)
library(fields)

##############
### Matern covariance function
d<- seq( 0, 15,length.out = 200)
y<- cbind(  Matern( d, nu=.5),
            Matern( d, nu=1.0),
            Matern( d, nu=2.0))
z<- cbind(  Matern( d, range=1),
            Matern( d, range=2),
            Matern( d, range=3))

set.panel(1,2)
matplot( d, y, type="l", lty=1, lwd=2, col=1:3, main='Matern with range 1')
legend('topright', legend=c('nu=0.5', 'nu=1', 'nu=2'), col=1:3, lty=1, lwd=2)

matplot( d, z, type="l", lty=1, lwd=2, col=1:3, main='Matern with smoothness 0.5')
legend('topright', legend=c('range=1', 'range=2', 'range=3'), col=1:3, lty=1, lwd=2)

dev.off()


##############
### 1D GP with Matern
set.seed(123)

d<- seq( 0, 10,length.out =  200)
nu=.5 # 1 1.5
range=1

my.Matern=function(d1,d2){
    Matern(abs(d1-d2), range=range, nu=nu)
}
my.cov=outer (d,d, my.Matern)

gp = mvrnorm(3, mu=rep(0, length(d)), Sigma = my.cov) 
my.gp = data.frame(x=d,y=t(gp))

ggplot(my.gp, aes(x=x))+
    geom_line(aes(y = y.1, colour = "y.1")) + 
    geom_line(aes(y = y.2, colour = "y.2")) + 
    geom_line(aes(y = y.3, colour = "y.3")) +
    ylab('y') + ggtitle(paste("Matern(nu=", nu,")"))
dev.off()
##############
### 2D GP with Matern

grid<- list( x= seq( 0,10,length.out = 200), y= seq(0,10,length.out = 200)) 

obj1<-matern.image.cov( grid=grid, theta=.5, smoothness=.5, setup=TRUE)
look1<- sim.rf( obj1)

obj2<-matern.image.cov( grid=grid, theta=.5, smoothness=1, setup=TRUE)
look2<- sim.rf( obj2)

obj3<-matern.image.cov( grid=grid, theta=.5, smoothness=2, setup=TRUE)
look3<- sim.rf( obj3)


# take a look at first two  
set.panel(1,3)
image.plot( grid$x, grid$y, look1) 
title("Matern, nu=0.5")
image.plot( grid$x, grid$y, look2) 
title("Matern, nu=1")
image.plot( grid$x, grid$y, look3) 
title("Matern, nu=2 ")

dev.off()




