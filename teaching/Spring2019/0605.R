library(MASS)
library(ggplot2)
library(fields)

set.seed(123)

d<- seq( 0, 10,,200)
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

                        