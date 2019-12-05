library(mlbench)
library(energy)
alpha <- .1
n <- c(20, 40, 60, 80, 100)
m <- 100 #try small m for a trial run
d <- c(2, 5, 8)
test.ku <- test.sk <- test.eng <- matrix(NA, length(n), length(d))

# function to compute kurtosis statistic
ku <- function(x){
    b <- 0
    for (i in 1:n){
        b <- b + ((x[i,]-colMeans(x)) %*% solve(cov(x)) %*% t(t(x[i,]-colMeans(x))))^2
    }
    b <- b/n
    return (b)
}

# function to compute skewness statistic
sk <- function(x){
    b <- 0
    for (i in 1:n){
        for (j in 1:n){
            b <- b + ((x[i,]-colMeans(x)) %*% solve(cov(x)) %*% t(t(x[j,]-colMeans(x))))^3
        }
    }
    b <- b/(n*n)
    return (b)
}

#critical value for the skewness and kurtosis test
cvku <- qnorm(1-alpha/2)
cvsk <- function(d){qchisq(1-alpha,df=d*(d+1)*(d+2)/6)}


# estimate power
for(nn in 1:length(n)){
    for(dd in 1:length(d)){
        my.n=n[nn]; my.d=d[dd]
        my.ku=my.sk=my.eng=NA
        for (i in 1:m){
            x <- mlbench.twonorm(my.n, my.d)$x
            my.sk[i] <- as.integer(abs(my.n*sk(x)/6)>=cvsk(my.d))
            my.ku[i] <- as.integer(abs((ku(x)-my.d*(my.d+2))/sqrt(8*my.d*(my.d+2)/my.n))>=cvku)
            my.eng[i] <- as.integer(mvnorm.etest(x,R=200)$p.value<=alpha)
        }
        test.sk[nn,dd] <- mean(my.sk)
        test.ku[nn,dd] <- mean(my.ku)
        test.eng[nn,dd] <- mean(my.eng)
    }
}

# my.d=5; my.n=100

# list(mean(my.sk), mean(my.ku), mean(my.eng))

par(mfrow=c(1,3))
for(dd in 1:length(d)){
    plot(c(1,100), c(0,1), xlab='sample size', ylab='power', main=paste('d=',d[dd]), type='n')
    lines(n, test.sk[,dd], lwd=2, col=1)
    lines(n, test.ku[,dd], lwd=2, col=2)
    lines(n, test.eng[,dd], lwd=2, col=3)
}
