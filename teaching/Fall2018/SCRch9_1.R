### Example 9.1 (Metropolis-Hastings sampler)

f <- function(x, sigma) {
    if (any(x < 0)) return (0)
    stopifnot(sigma > 0)
    return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

m <- 10000
sigma <- 4
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)

for (i in 2:m) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den) x[i] <- y else {
         x[i] <- xt
         k <- k+1     #y is rejected
         }
    }

print(k)

index <- 5000:5500
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")


### Example 9.2 (Example 9.1, cont.)

b <- 2001      #discard the burnin sample
y <- x[b:m]
a <- ppoints(100)
QR <- sigma * sqrt(-2 * log(1 - a))  #quantiles of Rayleigh
Q <- quantile(x, a)

qqplot(QR, Q, main="",
    xlab="Rayleigh Quantiles", ylab="Sample Quantiles")
abline(0,1, col='red')
hist(y, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR, f(QR, 4))


### Example 9.3 (Random walk Metropolis) 

rw.Metropolis <- function(n, sigma, x0, N) {
    x <- numeric(N)
    x[1] <- x0
    u <- runif(N)
    k <- 0
    for (i in 2:N) {
        y <- rnorm(1, x[i-1], sigma)
            if (u[i] <= (dt(y, n) / dt(x[i-1], n)))
            x[i] <- y  else {
                x[i] <- x[i-1]
                k <- k + 1
            }
        }
    return(list(x=x, k=k))
    }


n <- 4  #degrees of freedom for target Student t dist.
N <- 2000
sigma <- c(.05, .5, 2,  16)

x0 <- 25
rw1 <- rw.Metropolis(n, sigma[1], x0, N)
rw2 <- rw.Metropolis(n, sigma[2], x0, N)
rw3 <- rw.Metropolis(n, sigma[3], x0, N)
rw4 <- rw.Metropolis(n, sigma[4], x0, N)

#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k))
## rejection rate
print(c(rw1$k, rw2$k, rw3$k, rw4$k)/N)

# paths
## individual
par(mfrow=c(2,2))
plot(rw1$x, type='l')
plot(rw2$x, type='l')
plot(rw3$x, type='l')
plot(rw4$x, type='l')

## comparative
y.lim=range(c(rw1$x, rw2$x, rw3$x, rw4$x))
plot(rw1$x, type='l', ylim=y.lim)
plot(rw2$x, type='l', ylim=y.lim)
plot(rw3$x, type='l', ylim=y.lim)
plot(rw4$x, type='l', ylim=y.lim)

dev.off()

# ### Code for Figure 9.3
# 
#     par(mfrow=c(2,2))  #display 4 graphs together
#     refline <- qt(c(.025, .975), df=n)
#     rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
#     for (j in 1:4) {
#         plot(rw[,j], type="l",
#              xlab=bquote(sigma == .(round(sigma[j],3))),
#              ylab="X", ylim=range(rw[,j]))
#         abline(h=refline)
#     }
#     par(mfrow=c(1,1)) #reset to default


### Example 9.4 (Example 9.3, cont.)

a <- c(.05, seq(.1, .9, .1), .95)
Q <- qt(a, n)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
mc <- rw[501:N, ]
Qrw <- apply(mc, 2, function(x) quantile(x, a))
print(round(cbind(Q, Qrw), 3))

plot(a,Q, type='l',lwd=2)
for(i in 1:4){
  lines(a,Qrw[,i],col=i+1)
}
legend('topleft',legend=1:4, lty=1, col=2:5)
# xtable::xtable(round(cbind(Q, Qrw), 3)) #latex format


### Example 9.6 (Independence sampler)

m <- 5000 #length of chain
xt <- numeric(m)
a <- 1          #parameter of Beta(a,b) proposal dist.
b <- 1          #parameter of Beta(a,b) proposal dist.

# alternative
a <- 5            #parameter of Beta(a,b) proposal dist.
b <- 2            #parameter of Beta(a,b) proposal dist.

p <- .2           #mixing parameter
n <- 30           #sample size
mu <- c(0, 5)     #parameters of the normal densities
sigma <- c(1, 1)

# generate the observed sample
i <- sample(1:2, size=n, replace=TRUE, prob=c(p, 1-p))
x <- rnorm(n, mu[i], sigma[i])

# generate the independence sampler chain
u <- runif(m)
y <- rbeta(m, a, b)      #proposal distribution
xt[1] <- .5

for (i in 2:m) {
  fy <- y[i] * dnorm(x, mu[1], sigma[1]) +
    (1-y[i]) * dnorm(x, mu[2], sigma[2])
  fx <- xt[i-1] * dnorm(x, mu[1], sigma[1]) +
    (1-xt[i-1]) * dnorm(x, mu[2], sigma[2])
  
  r <- prod(fy / fx) *
    (xt[i-1]^(a-1) * (1-xt[i-1])^(b-1)) /
    (y[i]^(a-1) * (1-y[i])^(b-1))
  
  if (u[i] <= r) xt[i] <- y[i] else
    xt[i] <- xt[i-1]
}

plot(xt, type="l", ylab="p")
hist(xt[101:m], main="", xlab="p", prob=TRUE)
print(mean(xt[101:m]))


### Example 9.7 (Gibbs sampler: Bivariate distribution)

#initialize constants and parameters
N <- 5000               #length of chain
burn <- 1000            #burn-in length
X <- matrix(0, N, 2)    #the chain, a bivariate sample

rho <- -.75             #correlation
mu1 <- 0
mu2 <- 2
sigma1 <- 1
sigma2 <- .5
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

###### generate the chain #####

X[1, ] <- c(mu1, mu2)            #initialize

for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ]

# compare sample statistics to parameters
colMeans(x)
cov(x)
cor(x)

library(MASS)
S = matrix(c(sigma1^2, rho*sigma1*sigma2,
             rho*sigma1*sigma2, sigma2^2),2)
x_MASS = mvrnorm(n=dim(x)[1],mu=c(mu1, mu2), Sigma=S)

par(mfrow=c(1,2))
plot(x, cex=.5, xlab=bquote(X[1]),
     ylab=bquote(X[2]), ylim=range(x[,2]), main='MCMC')
plot(x_MASS, cex=.5, xlab=bquote(X[1]),
     ylab=bquote(X[2]), ylim=range(x[,2]), main='MASS')





