### Example 9.7 (Gibbs sampler: Bivariate distribution)
set.seed(123)
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
     ylab=bquote(X[2]), ylim=range(x[,2]), xlim=range(x[,1]), main='MCMC')
plot(x_MASS, cex=.5, xlab=bquote(X[1]),
     ylab=bquote(X[2]), ylim=range(x[,2]), xlim=range(x[,1]), main='MASS')



### Example 9.8 (Gelman-Rubin method of monitoring convergence)

Gelman.Rubin <- function(psi) {
    # psi[i,j] is the statistic psi(X[i,1:j])
    # for chain in i-th row of X
    psi <- as.matrix(psi)
    n <- ncol(psi)
    k <- nrow(psi)

    psi.means <- rowMeans(psi)     #row means
    B <- n * var(psi.means)        #between variance est.
    psi.w <- apply(psi, 1, "var")  #within variances
    W <- mean(psi.w)               #within est.
    v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
    r.hat <- v.hat / W             #G-R statistic
    return(r.hat)
    }

normal.chain <- function(sigma, N, X1) {
    #generates a Metropolis chain for Normal(0,1)
    #with Normal(X[t], sigma) proposal distribution
    #and starting value X1
    x <- rep(0, N)
    x[1] <- X1
    u <- runif(N)

    for (i in 2:N) {
        xt <- x[i-1]
        y <- rnorm(1, xt, sigma)     #candidate point
        r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
        r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
        r <- r1 / r2
        if (u[i] <= r) x[i] <- y else
             x[i] <- xt
        }
    return(x)
    }

sigma <- 4     #parameter of proposal distribution
k <- 4          #number of chains to generate
n <- 15000      #length of chains
b <- 1000       #burn-in length

#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
    X[i, ] <- normal.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
    psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
set.seed(123)
par(mfrow=c(2,2))
for (i in 1:k)
    plot(psi[i, (b+1):n], type="l",
        xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
set.seed(123)
rhat <- rep(0, n)
for (j in (b+1):n)
    rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)


### Example 9.9 (Coal mining disasters)

library(boot)     #for coal data
data(coal)
year <- floor(coal)
y <- table(year)
plot(y)  #a time plot

y <- floor(coal[[1]])  # convert data.frame to numeric
y <- tabulate(y)
y <- y[1851:length(y)]   # starting from year 1851

# Gibbs sampler for the coal mining change point

# initialization
n <- length(y)    #length of the data
m <- 1000         #length of the chain
mu <- lambda <- k <- numeric(m)
L <- numeric(n)
k[1] <- sample(1:n, 1)
mu[1] <- 1
lambda[1] <- 1
b1 <- 1
b2 <- 1

# run the Gibbs sampler
for (i in 2:m) {
    kt <- k[i-1]

    #generate mu
    r <- .5 + sum(y[1:kt])
    mu[i] <- rgamma(1, shape = r, rate = kt + b1)

    #generate lambda
    if (kt + 1 > n) r <- .5 + sum(y) else
        r <- .5 + sum(y[(kt+1):n])
    lambda[i] <- rgamma(1, shape = r, rate = n - kt + b2)

    #generate b1 and b2
    b1 <- rgamma(1, shape = .5, rate = mu[i]+1)
    b2 <- rgamma(1, shape = .5, rate = lambda[i]+1)

    for (j in 1:n) {
        L[j] <- exp((lambda[i] - mu[i]) * j) *
                      (mu[i] / lambda[i])^sum(y[1:j])
        }
    L <- L / sum(L)

    #generate k from discrete distribution L on 1:n
    k[i] <- sample(1:n, prob=L, size=1)
}

b <- 201
j <- k[b:m]

print(mean(k[b:m]))
print(mean(lambda[b:m]))
print(mean(mu[b:m]))


### Code for Figure 9.11

# plots of the chains for Gibbs sampler output

par(mfcol=c(3,1))
plot(mu, type="l", ylab="mu")
plot(lambda, type="l", ylab="lambda")
plot(k, type="l", ylab="change point = k")


### Code for Figure 9.12

# histograms from the Gibbs sampler output

par(mfrow=c(1,3))
labelk <- "changepoint"
label1 <- paste("mu", round(mean(mu[b:m]), 1))
label2 <- paste("lambda", round(mean(lambda[b:m]), 1))

hist(mu[b:m], main="", xlab=label1,
     breaks = "scott", prob=TRUE) #mu posterior
hist(lambda[b:m], main="", xlab=label2,
     breaks = "scott", prob=TRUE) #lambda posterior
hist(j, breaks=min(j):max(j), prob=TRUE, main="",
    xlab = labelk)

dev.off()
