### Example 5.1 (Simple Monte Carlo integration)
m <- 10000
x <- runif(m)
theta.hat <- mean(exp(-x))
print(theta.hat)
print(1 - exp(-1))


### Example 5.2 (Simple Monte Carlo integration, cont.)
m <- 10000
x <- runif(m, min=2, max=4)
theta.hat <- mean(exp(-x)) * 2
print(theta.hat)
print(exp(-2) - exp(-4))


### Example 5.3 (Monte Carlo integration, unbounded interval)
x <- seq(.1, 2.5, length = 10)
m <- 10000
u <- runif(m)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
    g <- x[i] * exp(-(u * x[i])^2 / 2)
    cdf[i] <- mean(g) / sqrt(2 * pi) + 0.5
}

Phi <- pnorm(x)
print(round(rbind(x, cdf, Phi), 3))


### Example 5.4 (Example 5.3, cont.)
x <- seq(.1, 2.5, length = 10)
m <- 10000
z <- rnorm(m)
dim(x) <- length(x)
p <- apply(x, MARGIN = 1,
         FUN = function(x, z) {mean(z < x)}, z = z)

Phi <- pnorm(x)
print(round(rbind(x, p, Phi), 3))

### Example 5.5 (Error bounds for MC integration)
x <- 2
m <- 10000
z <- rnorm(m)
g <- (z < x)  #the indicator function
v <- mean((g - mean(g))^2) / m
cdf <- mean(g)
c(cdf, v)
c(cdf - 1.96 * sqrt(v), cdf + 1.96 * sqrt(v))
my.cdf = pnorm(2)
my.v = my.cdf*(1-my.cdf)/m
c(my.cdf, my.v)
c(my.cdf - 1.96 * sqrt(my.v), my.cdf + 1.96 * sqrt(my.v))

