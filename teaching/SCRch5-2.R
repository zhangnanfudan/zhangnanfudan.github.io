##########################################
### R code for Chapter 5 Examples      ###
##########################################


### Example 5.7 (Control variate)

    m <- 10000
    a <- - 12 + 6 * (exp(1) - 1)
    U <- runif(m)
    T1 <- exp(U)                  #simple MC
    T2 <- exp(U) + a * (U - 1/2)  #controlled

    mean(T1)    
    mean(T2)    
    (var(T1) - var(T2)) / var(T1)


### Example 5.8 (MC integration using control variates)

    f <- function(u)
        exp(-.5)/(1+u^2)

    g <- function(u)
        exp(-u)/(1+u^2)

    set.seed(510) #needed later
    u <- runif(10000)
    B <- f(u)
    A <- g(u)

    cor(A, B)
    a <- -cov(A,B) / var(B)    #est of c*
    a

    m <- 100000
    u <- runif(m)
    T1 <- g(u)
    T2 <- T1 + a * (f(u) - exp(-.5)*pi/4)

    c(mean(T1), mean(T2))
    c(var(T1), var(T2))
    (var(T1) - var(T2)) / var(T1)


### Example 5.9 (Control variate and regression)

    set.seed(510)
    u <- runif(10000)
    f <- exp(-.5)/(1+u^2)
    g <- exp(-u)/(1+u^2)
    c.star <-  - lm(g ~ f)$coeff[2]   # beta[1]
    mu <- exp(-.5)*pi/4

    c.star

    u <- runif(10000)
    f <- exp(-.5)/(1+u^2)
    g <- exp(-u)/(1+u^2)
    L <- lm(g ~ f)
    theta.hat <- sum(L$coeff * c(1, mu))  #pred. value at mu

    theta.hat
    summary(L)$sigma^2
    summary(L)$r.squared


### Example 5.10 (Choice of the importance function)
    #code for plot is at the end of the file

    m <- 10000
    theta.hat <- se <- numeric(5)
    g <- function(x) {
        exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
        }

    x <- runif(m)     #using f0
    fg <- g(x)
    theta.hat[1] <- mean(fg)
    se[1] <- sd(fg)

    x <- rexp(m, 1)   #using f1
    fg <- g(x) / exp(-x)
    theta.hat[2] <- mean(fg)
    se[2] <- sd(fg)

    x <- rcauchy(m)   #using f2
    i <- c(which(x > 1), which(x < 0))
    x[i] <- 2  #to catch overflow errors in g(x)
    fg <- g(x) / dcauchy(x)
    theta.hat[3] <- mean(fg)
    se[3] <- sd(fg)

    u <- runif(m)     #f3, inverse transform method
    x <- - log(1 - u * (1 - exp(-1)))
    fg <- g(x) / (exp(-x) / (1 - exp(-1)))
    theta.hat[4] <- mean(fg)
    se[4] <- sd(fg)

    u <- runif(m)    #f4, inverse transform method
    x <- tan(pi * u / 4)
    fg <- g(x) / (4 / ((1 + x^2) * pi))
    theta.hat[5] <- mean(fg)
    se[5] <- sd(fg)

    rbind(theta.hat, se)
    

### Example 5.11 (Example 5.10, cont.)

    M <- 20   #number of replicates
    T2 <- numeric(4)
    estimates <- matrix(0, 10, 2)

    g <- function(x) {
        exp(-x - log(1+x^2)) * (x > 0) * (x < 1) }

    for (i in 1:10) {
        estimates[i, 1] <- mean(g(runif(M)))
        T2[1] <- mean(g(runif(M/4, 0, .25)))
        T2[2] <- mean(g(runif(M/4, .25, .5)))
        T2[3] <- mean(g(runif(M/4, .5, .75)))
        T2[4] <- mean(g(runif(M/4, .75, 1)))
        estimates[i, 2] <- mean(T2)
    }
    
    estimates
    apply(estimates, 2, mean)
    apply(estimates, 2, var)


### Example 5.12 (Examples 5.10-5.11, cont.)

    M <- 10000  #number of replicates
    k <- 10     #number of strata
    r <- M / k  #replicates per stratum
    N <- 50     #number of times to repeat the estimation
    T2 <- numeric(k)
    estimates <- matrix(0, N, 2)

    g <- function(x) {
        exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
        }

    for (i in 1:N) {
        estimates[i, 1] <- mean(g(runif(M)))
        for (j in 1:k)
            T2[j] <- mean(g(runif(M/k, (j-1)/k, j/k)))
        estimates[i, 2] <- mean(T2)
    }

    apply(estimates, 2, mean)
    apply(estimates, 2, var)
    


### Plot importance functions in Figures 5.1(a) and 5.1.(b)

    #par(ask = TRUE) #uncomment to pause between graphs
    
    x <- seq(0, 1, .01)
    w <- 2
    f1 <- exp(-x)
    f2 <- (1 / pi) / (1 + x^2)
    f3 <- exp(-x) / (1 - exp(-1))
    f4 <- 4 / ((1 + x^2) * pi)
    g <- exp(-x) / (1 + x^2)

    #for color change lty to col

    #figure (a)
    plot(x, g, type = "l", main = "", ylab = "",
         ylim = c(0,2), lwd = w)
    lines(x, g/g, lty = 2, lwd = w)
    lines(x, f1, lty = 3, lwd = w)
    lines(x, f2, lty = 4, lwd = w)
    lines(x, f3, lty = 5, lwd = w)
    lines(x, f4, lty = 6, lwd = w)
    legend("topright", legend = c("g", 0:4),
           lty = 1:6, lwd = w, inset = 0.02)

    #figure (b)
    plot(x, g, type = "l", main = "", ylab = "",
        ylim = c(0,3.2), lwd = w, lty = 2)
    lines(x, g/f1, lty = 3, lwd = w)
    lines(x, g/f2, lty = 4, lwd = w)
    lines(x, g/f3, lty = 5, lwd = w)
    lines(x, g/f4, lty = 6, lwd = w)
    legend("topright", legend = c(0:4),
           lty = 2:6, lwd = w, inset = 0.02)
