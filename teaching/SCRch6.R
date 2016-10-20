#####################################
### R code for Chapter 6 Examples ###
#####################################


### Example 6.1 (Basic Monte Carlo estimation)

    m <- 100000
    g <- numeric(m)
    for (i in 1:m) {
        x <- rnorm(2)
        g[i] <- abs(x[1] - x[2])
    }
    est <- mean(g)
    est
    
    # unbiased
    sd(g)/sqrt(m)
    # unbiased
    sqrt(sum((g - mean(g))^2)/(m*(m-1))) 
    # biased
    sqrt(sum((g - mean(g))^2)) / m
    # True
    sqrt((2-4/pi)/m)

### Example 6.2 (Estimating the MSE of a trimmed mean)

    n <- 20
    m <- 1000
    tmean <- numeric(m)
    for (i in 1:m) {
        x <- sort(rnorm(n))
        tmean[i] <- sum(x[2:(n-1)]) / (n-2)
        }
    mse <- mean(tmean^2)
    mse
    sqrt(sum((tmean - mean(tmean))^2)) / m    #se

    n <- 20
    m <- 1000
    tmean <- numeric(m)
    for (i in 1:m) {
        x <- sort(rnorm(n))
        tmean[i] <- median(x)
        }
    mse <- mean(tmean^2)
    mse
    sqrt(sum((tmean - mean(tmean))^2)) / m    #se

    # Different levels of trimming
    n <- 20
    m <- 1000
    mse=numeric(10)
    
    for(k in 0:9)
    {
      tmean <- numeric(m) 
      for (i in 1:m) 
      {
        x <- sort(rnorm(n))
        tmean[i] <- sum(x[(k+1):(n-k)]) / (n-2*k) 
      }
      mse[k+1] <- mean(tmean^2)
    }
    
    plot(mse)
    
### Example 6.3 (MSE of a trimmed mean, cont.)

     set.seed(522)
     n <- 20
     K <- n/2 - 1
     m <- 1000
     mse <- matrix(0, n/2, 6)

     trimmed.mse <- function(n, m, k, p) {
         #MC est of mse for k-level trimmed mean of
         #contaminated normal pN(0,1) + (1-p)N(0,100)
         tmean <- numeric(m)
         for (i in 1:m) {
             sigma <- sample(c(1, 10), size = n,
                 replace = TRUE, prob = c(p, 1-p))
             x <- sort(rnorm(n, 0, sigma))
             tmean[i] <- sum(x[(k+1):(n-k)]) / (n-2*k)
             }
         mse.est <- mean(tmean^2)
         se.mse <- sqrt(mean((tmean-mean(tmean))^2)) / sqrt(m)
         return(c(mse.est, se.mse))
     }

    for (k in 0:K) {
        mse[k+1, 1:2] <- trimmed.mse(n=n, m=m, k=k, p=1.0)
        mse[k+1, 3:4] <- trimmed.mse(n=n, m=m, k=k, p=.95)
        mse[k+1, 5:6] <- trimmed.mse(n=n, m=m, k=k, p=.9)
    }

### Example 6.4 (Confidence interval for variance)

    n <- 20
    alpha <- .05
    x <- rnorm(n, mean=0, sd=2)
    UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)
    UCL


### Example 6.5 (MC estimate of confidence level)
    
    #set.seed(123)
    n <- 20
    alpha <- .05
    UCL <- replicate(1000, expr = {
        x <- rnorm(n, mean = 0, sd = 2)
        (n-1) * var(x) / qchisq(alpha, df = n-1)
        } )
    #count the number of intervals that contain sigma^2=4
    sum(UCL > 4)
    #or compute the mean to get the confidence level
    mean(UCL > 4)


### Example 6.6 (Empirical confidence level)

    n <- 20
    alpha <- .05
    UCL <- replicate(1000, expr = {
        x <- rchisq(n, df = 2)
        (n-1) * var(x) / qchisq(alpha, df = n-1)
        } )
    sum(UCL > 4)
    mean(UCL > 4)


### Example 6.7 (Empirical Type I error rate)

    n <- 20
    alpha <- .05
    mu0 <- 500
    sigma <- 100

    m <- 10000          #number of replicates
    p <- numeric(m)     #storage for p-values
    for (j in 1:m) {
        x <- rnorm(n, mu0, sigma)
        ttest <- t.test(x, alternative = "greater", mu = mu0)
        p[j] <- ttest$p.value
        }

    p.hat <- mean(p < alpha)
    se.hat <- sqrt(p.hat * (1 - p.hat) / m)
    print(c(p.hat, se.hat))


### Example 6.8 (Skewness test of normality)

    n <- c(10, 20, 30, 50, 100, 500) #sample sizes
    cv <- qnorm(.975, 0, sqrt(6/n))  #crit. values for each n

    sk <- function(x) {
        #computes the sample skewness coeff.
        xbar <- mean(x)
        m3 <- mean((x - xbar)^3)
        m2 <- mean((x - xbar)^2)
        return( m3 / m2^1.5 )
    }

    #n is a vector of sample sizes
    #we are doing length(n) different simulations

    p.reject <- numeric(length(n)) #to store sim. results
    m <- 10000                     #num. repl. each sim.

    for (i in 1:length(n)) {
        sktests <- numeric(m)       #test decisions
        for (j in 1:m) {
            x <- rnorm(n[i])
            #test decision is 1 (reject) or 0
            sktests[j] <- as.integer(abs(sk(x)) >= cv[i] )
            }
        p.reject[i] <- mean(sktests) #proportion rejected
    }

    p.reject


### Example 6.9 (Empirical power)

    n <- 20
    m <- 1000
    mu0 <- 500
    sigma <- 100
    mu <- c(seq(450, 650, 10))  #alternatives
    M <- length(mu)
    power <- numeric(M)
    for (i in 1:M) {
        mu1 <- mu[i]
        pvalues <- replicate(m, expr = {
            #simulate under alternative mu1
            x <- rnorm(n, mean = mu1, sd = sigma)
            ttest <- t.test(x,
                     alternative = "greater", mu = mu0)
            ttest$p.value  } )
        power[i] <- mean(pvalues <= .05)
    }

    #par(ask = TRUE)
    library(Hmisc)  #for errbar
    plot(mu, power)
    abline(v = mu0, lty = 1)
    abline(h = .05, lty = 1)

    #add standard errors
    se <- sqrt(power * (1-power) / m)
    errbar(mu, power, yplus = power+se, yminus = power-se,
        xlab = bquote(theta))
    lines(mu, power, lty=3)
    detach(package:Hmisc)
    #par(ask = FALSE)


### Example 6.10 (Power of the skewness test of normality)

    alpha <- .1
    n <- 30
    m <- 2500
    epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
    N <- length(epsilon)
    pwr <- numeric(N)
    #critical value for the skewness test
    cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

    for (j in 1:N) {           #for each epsilon
        e <- epsilon[j]
        sktests <- numeric(m)
        for (i in 1:m) {       #for each replicate
            sigma <- sample(c(1, 10), replace = TRUE,
                size = n, prob = c(1-e, e))
            x <- rnorm(n, 0, sigma)
            sktests[i] <- as.integer(abs(sk(x)) >= cv)
            }
        pwr[j] <- mean(sktests)
        }
    #plot power vs epsilon
    plot(epsilon, pwr, type = "b",
         xlab = bquote(epsilon), ylim = c(0,1))
    abline(h = .1, lty = 3)
    se <- sqrt(pwr * (1-pwr) / m)  #add standard errors
    lines(epsilon, pwr+se, lty = 3)
    lines(epsilon, pwr-se, lty = 3)

