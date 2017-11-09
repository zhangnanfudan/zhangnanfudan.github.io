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
    sd(g)/sqrt(m) # sd uses m-1 like var
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
    sqrt(sum((tmean - mean(tmean))^2)) / m    #standard error

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
    
    mse # see errata
    
    
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
    
    