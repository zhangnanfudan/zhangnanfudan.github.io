### Proj.6D
library(tidyverse)
library(mlbench)
library(energy)

skew.test <- function(X, alpha = 0.1) {
  # multivariate skewness test
  # return 1 if reject, 0 if not reject
  n <- nrow(X)
  d <- ncol(X)
  X.bar <- colMeans(X)
  
  b1d <- function(X) {
    result <- 0
    Sigma.inv <- solve(cov(X))
    for (i in 1:n) {
      for (j in 1:n) {
        result <- result + (matrix(X[i,]-X.bar, nrow = 1) %*%
                              Sigma.inv %*% matrix(X[j,]-X.bar, ncol = 1))^3
      }
    }
    return(result/n^2)
  }
  statistic <- n*b1d(X)/6
  return(as.integer(statistic >= qchisq(1-alpha, d*(d+1)*(d+2)/6)))
}

kurtosis.test <- function(X, alpha = 0.1) {
  # multivariate kurtosis test
  # return 1 if reject, 0 if not reject
  n <- nrow(X)
  d <- ncol(X)
  X.bar <- colMeans(X)
  
  b2d <- function(X) {
    result <- 0
    Sigma.inv <- solve(cov(X))
    for (i in 1:n) {
      result <- result + (matrix(X[i,]-X.bar, nrow = 1) %*%
                            Sigma.inv %*% matrix(X[i,]-X.bar, ncol = 1))^2
    }
    return(result / n)
  }
  
  statistic <- abs((b2d(X) - d*(d+2)) / sqrt(8*d*(d+2)/n))
  return(as.integer(statistic >= qnorm(1-0.5*alpha))) 
}

energy.test <- function(X, alpha = 0.1) {
  # energy test
  # return 1 if reject, 0 if not reject
  p <- mvnorm.etest(X, R=200)$p.value
  return(as.integer(p <= alpha))
}

# calculate empirical power for the three tests
m <- 100 # number of trials
alpha <- 0.1


# trial 1) 不同dimension的mlbench.twonorm()数据
emp.power1 <- tibble(dim = numeric(),
                     skew = numeric(),
                     kurt = numeric(),
                     energy = numeric())

skew.test.rejects <- kurtosis.test.rejects <- energy.test.rejects <- numeric(m)

for (d in c(2,4,6,8)) {
  for (i in 1:m) {
    X <- mlbench.twonorm(n = 30, d = d)$x
    skew.test.rejects[i] <- skew.test(X, alpha = alpha)
    kurtosis.test.rejects[i] <- kurtosis.test(X, alpha = alpha)
    energy.test.rejects[i] <- energy.test(X, alpha = alpha)
  }
  emp.power1 <- add_row(emp.power1,
                        dim = d, skew = mean(skew.test.rejects),
                        kurt = mean(kurtosis.test.rejects),
                        energy = mean(energy.test.rejects))
}

# trial 2) 不同样本量的mlbench.twonorm()数据
emp.power2 <- tibble(num = numeric(),
                     skew = numeric(),
                     kurt = numeric(),
                     energy = numeric())
skew.test.rejects <- kurtosis.test.rejects <- energy.test.rejects <- numeric(m)
for (n in c(20, 40, 60, 80, 100)) {
  for (i in 1:m) {
    X <- mlbench.twonorm(n, d = 5)$x
    skew.test.rejects[i] <- skew.test(X, alpha = alpha)
    kurtosis.test.rejects[i] <- kurtosis.test(X, alpha = alpha)
    energy.test.rejects[i] <- energy.test(X, alpha = alpha)
  }
  emp.power2 <- add_row(emp.power2,
                        num = n, skew = mean(skew.test.rejects),
                        kurt = mean(kurtosis.test.rejects),
                        energy = mean(energy.test.rejects))
}
