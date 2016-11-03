#####################################
### R code for Chapter 7 Examples ###
#####################################

### Example 7.1
set.seed(8)
x=rpois(10,2)
hist(x)


### Example 7.2 (Bootstrap estimate of standard error)

library(bootstrap)    #for the law data
print(cor(law$LSAT, law$GPA))
print(cor(law82$LSAT, law82$GPA))

#set up the bootstrap
B <- 200            #number of replicates
n <- nrow(law)      #sample size
R <- numeric(B)     #storage for replicates

#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  LSAT <- law$LSAT[i]       #i is a vector of indices
  GPA <- law$GPA[i]
  R[b] <- cor(LSAT, GPA)
}
#output
print(se.R <- sd(R))
hist(R, prob = TRUE)


### Example 7.3 (Bootstrap estimate of standard error: boot function)

r <- function(x, i) {
  #want correlation of columns 1 and 2
  cor(x[i,1], x[i,2])
}

library(boot)       #for boot function
obj <- boot(data = law, statistic = r, R = 2000)
obj
y <- obj$t
sd(y)


### Example 7.4 (Bootstrap estimate of bias)

#sample estimate for n=15
theta.hat <- cor(law$LSAT, law$GPA)

#bootstrap estimate of bias
B <- 2000   #larger for estimating bias
n <- nrow(law)
theta.b <- numeric(B)

for (b in 1:B) {
  i <- sample(1:n, size = n, replace = TRUE)
  LSAT <- law$LSAT[i]
  GPA <- law$GPA[i]
  theta.b[b] <- cor(LSAT, GPA)
}
bias <- mean(theta.b - theta.hat)
bias
# compare
obj


### Example 7.5 (Bootstrap estimate of bias of a ratio estimate)


data(patch, package = "bootstrap")
patch

n <- nrow(patch)  #in bootstrap package
B <- 2000
theta.b <- numeric(B)
theta.hat <- mean(patch$y) / mean(patch$z)

#bootstrap
for (b in 1:B) {
  i <- sample(1:n, size = n, replace = TRUE)
  y <- patch$y[i]
  z <- patch$z[i]
  theta.b[b] <- mean(y) / mean(z)
}
bias <- mean(theta.b) - theta.hat
se <- sd(theta.b)
print(list(est=theta.hat, bias = bias,
           se = se, cv = bias/se))


#################
### Example 7.10 (Bootstrap confidence intervals for patch ratio statistic)

library(boot)       #for boot and boot.ci
data(patch, package = "bootstrap")

theta.boot <- function(dat, ind) {
  #function to compute the statistic
  y <- dat[ind, 1]
  z <- dat[ind, 2]
  mean(y) / mean(z)
}

y <- patch$y
z <- patch$z
dat <- cbind(y, z)
boot.obj <- boot(dat, statistic = theta.boot, R = 2000)

print(boot.obj)
print(boot.ci(boot.obj,
              type = c("basic", "norm", "perc")))


#calculations for bootstrap confidence intervals
alpha <- c(.025, .975)

#normal
print(boot.obj$t0 + qnorm(alpha) * sd(boot.obj$t))

#basic
print(2*boot.obj$t0 - 
        quantile(boot.obj$t, rev(alpha), type=1))

#percentile
print(quantile(boot.obj$t, alpha, type=6))


### Example 7.11 (Bootstrap confidence intervals for the correlation statistic)

library(boot)
data(law, package = "bootstrap")
boot.obj <- boot(law, R = 2000,
                 statistic = function(x, i){cor(x[i,1], x[i,2])})
print(boot.ci(boot.obj, type=c("basic","norm","perc")))


### Example 7.12 (Bootstrap t confidence interval)

boot.t.ci <-
  function(x, B = 500, R = 100, level = .95, statistic){
    #compute the bootstrap t CI
    x <- as.matrix(x);  n <- nrow(x)
    stat <- numeric(B); se <- numeric(B)
    
    boot.se <- function(x, R, f) {
      #local function to compute the bootstrap
      #estimate of standard error for statistic f(x)
      x <- as.matrix(x); m <- nrow(x)
      th <- replicate(R, expr = {
        i <- sample(1:m, size = m, replace = TRUE)
        f(x[i, ])
      })
      return(sd(th))
    }
    
    for (b in 1:B) {
      j <- sample(1:n, size = n, replace = TRUE)
      y <- x[j, ]
      stat[b] <- statistic(y)
      se[b] <- boot.se(y, R = R, f = statistic)
    }
    stat0 <- statistic(x)
    t.stats <- (stat - stat0) / se
    se0 <- sd(stat)
    alpha <- 1 - level
    Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
    names(Qt) <- rev(names(Qt))
    CI <- rev(stat0 - Qt * se0)
  }


### Example 7.13 (Bootstrap t confidence interval for patch ratio statistic)

#boot package and patch data were loaded in Example 7.10
library(boot)       #for boot and boot.ci
data(patch, package = "bootstrap")

dat <- cbind(patch$y, patch$z)
stat <- function(dat) {
  mean(dat[, 1]) / mean(dat[, 2]) }
ci <- boot.t.ci(dat, statistic = stat, B=2000, R=200)
print(ci)

