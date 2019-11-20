#####################
### The Jackknife ###
#####################
### Example 7.6 (Jackknife estimate of bias)

data(patch, package = "bootstrap")
n <- nrow(patch)
y <- patch$y
z <- patch$z
theta.hat <- mean(y) / mean(z)
print (theta.hat)

#compute the jackknife replicates, leave-one-out estimates
theta.jack <- numeric(n)
for (i in 1:n)
  theta.jack[i] <- mean(y[-i]) / mean(z[-i])
bias <- (n - 1) * (mean(theta.jack) - theta.hat)

print(bias)  #jackknife estimate of bias


### Example 7.7 (Jackknife estimate of standard error)

se <- sqrt((n-1) *
             mean((theta.jack - mean(theta.jack))^2))
print(se)


### Example 7.8 (Failure of jackknife)
#for the specific example given
set.seed(123)
#change the seed to see other examples

n <- 10
x <- sample(1:100, size = n)

#jackknife estimate of se
M <- numeric(n)
for (i in 1:n) {        #leave one out
  y <- x[-i]
  M[i] <- median(y)
}

Mbar <- mean(M)
print(sqrt((n-1)/n * sum((M - Mbar)^2)))

#bootstrap estimate of se
Mb <- replicate(1000, expr = {
  y <- sample(x, size = n, replace = TRUE)
  median(y) })
print(sd(Mb)) 
print(x)
print(M)
print(Mb)
