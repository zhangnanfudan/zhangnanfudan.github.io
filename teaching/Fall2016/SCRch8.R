########################
# R code for Chapter 8 #
########################

### Example 8.1 (Permutation distribution of a statistic)

    attach(chickwts)
    x <- sort(as.vector(weight[feed == "soybean"]))
    y <- sort(as.vector(weight[feed == "linseed"]))
    detach(chickwts)

    R <- 999              #number of replicates
    z <- c(x, y)          #pooled sample
    K <- 1:26
    reps <- numeric(R)   #storage for replicates
    t0 <- t.test(x, y)$statistic

    for (i in 1:R) {
        #generate indices k for the first sample
        k <- sample(K, size = 14, replace = FALSE)
        x1 <- z[k]
        y1 <- z[-k]      #complement of x1
        reps[i] <- t.test(x1, y1)$statistic
        }
    p <- mean(c(t0, reps) >= t0)
    p

    hist(reps, main = "", freq = FALSE, xlab = "T (p = 0.202)",
        breaks = "scott")
    points(t0, 0, cex = 1, pch = 16)      #observed T


### Example 8.2 (Permutation distribution of the K-S statistic)    

    # continues Example 8.1
    R <- 999             #number of replicates
    z <- c(x, y)         #pooled sample
    K <- 1:26
    D <- numeric(R)      #storage for replicates
    options(warn = -1)
    D0 <- ks.test(x, y, exact = FALSE)$statistic
    for (i in 1:R) {
        #generate indices k for the first sample
        k <- sample(K, size = 14, replace = FALSE)
        x1 <- z[k]
        y1 <- z[-k]      #complement of x1
        D[i] <- ks.test(x1, y1, exact = FALSE)$statistic
        }
    p <- mean(c(D0, D) >= D0)
    options(warn = 0)
    p

    hist(D, main = "", freq = FALSE, xlab = "D (p = 0.46)",
        breaks = "scott")
    points(D0, 0, cex = 1, pch = 16)      #observed D


### Example 8.3 (Example 8.2, cont.)

    attach(chickwts)
    x <- sort(as.vector(weight[feed == "sunflower"]))
    y <- sort(as.vector(weight[feed == "linseed"]))
    detach(chickwts)
    summary(cbind(x, y))
    options(warn = -1)
    D0 <- ks.test(x, y, exact = FALSE)$statistic
    for (i in 1:R) {
        #generate indices k for the first sample
        k <- sample(K, size = 14, replace = FALSE)
        x1 <- z[k]
        y1 <- z[-k]      #complement of x1
        D[i] <- ks.test(x1, y1, exact = FALSE)$statistic
        }
    p <- mean(c(D0, D) >= D0)
    options(warn = 0)
    p

