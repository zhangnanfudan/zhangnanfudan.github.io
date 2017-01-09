#########################
# R code for Chapter 11 #
#########################

### Example 11.1 (Identical and nearly equal)
    .Machine
    
    (.2-.1)==.1
    (0.3-0.1)==0.2
    0.2-(0.3-0.1)

    isTRUE(all.equal(.2, .3 - .1))
    all.equal(.2, .3)          #not a logical value
    isTRUE(all.equal(.2, .3))  #always a logical value

    x <- 1:4
    y <- 2
    y == 2
    x == y  #not necessarily a single logical value
    identical(x, y)  #always a single logical value
    identical(y, 2)


### Example 11.2 (Ratio of two large numbers)

    n <- 400
    factorial(n)
    (gamma((n-1)/2) / (sqrt(pi) * gamma((n-2)/2)))
    exp(lgamma((n-1)/2) - lgamma((n-2)/2)) / sqrt(pi)


### Example 11.3 (Taylor expansion)

    system.time({
        for (i in 1:1000) {
            a <- rep(0, 24)
            a0 <- pi / 6
            a2 <- a0 * a0
            a[1] <- -a0^3 / 6
            for (i in 2:24)
                a[i] <- - a2 * a[i-1] / ((2*i+1)*(2*i))
            a0 + sum(a)}
        })

    system.time({
        for (i in 1:1000) {
            K <- 2 * (0:24) + 1
            i <- rep(c(1, -1), length=25)
            sum(i * (pi/6)^K / factorial(K))}
        })


    
### Example 11.6 (Solving f(x)=0)

    f <- function(y, a, n) {
        a^2 + y^2 + 2*a*y/(n-1) - (n-2)
    }

    a <- 0.5
    n <- 20
    b0 <- 0
    b1 <- 5*n

    #solve using bisection
    it <- 0
    eps <- .Machine$double.eps^0.25
    r <- seq(b0, b1, length=3)
    y <- c(f(r[1], a, n), f(r[2], a, n), f(r[3], a, n))
    if (y[1] * y[3] > 0)
        stop("f does not have opposite sign at endpoints")

    while(it < 1000 && abs(y[2]) > eps) {
        it <- it + 1
        if (y[1]*y[2] < 0) {
            r[3] <- r[2]
            y[3] <- y[2]
        } else {
            r[1] <- r[2]
            y[1] <- y[2]
        }
        r[2] <- (r[1] + r[3]) / 2
        y[2] <- f(r[2], a=a, n=n)
        print(c(r[1], y[1], y[3]-y[2]))
    }
    it

### Example 11.7 (Solving f(x)=0 with Brent's method: uniroot)

    a <- 0.5
    n <- 20
    out <- uniroot(function(y) {
               a^2 + y^2 + 2*a*y/(n-1) - (n-2) },
               lower = 0, upper = n*5)
    unlist(out)
    uniroot(function(y) {a^2 + y^2 + 2*a*y/(n-1) - (n-2)},
            interval = c(-n*5, 0))$root
            

### Example 11.8 (Numerical integration with integrate)

    f <- function(y, N, r, rho) {
        (cosh(y) - rho * r)^(1 - N)
    }
    integrate(f, lower=0, upper=Inf,
              rel.tol=.Machine$double.eps^0.25,
              N=10, r=0.5, rho=0.2)

    ro <- seq(-.99, .99, .01)
    v <- rep(0, length(ro))
    for (i in 1:length(ro)) {
        v[i] <- integrate(f, lower=0, upper=Inf,
                  rel.tol=.Machine$double.eps^0.25,
                  N=10, r=0.5, rho=ro[i])$value
        }
    plot(ro, v, type="l", xlab=expression(rho),
         ylab="Integral Value (n=10, r=0.5)")


### Example 11.9 (Density of sample correlation coefficient)

    .dcorr <- function(r, N, rho=0) {
        # compute the density function of sample correlation
        if (abs(r) > 1 || abs(rho) > 1) return (0)
        if (N < 4) return (NA)

        if (isTRUE(all.equal(rho, 0.0))) {
            a <- exp(lgamma((N - 1)/2) - lgamma((N - 2)/2)) /
                     sqrt(pi)
            return (a * (1 - r^2)^((N - 4)/2))
        }

        # if rho not 0, need to integrate
        f <- function(w, R, N, rho)
            (cosh(w) - rho * R)^(1 - N)

        #need to insert some error checking here
        i <- integrate(f, lower=0, upper=Inf,
                R=r, N=N, rho=rho)$value
        c1 <- (N - 2) * (1 - rho^2)^((N - 1)/2)
        c2 <- (1 - r^2)^((N - 4) / 2) / pi
        return(c1 * c2 * i)
    }

    r <- as.matrix(seq(-1, 1, .01))
    d1 <- apply(r, 1, .dcorr, N=10, rho=.0)
    d2 <- apply(r, 1, .dcorr, N=10, rho=.5)
    d3 <- apply(r, 1, .dcorr, N=10, rho=-.5)
    plot(r, d2, type="l", lty=2, lwd=2, ylab="density")
    lines(r, d1, lwd=2)
    lines(r, d3, lty=4, lwd=2)
    legend("top", inset=.02,
           c("rho = 0", "rho = 0.5", "rho = -0.5"), lty=c(1,2,4), lwd=2)


 ### Example 11.10 (MLE using mle)

    #the observed sample
    y <- c(0.04304550, 0.50263474)

    mlogL <- function(theta=1) {
        #minus log-likelihood of exp. density, rate 1/theta
        return( - (length(y) * log(theta) - theta * sum(y)))
    }

    library(stats4)
    fit <- mle(mlogL)
    summary(fit)
    1/mean(y)

    # Alternately, the initial value for the optimizer could 
    # be supplied in the call to mle; two examples are

    mle(mlogL, start=list(theta=1))
    mle(mlogL, start=list(theta=mean(y)))


### Example 11.11 (One-dimensional optimization with optimize)

    x <- seq(2, 8, .001)
    y <- log(x + log(x))/(log(1+x))
    plot(x, y, type = "l")

    f <- function(x)
        log(x + log(x))/log(1+x)

    optimize(f, lower = 4, upper = 8, maximum = TRUE)


### Example 11.12 (MLE: Gamma distribution)

    m <- 20000
    est <- matrix(0, m, 2)
    n <- 200
    r <- 5
    lambda <- 2

    obj <- function(lambda, xbar, logx.bar) {
        digamma(lambda * xbar) - logx.bar - log(lambda)
        }

    for (i in 1:m) {
        x <- rgamma(n, shape=r, rate=lambda)
        xbar <- mean(x)
        u <- uniroot(obj, lower = .001, upper = 10e5,
                xbar = xbar, logx.bar = mean(log(x)))
        lambda.hat <- u$root
        r.hat <- xbar * lambda.hat
        est[i, ] <- c(r.hat, lambda.hat)
    }

    ML <- colMeans(est)
    ML

    hist(est[, 1], breaks="scott", freq=FALSE,
        xlab="r", main="")
    points(ML[1], 0, cex=1.5, pch=20)
    hist(est[, 2], breaks="scott", freq=FALSE,
         xlab=bquote(lambda), main="")
    points(ML[2], 0, cex=1.5, pch=20)


### Example 11.13 (Two-dimensional optimization with optim)

    LL <- function(theta, sx, slogx, n) {
        r <- theta[1]
        lambda <- theta[2]
        loglik <- n * r * log(lambda) + (r - 1) * slogx -
            lambda * sx - n * log(gamma(r))
        - loglik
        }

    n <- 200
    r <- 5;    lambda <- 2
    x <- rgamma(n, shape=r, rate=lambda)

    optim(c(1,1), LL, sx=sum(x), slogx=sum(log(x)), n=n)

    mlests <- replicate(20000, expr = {
      x <- rgamma(200, shape = 5, rate = 2)
      optim(c(1,1), LL, sx=sum(x), slogx=sum(log(x)), n=n)$par
      })
    colMeans(t(mlests))


### Example 11.16 (Simplex algorithm)

    library(boot)   #for simplex function
    A1 <- rbind(c(-2, 1, 1), c(4, -1, 3))
    b1 <- c(1, 3)
    a <- c(2, 2, 3)
    simplex(a = a, A1 = A1, b1 = b1, maxi = TRUE)
    detach(package:boot)


