#########################################
###   R code for Chapter 9 Examples   ###
#########################################

### Example 9.1 (Metropolis-Hastings sampler)

    f <- function(x, sigma) {
        if (any(x < 0)) return (0)
        stopifnot(sigma > 0)
        return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
    }

    m <- 10000
    sigma <- 4
    x <- numeric(m)
    x[1] <- rchisq(1, df=1)
    k <- 0
    u <- runif(m)

    for (i in 2:m) {
        xt <- x[i-1]
        y <- rchisq(1, df = xt)
        num <- f(y, sigma) * dchisq(xt, df = y)
        den <- f(xt, sigma) * dchisq(y, df = xt)
        if (u[i] <= num/den) x[i] <- y else {
             x[i] <- xt
             k <- k+1     #y is rejected
             }
        }

    print(k)

    index <- 5000:5500
    y1 <- x[index]
    plot(index, y1, type="l", main="", ylab="x")


### Example 9.2 (Example 9.1, cont.)

    b <- 2001      #discard the burnin sample
    y <- x[b:m]
    a <- ppoints(100)
    QR <- sigma * sqrt(-2 * log(1 - a))  #quantiles of Rayleigh
    Q <- quantile(x, a)

    qqplot(QR, Q, main="",
        xlab="Rayleigh Quantiles", ylab="Sample Quantiles")

    hist(y, breaks="scott", main="", xlab="", freq=FALSE)
    lines(QR, f(QR, 4))


### Example 9.3 (Random walk Metropolis) 

    rw.Metropolis <- function(n, sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (dt(y, n) / dt(x[i-1], n)))
                x[i] <- y  else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
        }

    n <- 4  #degrees of freedom for target Student t dist.
    N <- 2000
    sigma <- c(.05, .5, 2,  16)

    x0 <- 25
    rw1 <- rw.Metropolis(n, sigma[1], x0, N)
    rw2 <- rw.Metropolis(n, sigma[2], x0, N)
    rw3 <- rw.Metropolis(n, sigma[3], x0, N)
    rw4 <- rw.Metropolis(n, sigma[4], x0, N)

    #number of candidate points rejected
    print(c(rw1$k, rw2$k, rw3$k, rw4$k))
    

### Code for Figure 9.3

    par(mfrow=c(2,2))  #display 4 graphs together
    refline <- qt(c(.025, .975), df=n)
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
    par(mfrow=c(1,1)) #reset to default


### Example 9.4 (Example 9.3, cont.)

    a <- c(.05, seq(.1, .9, .1), .95)
    Q <- qt(a, n)
    rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
    mc <- rw[501:N, ]
    Qrw <- apply(mc, 2, function(x) quantile(x, a))
    print(round(cbind(Q, Qrw), 3))
    
    plot(a,Q, type='l',lwd=2)
    for(i in 1:4){
      lines(a,Qrw[,i],col=i+1)
    }
    legend('topleft',legend=1:4, lty=1, col=2:5)
    # xtable::xtable(round(cbind(Q, Qrw), 3)) #latex format

### Code for Figures 9.4(a) and 9.4(b)

    plot(x, type="l")
    abline(h=b, v=501, lty=3)
    xb <- x[- (1:501)]
    hist(xb, prob=TRUE, xlab=bquote(beta), ylab="X", main="")
    z <- seq(min(xb), max(xb), length=100)
    lines(z, dnorm(z, mean(xb), sd(xb)))


### Example 9.5 (Bayesian inference: A simple investment model)

    b <- .2          #actual value of beta
    w <- .25         #width of the uniform support set
    m <- 5000        #length of the chain
    burn <- 1000     #burn-in time
    days <- 250
    x <- numeric(m)  #the chain

    # generate the observed frequencies of winners
    i <- sample(1:5, size=days, replace=TRUE,
            prob=c(1, 1-b, 1-2*b, 2*b, b))
    win <- tabulate(i)
    print(win)

    prob <- function(y, win) {
        # computes (without the constant) the target density
        if (y < 0 || y >= 0.5)
            return (0)
        return((1/3)^win[1] *
            ((1-y)/3)^win[2] * ((1-2*y)/3)^win[3] *
                ((2*y)/3)^win[4] * (y/3)^win[5])
    }

    u <- runif(m)         #for accept/reject step
    v <- runif(m, -w, w)  #proposal distribution
    x[1] <- .25
    for (i in 2:m) {
        y <- x[i-1] + v[i]
        if (u[i] <= prob(y, win) / prob(x[i-1], win))
            x[i] <- y  else
                x[i] <- x[i-1]
    }

    print(win)
    print(round(win/days, 3))
    print(round(c(1, 1-b, 1-2*b, 2*b, b)/3, 3))
    xb <- x[(burn+1):m]
    print(mean(xb))
    
