#####################################
### R code for Chapter 4 Examples ###
#####################################

### Example 4.1 (Scatterplot matrix)

data(iris)
#virginica data in first 4 columns of the last 50 obs.

# not shown in text
pairs(iris[101:150, 1:4])

panel.d <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, .5))
  lines(density(x))
}

# Fig. 4.1
x <- scale(iris[101:150, 1:4])
r <- range(x)
pairs(x, diag.panel = panel.d, xlim = r, ylim = r)

library(lattice)
splom(iris[101:150, 1:4])    #plot 1

#for all 3 at once, in color, plot 2
splom(iris[,1:4], groups = iris$Species)

# Fig. 4.2
#for all 3 at once, black and white, plot 3
splom(~iris[1:4], groups = Species, data = iris,
      col = 1, pch = c(1, 2, 3),  cex = c(.5,.5,.5))

# same plot
splom(iris[,1:4], groups = iris$Species,
      col = 1, pch = c(1, 2, 3),  cex = c(.5,.5,.5))


### Example 4.2 (Plot bivariate normal density)

#the standard BVN density
f <- function(x,y) {
  z <- (1/(2*pi)) * exp(-.5 * (x^2 + y^2))
}

y <- x <- seq(-3, 3, length= 50)
z <- outer(x, y, f)   #compute density for all (x,y)

persp(x, y, z)        #the default plot

persp(x, y, z, theta = 45, phi = 30, expand = 0.6,
      ltheta = 120, shade = 0.75, ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "f(x, y)")
# More attractive
persp(x, y, z, theta = 45, phi = 30, expand = 0.6,
      ltheta = 120, shade = 0.75, ticktype = "simple",
      xlab = "X", ylab = "Y", zlab = "", box=F, 
      col='lightblue')

### Example 4.3 (Add elements to perspective plot)

#store viewing transformation in M
M <- persp(x, y, z, theta = 45, phi = 30,
           expand = .4, box = FALSE)

#add some points along a circle
a <- seq(-pi, pi, pi/16)
newpts <- cbind(cos(a), sin(a)) * 2
newpts <- cbind(newpts, 0, 1)  #z=0, t=1
N <- newpts %*% M
points(N[,1]/N[,4], N[,2]/N[,4], col=2)

#add lines
x2 <- seq(-3, 3, .1)
y2 <- -x2^2 / 3
z2 <- dnorm(x2) * dnorm(y2)
N <- cbind(x2, y2, z2, 1) %*% M
lines(N[,1]/N[,4], N[,2]/N[,4], col=4)

#add text
x3 <- c(0, 3.1)
y3 <- c(0, -3.1)
z3 <- dnorm(x3) * dnorm(y3) * 1.1
N <- cbind(x3, y3, z3, 1) %*% M
text(N[1,1]/N[1,4], N[1,2]/N[1,4], "f(x,y)")
text(N[2,1]/N[2,4], N[2,2]/N[2,4], bquote(y==-x^2/3))


### Example 4.4 (Surface plot using wireframe(lattice))

library(lattice)
x <- y <- seq(-3, 3, length= 50)

xy <- expand.grid(x, y)
z <- (1/(2*pi)) * exp(-.5 * (xy[,1]^2 + xy[,2]^2))
wireframe(z ~ xy[,1] * xy[,2])


### Example 4.5 (3D scatterplot)

library(lattice)
attach(iris)
#basic 3 color plot with arrows along axes
cloud(Petal.Length ~ Sepal.Length * Sepal.Width,
      data=iris, groups=Species)

print(cloud(Sepal.Length ~ Petal.Length * Petal.Width,
            data = iris, groups = Species, main = "1", pch=1:3,
            scales = list(draw = FALSE), zlab = "SL",
            screen = list(z = 30, x = -75, y = 0)),
      split = c(1, 1, 2, 2), more = TRUE)

print(cloud(Sepal.Width ~ Petal.Length * Petal.Width,
            data = iris, groups = Species, main = "2", pch=1:3,
            scales = list(draw = FALSE), zlab = "SW",
            screen = list(z = 30, x = -75, y = 0)),
      split = c(2, 1, 2, 2), more = TRUE)

print(cloud(Petal.Length ~ Sepal.Length * Sepal.Width,
            data = iris, groups = Species, main = "3", pch=1:3,
            scales = list(draw = FALSE), zlab = "PL",
            screen = list(z = 30, x = -55, y = 0)),
      split = c(1, 2, 2, 2), more = TRUE)

print(cloud(Petal.Width ~ Sepal.Length * Sepal.Width,
            data = iris, groups = Species, main = "4", pch=1:3,
            scales = list(draw = FALSE), zlab = "PW",
            screen = list(z = 30, x = -55, y = 0)),
      split = c(2, 2, 2, 2))
detach(iris)


### Example 4.6 (Contour plot)


#contour plot with labels
contour(volcano, asp = 1, labcex = 1)

#another version from lattice package
library(lattice)
contourplot(volcano) #similar to above


### Example 4.7 (Filled contour plots)

image(volcano, col = terrain.colors(100), axes = FALSE)
contour(volcano, levels = seq(100,200,by = 10), add = TRUE)

filled.contour(volcano, color = terrain.colors, asp = 1)
levelplot(volcano, scales = list(draw = FALSE),
          xlab = "", ylab = "")


### Example 4.8 (2D histogram)

library(hexbin)
x <- matrix(rnorm(4000), 2000, 2)
plot(hexbin(x[,1], x[,2]))


