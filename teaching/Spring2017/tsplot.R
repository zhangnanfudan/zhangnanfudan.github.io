# http://www.stat.pitt.edu/stoffer/tsa4/R_toot.htm

library(astsa)

data()
jj # not a vector
class(jj)

jj[1]
jj[5]

jjm = as.matrix(jj)
dim(jjm)  

# make a monthly time series object that starts in June of the year 2293
zardoz = ts(rnorm(48), start=c(2293,6), frequency=12)
zardoz

# use window() if you want a part of a ts object
oz = window(zardoz, start=2294, end=c(2295,12))
oz

time(jj)   
cycle(jj)

plot(jj, ylab="Earnings per Share", main="J & J")
plot(jj, type="o", col="blue", lty="dashed")
plot(diff(log(jj)), main="logged and diffed") 

x = -5:5                  # sequence of integers from -5 to 5
y = 5*cos(x)              # guess
#---  plot:
plot(x, main="plot(x)")
plot(x, y, main="plot(x,y)")
#---  plot.ts:
plot.ts(x, main="plot.ts(x)")
plot.ts(x, y, main="plot.ts(x,y)")
#---  ts.plot:
ts.plot(x, main="ts.plot(x)")
# note- x and y are ts objects 
ts.plot(ts(x), ts(y), col=1:2, main="ts.plot(x,y)")  

dev.off()











