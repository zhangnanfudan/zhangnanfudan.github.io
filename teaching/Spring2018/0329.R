library(astsa)

####################
# March 29, 2018
####################
# Example 2.1 Linear trend
summary(fit <- lm(chicken~time(chicken))) # regress price on time
tsplot(chicken, ylab="cents per pound", col=4, lwd=2)
abline(fit)           # add the fitted regression line to the plot          

####################
# Example 2.2 LA pollution, temperature and mortality
par(mfrow=c(3,1))
tsplot(cmort, main="Cardiovascular Mortality", ylab="")
tsplot(tempr, main="Temperature",  ylab="")
tsplot(part, main="Particulates", ylab="")

pairs(cbind(Mortality=cmort, Temperature=tempr, Particulates=part))

## Regression
temp  = tempr-mean(tempr)  # center temperature    
temp2 = temp^2             # square it  
trend = time(cmort)        # time

fit = lm(cmort~ trend + temp + temp2 + part, na.action=NULL)

summary(fit)       # regression results
summary(aov(fit))  # ANOVA table   (compare to next line)
summary(aov(lm(cmort~cbind(trend, temp, temp2, part)))) # Table 2.1

num = length(cmort)                                     # sample size
AIC(fit)/num - log(2*pi)                                # AIC 
BIC(fit)/num - log(2*pi)                                # BIC 
# AIC(fit, k=log(num))/num - log(2*pi)                  # BIC (alt method) 
# (AICc = log(sum(resid(fit)^2)/num) + (num+5)/(num-5-2)) # AICc
dev.off()

####################
# Examples 2.3 Regression with lagged variables
fish = ts.intersect(rec, soiL6=lag(soi,-6), dframe=TRUE)   
summary(fit <- lm(rec~soiL6, data=fish, na.action=NULL))
tsplot(fish$rec, ylim=c(0,111))  # plot the data and the fitted values (not shown in text) 
lines(fitted(fit), col=2)
dev.off()

####################
# Examples 2.4 and 2.5
fit = lm(chicken~time(chicken), na.action=NULL) # regress chicken on time
par(mfrow=c(2,1))
tsplot(resid(fit), main="detrended")
tsplot(diff(chicken), main="first difference")

par(mfrow=c(3,1))     # plot ACFs
acf1(chicken, 48, main="chicken")
acf1(resid(fit), 48, main="detrended")
acf1(diff(chicken), 48, main="first difference")
dev.off()

####################
# Example 2.6
par(mfrow=c(2,1))
tsplot(diff(globtemp), type="o")
mean(diff(globtemp))     # drift estimate = .008
acf1(diff(gtemp), 48, main="")

####################
# Example 2.7
par(mfrow=c(2,1))
tsplot(varve, main="varve", ylab="")
tsplot(log(varve), main="log(varve)", ylab="" )
