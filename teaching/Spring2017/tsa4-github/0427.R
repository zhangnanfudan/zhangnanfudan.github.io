library('astsa')
source('grid.r')

#######################
# Estimation for Recruitment Series
### Yule-Walker 
rec.yw = ar.yw(rec, order=2)
rec.yw
rec.yw$x.mean # = 62.26 (mean estimate)
rec.yw$ar # parameter estimates
sqrt(diag(rec.yw$asy.var.coef)) # standard errors
rec.yw$var.pred # error variance estimate

rec.pr = predict(rec.yw, n.ahead=24)
U = rec.pr$pred + rec.pr$se
L = rec.pr$pred - rec.pr$se
minx = min(rec,L); maxx = max(rec,U)
par(mar=c(2.5,2.5,0,0)+.5, mgp=c(1.6,.6,0))
ts.plot(rec, rec.pr$pred, xlim=c(1980,1990), ylim=c(minx,maxx))
grid(lty=1); par(new=TRUE)
lines(rec.pr$pred, col="red", type="o")
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")
dev.off()

### MLE 
rec.mle = ar.mle(rec, order=2)
rec.mle$x.mean
rec.mle$ar 
sqrt(diag(rec.mle$asy.var.coef))


#######################
# Gauss-Newton optimatization of glacial varve data
x=diff(log(varve))
r=acf(x, lag=1, plot=FALSE)$acf[-1]
rstart = (1-sqrt(1-4*(r^2)))/(2*r)    #example 3.29 (e2.27)
c(0) -> w 
c() -> Sc 
num = length(x)
th = seq(-.3,-.94,-.01)
for (p in 1:length(th)){
  for (i in 2:num){w[i]=x[i]-th[p]*w[i-1]}
  Sc[p] = sum(w^2)
}		
par(mar=c(2,2.5,0,0)+.5, mgp=c(1.6,.6,0))	
plot(th, Sc, type="l",ylab=expression(S[c](theta)), xlab=expression(theta),lwd=2, panel.first=grid(NA, NULL,lty=1)) 
# estimation
c(0) -> w -> z
c() -> Sc -> Sz -> Szw
para = c()
niter = 15
para[1]=rstart
for (p in 1:niter){
  for (i in 2:num){w[i]=x[i]-para[p]*w[i-1]
  z[i]=w[i-1]-para[p]*z[i-1]
  }
  Sc[p] = sum(w^2)				   
  Sz[p]=sum(z^2)
  Szw[p]=sum(z*w)
  para[p+1] = para[p] + Szw[p]/Sz[p]
}  
#round(cbind(iteration=0:(niter-1), thetahat=para[1:niter] , Sc , Sz ), 3)
abline(v=para[1:12], lty=2)
points(para[1:12], Sc[1:12], pch=16)
dev.off()
