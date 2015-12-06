###########2222222222(1)############
rm(list = ls())

h<-matrix(0,nrow=7,ncol=20)
n<-c(25,50,100,200,400,800,1600)

for(i in 1:7)
{
for(j in 1:20)
{
q<-seq(1,n[i],by=1)
x<-(q-0.5)/n[i]
z<-rnorm(n[i])
y<-6*x+1+z
ss<-smooth.spline(x,y,all.knots=T,cv=F)
h[i,j]<-ss$df
}
}

plot(log(n)/log(10),apply(h, 1, median),type='l',ylim=c(0,10),ylab='Median of estimated df',main='estimated df vs sample size')
lines(log(n)/log(10),apply(h, 1, median)+4/(20)^0.5*apply(h, 1,mad))
lines(log(n)/log(10),apply(h, 1, median)-4/(20)^0.5*apply(h, 1, mad))

############22222222222222(2)###############
rm(list = ls())
h<-matrix(0,nrow=7,ncol=20)
n<-c(25,50,100,200,400,800,1600)

for(i in 1:7)
{
for(j in 1:20)
{
q<-seq(1,n[i],by=1)
x<-(q-0.5)/n[i]
z<-rnorm(n[i])
y<-8*exp(-20*(x-0.3)^2)+5*exp(-40*(x-0.8)^2)+z
ss<-smooth.spline(x,y,all.knots=T,cv=F)
h[i,j]<-ss$df
}
}

plot(log(n)/log(10),apply(h, 1, median),type='l',ylim=c(4,16),ylab='Median of estimated df',main='estimated df vs sample size')
lines(log(n)/log(10),apply(h, 1, median)+4/(20)^0.5*apply(h, 1,mad))
lines(log(n)/log(10),apply(h, 1, median)-4/(20)^0.5*apply(h, 1, mad))


#######3(gcv)#################################
library(gss)
rm(list = ls())
h<-matrix(0,nrow=20,ncol=100)
n<-100
sigma2<-rep(0,20)
for(i in 1:20)
{
q<-seq(1,n,by=1)
x<-(q-0.5)/n
z<-rnorm(n)
f0<-1+3*sin(2*pi*x-pi)
y<-f0+z
ss<-ssanova(y~x,type='cubic',method='v')
h[i,]<-predict(ss,data.frame(x=x),se=FALSE)
sigma2[i]<-(summary(ss)$sigma)^2
}

f<-matrix(rep(f0,20),nrow=20,byrow=TRUE)

R<-diag((h-f)%*%t(h-f))/100

####################3(gml)###############3

h<-matrix(0,nrow=20,ncol=100)
n<-100
sigma2m<-rep(0,20)
for(i in 1:20)
{
q<-seq(1,n,by=1)
x<-(q-0.5)/n
z<-rnorm(n)
f0<-1+3*sin(2*pi*x-pi)
y<-f0+z
ss<-ssanova(y~x,type='cubic',method='m')
h[i,]<-predict(ss,data.frame(x=x),se=FALSE)
sigma2m[i]<-(summary(ss)$sigma)^2
}

f<-matrix(rep(f0,20),nrow=20,byrow=TRUE)

RM<-diag((h-f)%*%t(h-f))/100



par(mfrow = c(2, 2))
boxplot(R,ylim=c(0.01,0.20), main='R for GCV')
boxplot(RM,ylim=c(0.01,0.20), main='R for GML')
boxplot(sigma2,ylim=c(0.7,1.5), main='sigma^2 for GCV')
boxplot(sigma2m,ylim=c(0.7,1.5), main='sigma^2 for GML')
summary(R)
summary(RM)
summary(sigma2)
summary(sigma2m)

##########
> summary(R)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01843 0.04543 0.06754 0.06771 0.08135 0.13380 
> summary(RM)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01910 0.03237 0.05708 0.06673 0.08772 0.19970 
> summary(sigma2)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7701  0.8848  0.9764  0.9869  1.0690  1.2940 
> summary(sigma2m)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7244  0.7878  1.0010  0.9762  1.0890  1.2620 
#################4444444444444444(just try with original data)########
rm(list = ls())
library(MASS)
library(gss)
attach(mcycle)
z<-mcycle
for(i in 2:length(z$times))
{
if(z$times[i]==z$times[i-1]){z$accel[i]<-z$accel[i-1]}
}
uz<-unique(z)

ss<-ssanova(uz$accel~uz$times,type='cubic',method='v')
xgrid = seq(min(times),max(times),length=100)
pre<-predict(ss,data.frame(times=xgrid),se=TRUE)

plot(times,accel,xlab="Time (ms)",ylab="Acceleration  (g)",cex=.5,main='Acceleration vs Time')

# plot cubic smoothing splines for 2 to 12 df

lines(xgrid,pre$fit)
lines(xgrid,pre$fit+1.96*pre$se,lty=2,col=3,lwd=2)
lines(xgrid,pre$fit-1.96*pre$se,lty=2,col=3,lwd=2)
t<-summary(ss)
plot(residuals(ss)/(t$sigma*sqre(rep(1,length(xgrid))-t$)))


#######smooth spline(cleaned data)####
rm(list = ls())
library(MASS)
library(gss)
attach(mcycle)
z<-mcycle
for(i in 2:length(z$times))
{
if(z$times[i]==z$times[i-1]){z$accel[i]<-z$accel[i-1]}
}
uz<-unique(z)
ss<-smooth.spline(uz$times, uz$accel,all.knots=T)
sigma=sqrt(ss$pen.crit/(100-ss$df))
se=sigma*sqrt(ss$lev)


plot(times,accel,xlab="Time (ms)",ylab="Acceleration  (g)",cex=.5,main='Acceleration vs Time')

# plot cubic smoothing splines for 2 to 12 df
a<-c(predict(ss,data.frame(times=uz$times))$y$times)
lines(uz$times,a)
lines(uz$times,a+1.96*se,lty=2,col=3,lwd=2)
lines(uz$times,a-1.96*se,lty=2,col=3,lwd=2)


resid<-residuals(ss)
plot(resid/(sigma*sqrt(rep(1,length(uz$times))-ss$lev)),ylab='standard residuals',main='standard residule plot')
abline(h=2)
abline(h=-2)
abline(h=0)

t(a-uz$accel)%*%(a-uz$accel)/length(uz$times)

#######gss(weighted)(cleaned data)###############
library(gss)
plot(uz$times,log(resid^2))
log.var=ssanova(log(resid^2)~uz$times)
varhat=exp(fitted(log.var))
x<-uz$times
y<-uz$accel
weighted.fit=ssanova(y~x, weight=1/varhat)

plot(times,accel,xlab="Time (ms)",ylab="Acceleration  (g)",cex=.5,main='Acceleration vs Time(weighted)')

# plot cubic smoothing splines for 2 to 12 df
est<-predict(weighted.fit,data.frame(x=x),se=T)

lines(x,est$fit)
lines(x,est$fit+1.96*est$se.fit,lty=2,col=3,lwd=2)
lines(x,est$fit-1.96*est$se.fit,lty=2,col=3,lwd=2)

###############gss(cleaned data)##########
fit=ssanova(y~x)

plot(times,accel,xlab="Time (ms)",ylab="Acceleration  (g)",cex=.5,main='Acceleration vs Time')

# plot cubic smoothing splines for 2 to 12 df
est<-predict(fit,data.frame(x=x),se=T)

lines(x,est$fit)
lines(x,est$fit+1.96*est$se.fit,lty=2,col=3,lwd=2)
lines(x,est$fit-1.96*est$se.fit,lty=2,col=3,lwd=2)


resid<-residuals(fit)
plot(x,resid/est$se.fit,ylab='standard residuals',main='standard residule plot')
abline(h=2)
abline(h=-2)
abline(h=0)
























