library(MCMCpack)
library(maps)
library(MASS)
library(vars)
library(mvtnorm)
#####################################################################
##      Get the range of theta1 by plotting the target function    ##
#####################################################################

###Set the fixed parameter
p=3  #vector length
beta0=5
beta1=0.25 #beta1^2
tau=1   #tau^2
sigma=2.25 #sigma^2
one=t(t(array(1,p)))  #(1,...,1)t
identity=diag(p)
epsilon=0.1
X=c(2,4,9)
Xbar=mean(X)
theta=array(NA,p)
lambda=NA

i=1 #for theta1


#how we generate data X
#lambda=rnorm(1,beta0,beta1)
#theta=mvrnorm(1,lambda*one,tau*identity)
#X=mvrnorm(1,theta,sigma*identity)

###calculate a0 and a, as well as ABCD


a=(1-epsilon)/epsilon*exp(-1/2*(sum((X-beta0*one)^2)-
p^2*beta1/(sigma+tau+p*beta1))/(sigma+tau))

a0=a*(tau/(tau+sigma)*X[i]+beta0/(tau+p*beta1)+
p*sigma*beta1*Xbar/((sigma+tau)*(sigma+tau+p*beta1)))


A=1/((2*pi*(sigma+tau))^((p-1)/2)*sqrt(p))*exp(-sum((X-Xbar*one)^2)/2/(sigma
+tau))
B=(tau*X[i]+sigma*Xbar)/(sigma+tau)
C=B*A
D=sigma*A/sqrt(p*(sigma+tau))

###Get the function T,u,v and Hg,H0
T=function(z)
{T=-dnorm((beta0+z-Xbar)/sqrt((sigma+tau)/p))+
dnorm((beta0-Xbar)/sqrt((sigma+tau)/p))
return(T)}

u=function(z)
{u=pnorm((beta0+z-Xbar)/sqrt((sigma+tau)/p))-
pnorm((beta0-Xbar)/sqrt((sigma+tau)/p))
return(u)}

v=function(z)
{v=dnorm((beta0+z-Xbar)/sqrt((sigma+tau)/p))
return(v)}


Hg=function(z)
{if(z!=0) Hg=C*u(z)/z+D*T(z)/z
 if(z==0) Hg=1/(2*pi*(tau+sigma))^(p/2)*exp(-sum((X-beta0)^2)/2/(sigma+tau))*
(tau*X[i]+sigma*beta0)/(sigma+tau)
return(Hg)}

H0=function(z)
{if(z!=0) H0=A*u(z)/z
 if(z==0) H0=1/(2*pi*(tau+sigma))^(p/2)*exp(-sum((X-beta0)^2)/2/(sigma+tau))
return(H0)}


target=function(z) 
{target=(a0+Hg(z))/(a+H0(z))
return(target)}


z=seq(0.01,20,by=0.01)
targety=c(target(0),target(z))
targetx=c(0,z)
plot(targetx,targety,main="function(z) versus z")



#look for z at minimum and maximum value of target function
target(0)
min(target(z)) 
z[which.min(target(z))]

max(target(z))
z[which.max(target(z))]

#In this case, target(z) ranges from 4.337593 to 4.337607. 
#where we get the minimum at z=0.72, and maximum at z=Inf. And the maximum
#value of the function equals the posterior mean of theta given pi0.

########just some check code
#posterior mean of theta under pi0.
post0.theta=tau/(tau+sigma)*X+beta0/(tau+p*beta1)*one+
p*sigma*beta1*Xbar*one/((tau+sigma)*(tau+sigma+p*beta1))
#check m(x|pi0)---small
dmvnorm(X,beta0*one,(tau+sigma)*identity+beta1*matrix(1,p,p),log=FALSE)



#####################################################################
##       Try to calculate z value according to my derivation       ##
#####################################################################
length=100
zz=array(NA,length)
zz[1]=0.1#0.1 or 0.2

for(i in 2:length){
numerator=D*T(zz[i-1])*(a+A*v(zz[i-1])/sqrt((sigma+tau)/p))+D*v(zz[i-1])*Xbar*(a*zz[i-1]
+A*u(zz[i-1]))/((tau+sigma)/p)-u(zz[i-1])*(a0*A-a*C+A*D*v(zz[i-1])*beta0/
((sigma+tau)/p))-C*a*zz[i-1]*v(zz[i-1])/sqrt((sigma+tau)/p)

denominator=v(zz[i-1])/((tau+sigma)/p)*(a*D*zz[i-1]+a*D*beta0+A*D*u(zz[i-1])
-A*a0*sqrt((tau+sigma)/p))

zz[i]=numerator/denominator }
plot(zz,type="l",col="blue",main="Find z using my derivation")



#####################################################################
##    Try to calculate z value according newton-raphson method     ##
#####################################################################

f=function(z)
{f=(a0+C*v(z)/sqrt((tau+sigma)/p)+D*(beta0+z-Xbar)*v(z)/((sigma+tau)/p))*
(a*z+A*u(z))-(a0*z+C*u(z)+D*T(z))*(a+A*v(z)/sqrt((tau+sigma)/p))
return (f)
}

fdev=function(z)
{fdev=(-C/sqrt((tau+sigma)/p)*(beta0+z-Xbar)*v(z)/((sigma+tau)/p)-D*
(beta0+z-Xbar)^2*v(z)/((sigma+tau)/p)^2+D*v(z)/((sigma+tau)/p))*(a*z+A*u(z))+
(a0+C*v(z)/sqrt((tau+sigma)/p)+D*(beta0+z-Xbar)*v(z)/((sigma+tau)/p))*
(a+A*v(z)/sqrt((sigma+tau)/p))-
(a0+C*v(z)/sqrt((tau+sigma)/p)+D*(beta0+z-Xbar)/((sigma+tau)/p)*v(z))*
(a+A*v(z)/sqrt((tau+sigma)/p))-
(a0*z+C*u(z)+D*T(z))*(-A/sqrt((tau+sigma)/p)*v(z)*(beta0+z-Xbar)/((sigma+tau)/p))

return(fdev)
}

z=seq(0.01,4,by=0.01)
par(mfrow=c(2,2))
plot(z,f(z),type="l")
abline(h=0, v=0,col="red")
#dev.new()
plot(z,fdev(z),type="l")
abline(h=0, v=0,col="red")

#dev.new()
plot(z[1:200],f(z)[1:200]/fdev(z)[1:200])
abline(h=0, v=0,col="red")

plot(z,f(z)/fdev(z))
abline(h=0, v=0,col="red")



length=100
zzz=c(NA,length)
zzz[1]=1
for(i in 2:length){
zzz[i]=zzz[i-1]-f(zzz[i-1])/fdev(zzz[i-1])}

dev.new()
plot(zzz,type="l",col="blue",main="Find z using Newton-Raphson")



######################################################################
#Check my derivation directly, they are right
haha=function(z){
numer=D*T(z)*(a+A*v(z)/sqrt((sigma+tau)/p))+D*v(z)*Xbar*(a*z+A*u(z))/((tau+sigma)/p)-u(z)*(a0*A-a*C+A*D*v(z)*beta0/
((sigma+tau)/p))-C*a*z*v(z)/sqrt((sigma+tau)/p)
denom=v(z)/((tau+sigma)/p)*(a*D*z+a*D*beta0+A*D*u(z)-A*a0*sqrt((tau+sigma)/p))
haha=numer/denom
return(haha)
}

#haha(1.075)  #RIGHT!!!



######################################################################
##Try to calculate z value according newton-raphson method ORIGINAL!##
######################################################################
x = function(z){
x = a0*z+C*u(z)+D*T(z)
return(x)}

y = function(z){
y = a*z+A*u(z)
return(y)}

x1st = function(z){
x1st = a0+C*v(z)/((sigma+tau)/p)^(1/2)+D*v(z)*(beta0+z-Xbar)/((sigma+tau)/p)
return(x1st)}

x2nd = function(z){
x2nd = -C*v(z)*(beta0+z-Xbar)/((sigma+tau)/p)^(3/2)+D*v(z)/((sigma+tau)/p)-
       D*v(z)*(beta0+z-Xbar)^2/((sigma+tau)/p)^2
return(x2nd)}

y1st = function(z){
y1st = a+A*v(z)/((sigma+tau)/p)^(1/2)
return(y1st)}

y2nd = function(z){
y2nd = -A*v(z)*(beta0+z-Xbar)/((sigma+tau)/p)^(3/2)
return(y2nd)}

F = function(z){
F = (x1st(z)*y(z)-x(z)*y1st(z))/(y(z)^2)
return (F)}

Fdev = function(z){
Fdev = x2nd(z)/y(z)-(x(z)*y2nd(z)-2*x1st(z)*y1st(z))/(y(z)^2)+
2*x(z)*(y1st(z))^2/(y(z)^3) 
return (Fdev)}

z=seq(0.01,5,by=0.01)
par(mfrow=c(2,2))
plot(z,F(z),type="l")
abline(h=0, v=0,col="red")
#dev.new()
plot(z,Fdev(z),type="l")
abline(h=0, v=0,col="red")
#dev.new()
plot(z,F(z)/Fdev(z))
abline(h=0, v=0,col="red")

F(0.40)/Fdev(0.40)
F(0.41)/Fdev(0.41)
F(0.42)/Fdev(0.42)


length=1000
zaa=c(NA,length)
zaa[1]=0.415
for(i in 2:length){
zaa[i]=zaa[i-1]-F(zaa[i-1])/Fdev(zaa[i-1])}

par(mfrow=c(2,1))
plot(zaa,type="l",col="blue",main="Find z using ORIGINAL Newton-Raphson starting from 0.415")

zaa=c(NA,length)
zaa[1]=1.420
for(i in 2:length){
zaa[i]=zaa[i-1]-F(zaa[i-1])/Fdev(zaa[i-1])}

plot(zaa,type="l",col="blue",main="Find z using ORIGINAL Newton-Raphson starting from 0.420")





