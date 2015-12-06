
######################################
##  2 groups and mean are different ##
######################################
library(mvtnorm)
library(MASS)
q=2   

sigma=1    #sigma^2=1
theta=list()
p=array(10,q) #100 numbers
norm_vec <- function(x) sum(x^2)

b=0.1
c=4*(q-1)     #b>=0

tau=seq(0.01,5.01,by=0.5)
mean=list()
mean[[1]]=t(t(array(1,p[1])))
mean[[2]]=t(t(array(3,p[2])))
mean0=t(t(array(2,p[1])))
LofRisk=length(tau)

K=20000
X=array(NA,dim=c(K,sum(p)))
JS_Sep1=array(NA,dim=c(K,sum(p)))
JS_Sep2=array(NA,dim=c(K,sum(p)))
JS_Comb1=array(NA,dim=c(K,sum(p)))
JS_Comb2=array(NA,dim=c(K,sum(p)))
OUR1=array(NA,dim=c(K,sum(p)))
OUR2=array(NA,dim=c(K,sum(p)))
OURNEW=array(NA,dim=c(K,sum(p)))


L_JS_Sep1=array(NA,K)
L_JS_Sep2=array(NA,K)
L_JS_Comb1=array(NA,K)
L_JS_Comb2=array(NA,K)
L_OUR1=array(NA,K)
L_OUR2=array(NA,K)
L_MLE=array(NA,K)
L_OURNEW=array(NA,K)


R_JS_Sep1=array(NA,LofRisk)
R_JS_Sep2=array(NA,LofRisk)
R_JS_Comb1=array(NA,LofRisk)
R_JS_Comb2=array(NA,LofRisk)
R_OUR1=array(NA,LofRisk)
R_OUR2=array(NA,LofRisk)
R_OURNEW=array(NA,LofRisk)
R_MLE=array(NA,LofRisk)
theta_norm=array(NA,LofRisk)
R_mean=array(NA,dim=c(LofRisk,8))
R_sd=array(NA,dim=c(LofRisk,8))
R_meanSORT=array(NA,dim=c(LofRisk,8))
R_sdSORT=array(NA,dim=c(LofRisk,8))

print('break pint 1')
gamma=0.5

ones=list()
for(i in 1:q){
  ones[[i]]=array(1,p[i])}


oneALL=array(1,sum(p))

for(l in 1:LofRisk){
  print(l)
  
  for(k in 1:K){
    for(i in 1:q){
      u=runif(1)
      if(u<gamma)
        theta[[i]]=mvrnorm(1,mean[[i]],tau[l]*diag(p[i]))
      else
        theta[[i]]=mvrnorm(1,mean0,tau[l]*diag(p[i]))
    }
    for(i in 1:q){
      X[k,(i*p[[i]]-p[[i]]+1):(i*p[[i]])]=mvrnorm(1,theta[[i]],sigma*diag(p[i]))
    }
    for(j in 1:q){
      if((p[j]-2)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])<=1)
        JS_Sep1[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=(1-(p[j]-2)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]))*X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]
      else
        JS_Sep1[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=array(0,p[j])
    }
    for(j in 1:q){
      if((p[j]-3)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])<=1)
        JS_Sep2[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=(1-(p[j]-3)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]]))*(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])+mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]]
      else
        JS_Sep2[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]]
    }
    
    if((sum(p)-2)/norm_vec(X[k,])<=1) 
      JS_Comb1[k,]=(1-(sum(p)-2)/norm_vec(X[k,]))*X[k,]
    if((sum(p)-2)/norm_vec(X[k,])>1)
      JS_Comb1[k,]=array(0,sum(p))
    
    if((sum(p)-3)/norm_vec(X[k,]-mean(X[k,])*unlist(ones))<=1) 
      JS_Comb2[k,]=(1-(sum(p)-3)/norm_vec(X[k,]-mean(X[k,])*unlist(ones)))*(X[k,]-mean(X[k,])*unlist(ones))+mean(X[k,])*unlist(ones)
    if((sum(p)-3)/norm_vec(X[k,]-mean(X[k,])*unlist(ones))>1)
      JS_Comb2[k,]=mean(X[k,])*unlist(ones)
    
    for(j in 1:q){
      if(((p[j]-2)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])+c/(b+norm_vec(X[k,])))<=1)
        OUR1[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=(1-((p[j]-2)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])+c/(b+norm_vec(X[k,]))))*X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]
      else  
        OUR1[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=array(0,p[[j]])
    }
    for(j in 1:q){
      if((p[j]-3)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])+c/(b+norm_vec(X[k,]-mean(X[k,])*unlist(ones)))<=1)
        OUR2[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=(1-((p[j]-3)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])+c/(b+norm_vec(X[k,]-mean(X[k,])*unlist(ones)))))*(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])+mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]]
      else
        OUR2[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]]
    }
    prod=1
    for(j in 1:q){
      prod=prod*(q*norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])/norm_vec(X[k,]-mean(X[k,])*unlist(ones)))^(p[j]/2)
    }
    GammaX=(1+(1-gamma)/gamma*prod)^(-1)
    for(j in 1:q){
      OURNEW[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]=X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-
        (GammaX*(p[[j]]-3)/norm_vec(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])*(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])])*ones[[j]])+
           (1-GammaX)*(sum(p)-3)/norm_vec(X[k,]-mean(X[k,])*unlist(ones))*(X[k,(j*p[[j]]-p[[j]]+1):(j*p[[j]])]-mean(X[k,])*ones[[j]]))*sigma
      
    }
    L_JS_Comb1[k]=norm_vec(unlist(theta)-JS_Comb1[k,])
    L_JS_Comb2[k]=norm_vec(unlist(theta)-JS_Comb2[k,])
    L_JS_Sep1[k]=norm_vec(unlist(theta)-JS_Sep1[k,])
    L_JS_Sep2[k]=norm_vec(unlist(theta)-JS_Sep2[k,])
    L_OUR1[k]=norm_vec(unlist(theta)-OUR1[k,]) 
    L_OUR2[k]=norm_vec(unlist(theta)-OUR2[k,]) 
    L_OURNEW[k]=norm_vec(unlist(theta)-OURNEW[k,]) 
    L_MLE[k]=norm_vec(unlist(theta)-X[k,])}
  
  
  LL=cbind(L_JS_Comb1,L_JS_Comb2,L_JS_Sep1,L_JS_Sep2,L_OUR1,L_OUR2,L_MLE,L_OURNEW)
  R_mean[l,]=apply(LL,2,mean)
  R_sd[l,]=apply(LL,2,sd)
  theta_norm[l]=sqrt(norm_vec(unlist(theta)))  # add # when bayes risk
}

rank=rank(theta_norm) 
for(i in 1:LofRisk){
  R_meanSORT[rank[i],]=R_mean[i,]
  R_sdSORT[rank[i],]=R_sd[i,]}

print('break pint 2')

save(tau,rank,R_mean,R_sd, R_meanSORT,R_sdSORT,LofRisk,K,theta_norm,file="/home/wang.1964/BayesRiskComparison/SimMean.Rdata")