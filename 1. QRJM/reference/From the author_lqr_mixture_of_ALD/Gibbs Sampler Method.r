library(R2WinBUGS)

#create model
Gibsmodel<-function(){
  k1<-(1-2*tau)/tau*(1-tau)
  k2<-2/tau*(1-tau)
    for (i in 1:n){
      for (j in 1:t){
        my[i,j]~dnorm(mu[i,j], var[i,j])
        mu[i,j]<-inprod(mx[i,j,],beta[])+alpha[i]+k1*e[i,j]
        var[i,j]<-sigma/(k2*e[i,j]) # why is this: dexp(sigma) == exp(sigma)
        e[i,j]~dexp(sigma)
      }
    alpha[i]~dnorm(0,alpha.phi)
    }
  beta[1]~dnorm(0, 0.01)
  beta[2]~dnorm(0, 0.01)
  beta[3]~dnorm(0, 0.01)
  beta[4]~dnorm(0, 0.01)
  sigma~dgamma(0.1,0.1)
  alpha.phi~dgamma(0.1,0.1)
  phi.squared<-1/alpha.phi
  phi<-sqrt(phi.squared)
}  

if (is.R()){ # for R
    ## some temporary filename:
    filename <- file.path(tempdir(), "Gibsmodel.bug")}
## write model file:
write.model(Gibsmodel, filename)
## and let's take a look:
##file.show(filename)

## Generate data
n=5
t=30
Tbeta=c(5,6,7,8)
mx=array(0,dim=c(n,t,4))
for (p in 1:4)
  {mx[,,p]=rnorm(n*t,0,1)}
e=rnorm(n*t,0,1)
gamma=rnorm(n,0,2)
my=array(0,dim=c(n,t))
me=matrix(e,nrow=n,ncol=t,byrow=TRUE)
for (i in 1:n)
{
  for (j in 1:t)
   {
    my[i,j]=mx[i,j,]%*%Tbeta+gamma[i]+me[i,j]
    }
 }

#Run model
 tau=0.5
 datax=list("n","t","my","mx","tau")
 parameters<-c("beta[1:4]","alpha","sigma","phi")
 inits=function(){list(e<-matrix(rep(0.1,times=n*t),nrow=n,ncol=t),beta<-c(1,1,1,1),alpha<-rep(0.1,times=n),sigma<-0.1,alpha.phi<-0.1)}
 gib.sim<-bugs(datax,inits,parameters,model.file="Gibsmodel.bug",n.chains=2,n.iter=10000,n.burnin=5000)
 attach.bugs(gib.sim)
 gib.sim


