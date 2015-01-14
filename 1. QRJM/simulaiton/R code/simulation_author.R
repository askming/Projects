library(mvtnorm)
n=250
Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
u=rmvnorm(n,c(0,0),Sigma)
X=matrix(1,n,2)
X[,2]=rnorm(n)
## Z=c(1,t) 
H=matrix(rnorm(2*n),n,2)
W=matrix(rnorm(2*n),n,2)
alpha1=1
alpha2=1
delta = gamma = c(1,1)

surv=function(t) {
res=((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
exp(-res)}

rnd=runif(n) 
Ti=rep(Inf,n) 
w=which(1-surv(101)-rnd>0)
Ti[w]=sapply(w,function(j) uniroot(function(x) 1-surv(x)[j]-rnd[j],lower=0,upper=101)$root)
Ti[Ti==0]=min(Ti[Ti>0])/2
Ci=5*rbeta(n,4,1)
Delta=rep(1,n)
Delta[Ti>Ci]=0
Ti[Ti>Ci]=Ci[Ti>Ci]


Ti_sim = function(n=250, alpha1, alpha2, delta=c(1,1), gamma=c(1,1)){
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u=rmvnorm(n,c(0,0),Sigma)
	H = matrix(rnorm(2*n),n,2)
	W = matrix(rnorm(2*n),n,2)

	# define the survival function
	surv = function(t) {
		if(alpha1==0 & alpha2==0) { res = t*exp(W%*%gamma) }

		res = ((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
		
		exp(-res)
	}

	rnd =runif(n) 
	Ti = rep(Inf,n) 
	w = which(1-surv(101)-rnd>0)
	print(length(w))
	Ti[w] = sapply(w,function(j) uniroot(function(x) 1-surv(x)[j]-rnd[j],lower=0,upper=101)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = 5*rbeta(n,4,1)
	# Delta = rep(1,n)
	# Delta[Ti>Ci] = 0
	Delta2 = Ti < Ci
	print(sum(Ti>Ci))
	# Ti[Ti>Ci] = Ci[Ti>Ci]
	Ti2 = pmin(Ti, Ci)

	list(Ti=Ti, Ti2=Ti2, event=Delta, delta=Delta2, H=H, W=W)
}

sur_data = Ti_sim(alpha1=1, alpha2=1)

library(survival)
coxfit = coxph(Surv(Ti, event) ~ W + H, data=sur_data)
summary(coxfit) 




