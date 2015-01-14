rm(list=ls())
########################################################################
---------------------------- 1. simulate data --------------------------
########################################################################

###############################################################
############ function to simulate survival time ###############
###############################################################
# survival function is given by: S(t)= exp(- exp(B) * (exp(A*t) - 1) ) / A), where
# B = alpha[1] * delta[1] * H[,1] + alpha[2] * U[,1] + gamma %*% t(W)
# A = alpha[2] * U[,2] + alpha[1] * delta[2] * H[,2]

library(mvtnorm)
Ti_sim = function(n=250, alpha1, alpha2, delta=c(1,1), gamma=c(1,1)){
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u=rmvnorm(n,c(0,0),Sigma)
	H = matrix(rnorm(2*n),n,2)
	W = matrix(rnorm(2*n),n,2)

	# define the survival function
	surv = function(t) {
		if(alpha1==0 & alpha2==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha2*u[,1]+alpha2*t*u[,2]+alpha1*delta[1]*H[,1]+alpha1*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha1*H[,1]*delta[1]+W%*%gamma+alpha2*u[,1]))/(alpha2*u[,2]+alpha1*H[,2]*delta[2]))
		}
		
		return(exp(-res))
	}

	rnd =runif(n) 
	Ti = rep(Inf,n) 
	# w = which(surv(1e8)-rnd < 0)
	w = which(1-surv(101)-rnd > 0)
	Ti[w] = sapply(w,function(j) uniroot(function(x) 1-surv(x)[j]-rnd[j],lower=0,upper=101)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = 5*rbeta(n,4,1)
	Delta = as.numeric(Ti < Ci)
	Ti2 = pmin(Ti, Ci)

	list(Ti=Ti2, event=Delta, H=H, W=W)
}

# Ti_sim2 = function(n=250, gamma=c(1,1), k=1){
# 	W = matrix(rnorm(2*n), nrow=n)
# 	time<-numeric(0)
# 	event<-numeric(0)
#   	S<-runif(n)
#   	T<- as.vector(-log(S)/(k*exp(W%*%gamma)))
#   	C <- runif(n,2,5)
#   	time<-pmin(T, C)
#   	event<-as.numeric(T < C) 

# 	list(Ti=time, event=event, W=W)
# }

sur_data = Ti_sim(alpha1=1, alpha2=0)
# sur_data2 = Ti_sim2()

#### try using coxph(){survival}
# library(survival)
# coxfit = coxph(Surv(Ti, event) ~ W + H, data=sur_data)
# summary(coxfit) 

########################################################################
---------------------------- 2. run the model --------------------------
########################################################################
library(R2jags)

QRJM_jags = function(data, I=250){
	setwd("/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files")
	model = "/Users/askming/Documents/github/Projects/1.\ QRJM/simulaiton/model\ files/QRJM_survival01.txt"
	
	# load survival data
	Ti = data[['Ti']]
	H = data[['H']]
	W = data[['W']]
	event =  data[['event']]
	zeros = rep(0,I)

	jags.data = list(W=W, Ti=Ti, event=event, I=I, zeros=zeros, H=H)
 	jags.params = c("delta","gamma","alpha1","alpha2","c")
 	jags.inits = function(){ list(gamma=c(1,1),c=1,delta=c(1,1), alpha1=1, alpha2=1) }	
  	jags(data=jags.data, inits=jags.inits, jags.params, n.iter=10000, n.burnin=5000, model.file=model)	
}

surfit = QRJM_jags(data=sur_data)
print(surfit)




