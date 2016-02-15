### Joint model using INLA ###
### 2015-01-11 ###
### Ming Yang ###


# simulate survival data
sim_Ti = function(n=250, alpha=c(1,1,1), delta=c(1,1), gamma=c(1,1)){
	library(mvtnorm)
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u = rmvnorm(n,c(0,0),Sigma)	# random effects
	H = matrix(rnorm(2*n),n,2)	# fixed effects
	W = matrix(rnorm(2*n),n,2)	# fixed effects

	# define the survival function
	surv = function(t) {
		if(alpha[1]==0 & alpha[2]==0 & alpha[3]==0) {res=t*exp(W%*%gamma)}
		else{
			res = ((exp(alpha[2]*u[,1]+alpha[3]*u[,2]+alpha[1]*delta[1]*H[,1]+alpha[1]*delta[2]*H[,2]*t+W%*%gamma)-exp(alpha[1]*H[,1]*delta[1]+W%*%gamma+alpha[2]*u[,1]+alpha[3]*u[,2]))/(alpha[1]*H[,2]*delta[2]))
		}
		
		return(exp(-res))
	}

	rnd =runif(n) 
	Ti = rep(Inf,n) 
	w = which(1-surv(101)-rnd > 0)
	Ti[w] = sapply(w,function(j) uniroot(function(x) 1-surv(x)[j]-rnd[j],lower=0,upper=101)$root)
	Ti[Ti==0] = min(Ti[Ti>0])/2
	Ci = 5*rbeta(n,4,1) # censoring time, from beta dist.
	event = as.numeric(Ti < Ci) # event indicator
	Ti2 = pmin(Ti, Ci) # final survival time after censoring

	data.frame(Ti=Ti2, event=event, H=H, W=W, U=u) #output
}

# simulate quantile longitudinal data in long format
## longitudianl outcome have monotone missing pattern
sim_longitudinal_data = function(survival_data=surdata, n=250, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1), delta=c(1,1)){

	n_obs = length(time) # number of repeated measures per subject
	time_long = rep(time, n) 
	H = U = matrix(NA, n*n_obs, 2)
	Y = Y2= rep(NA, n*n_obs) # long format
	U[,1] = rep(survival_data$U.1, each=n_obs) # shared random effects, from survival data
	U[,2] = rep(survival_data$U.2, each=n_obs)
	H[,1] = rep(survival_data$H.1, each=n_obs) # shared fixed effect from survival data
	H[,2] = rep(survival_data$H.2, each=n_obs)*time_long
	X = cbind(1, rep(rnorm(n), each=n_obs)) # fixed effect only in longitudianl model
	
	# random variable generator: ALD(0,sigma,p)
	rald = function(n, location=0, scale, p){
		u = rnorm(n)
		z = rexp(n)
		v = scale*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		sample = theta * z + tau * sqrt(z) * u
	}

	# create outcome
	for (i in 1:length(Y)) {
		Y[i] = X[i,] %*% beta + U[i,] %*% c(1, time_long[i]) + rald(1,scale=sigma, p=tau)	

	}

	# Y2 use gaussian error
	for (i in 1:length(Y)) {
		Y2[i] = X[i,] %*% beta + U[i,] %*% c(1, time_long[i]) + rnorm(1)	

	}

	# create monotoe drop out missing
	Ti = survival_data$Ti # surival time 
	count = sapply(Ti, function(x) sum(x > time)) # number of observations for each subject after drop-outs
	end = seq(6, n*n_obs, 6)
	for(i in seq_along(end)){
		if(count[i]!=6){
					Y[min(count[i]+6*(i-1)+1,end[i]):end[i]] = NA
		}

	}
	obs_id = rep(seq(1:6), n)
	data.frame(obs_id, Y, Y2, X, time_long)
}


###############################################################
#----- function to simulate multiple joint data sets ---------#
###############################################################
sim_multiple_data = function(N, n, sur_fun=sim_Ti, longi_fun=sim_longitudinal_data, alpha, tau){
	# N - number of data sets to generate
	# sur_fun - function to simulate survival data
	# longi_fun - function to sumualte longitudinal data
	# alpha - association mechanism for JM
	# tau - quantile

	outdata = vector(mode='list', N)
	for (i in 1:N){
		sur_data = sur_fun(n=n, alpha)
		longi_data = longi_fun(n=n, sur_data, tau=tau)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)		
	}
	outdata	
}

# test: simulate single data
# surdata = sim_Ti()
# long_data = sim_longitudinal_data(surdata, tau=0.25)

library(INLA)
enable.model.likelihood.laplace = TRUE

# function to run simulation in INLA
run_inla = function(data){
	long_data = data$longitudinal_data
	surdata = data$survival_data

	## prepare data set
	N = 250
	ng =dim(long_data)[1]
	ns = N

	## prepare the response variable
	y_long = c(long_data$Y, rep(NA, ns))
	y_long2 = c(long_data$Y2, rep(NA, ns))
	y_surv = inla.surv(time=c(rep(NA, ng), surdata$Ti), event=c(rep(NA, ng), surdata$event))
	Yjoint1 = list(y_long, y_surv)
	Yjoint2 = list(y_long2, y_surv)

	# prepare the fixed covariate
	linear.covariate = data.frame(mu = as.factor(c(rep(1, ng), rep(2, ns))),
                               b1.X1 = c(long_data$X1, rep(0, ns)),
                               b2.X2 = c(long_data$X2, rep(0, ns)),
                               d1.H1 = c(long_data$H.1, rep(0, ns)),
                               d2.H2 =  c(long_data$H.2, rep(0, ns)), 
                               time = c(long_data$time_long, rep(0, ns)),
                               g1.W1 = c(rep(0, ng), surdata$W.1),
                               g2.W2 = c(rep(0, ng), surdata$W.2),
                               Ti = c(rep(0, ng), surdata$Ti),
                               d12.H1 = c(rep(0, ng), surdata$H.1),
                               d22.H2 = c(rep(0, ng), surdata$H.2*surdata$Ti)
                               )


	# prepare random effects
	random.covariate <- list(U11 = c(rep(1:N, each=6),rep(NA, ns)),
                         U21 = c(rep(N+(1:N), each=6),rep(NA, ns)),
                         U12 = c(rep(NA,ng), 1:N),
                         U22 = c(rep(NA,ng), N+(1:N))
                         )


	idh1 = idh2 = c(rep(1, ng), rep(NA, ns)) # longitudinal part
	idh12 = idh22 = c(rep(NA, ng), rep(1, ns)) # survival part

	joint.data <- c(linear.covariate,random.covariate)
	joint.data$Y <- Yjoint
	joint.data$Y2 <- Yjoint2
	joint.data$idh1 = idh1
	joint.data$idh2 = idh2
	joint.data$idh12 = idh12
	joint.data$idh22 = idh22

	
	### MODEL with ALD error
	formula = Y ~ b1.X1 + b2.X2 + g1.W1 + g2.W2 - 1 +
	f(idh1, d1.H1, hyper = list(prec = list(initial = -4, fixed=TRUE))) +
	f(idh2, d2.H2, hyper = list(prec = list(initial = -4, fixed=TRUE))) +
	f(idh12, d12.H1, copy="idh1", hyper = list(beta = list(fixed=FALSE))) +
	f(idh22, d22.H2, copy="idh2", hyper = list(beta = list(fixed=FALSE))) +
	f(U11, model="iid2d", n=2*N) + # random intercept
	f(U21, time, copy="U11") + # random slope in longitudinal model
	f(U12, copy="U11", fixed=FALSE) +
	f(U22, Ti, copy="U11", fixed=FALSE)  # random slope in surival model

	### MODEL with NORMAL error
	formula2 = Y2 ~ b1.X1 + b2.X2 + g1.W1 + g2.W2 - 1 +
	f(idh1, d1.H1, hyper = list(prec = list(initial = -4, fixed=TRUE))) +
	f(idh2, d2.H2, hyper = list(prec = list(initial = -4, fixed=TRUE))) +
	f(idh12, d12.H1, copy="idh1", hyper = list(beta = list(fixed=FALSE))) +
	f(idh22, d22.H2, copy="idh2", hyper = list(beta = list(fixed=FALSE))) +
	f(U11, model="iid2d", n=2*N) + # random intercept
	f(U21, time, copy="U11") + # random slope in longitudinal model
	f(U12, copy="U11", fixed=FALSE) +
	f(U22, Ti, copy="U11", fixed=FALSE)  # random slope in surival model


	# enable.model.likelihood.laplace=TRUE
	# mod1 = inla(formula, family=c("laplace", "exponential"), data = joint.data, verbose=TRUE, control.family=list(list(alpha=0.25), list()))

	mod2 = inla(formula2, family=c("gaussian", "exponential"), data = joint.data, verbose=TRUE)

	out = list(mod$summary.fixed, mod$summary.random, mod$summary.hyperpar)

}


# simulation experiment in multiple datasets
multi_data = sim_multiple_data(N=10, n=250, alpha=c(1,1,1), tau=0.25)

t1=Sys.time()
result = lapply(multi_data, run_inla)
Sys.time() - t1



# function to summarized the simulation result
sum = function(result){
		fe = sapply(result, function(x) x[[1]][,1])
		fe_sd = sapply(result, function(x) x[[1]][,2])
		idh1 = sapply(result, function(x) x[[2]]$idh1[,2])
		idh1_sd = sapply(result, function(x) x[[2]]$idh1[,3])
		idh2 = sapply(result, function(x) x[[2]]$idh2[,2])
		idh2_sd = sapply(result, function(x) x[[2]]$idh2[,3])
		hyper = sapply(result, function(x) x[[3]][5:8,1])
		hyper_sd = sapply(result, function(x) x[[3]][5:8,2])

		fe_mean = round(rowMeans(fe), 3)
		fe_sd_mean = round(rowMeans(fe_sd), 3)
		idh1_mean = round(mean(idh1), 3)
		idh1_sd_mean = round(mean(idh1_sd), 3)
		idh2_mean = round(mean(idh2), 3)
		idh2_sd_mean = round(mean(idh2_sd), 3)
		hyper_mean = round(rowMeans(hyper), 3)
		hyper_sd_mean = round(rowMeans(hyper_sd), 3)

		# list(fe_mean, fe_sd_mean, delta1=idh1_mean, delta2=idh2_mean, alpha=hyper_mean, alpha_sd=hyper_sd_mean)
		list(fe_mean, fe_sd_mean, delta1=list(mean=idh1_mean, sd=idh1_sd_mean),delta2=list(mean=idh2_mean, sd=idh2_sd_mean), hyper_mean, hyper_sd_mean)
}

sum(result)


