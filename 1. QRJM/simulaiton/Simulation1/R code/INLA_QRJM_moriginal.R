### Joint model using INLA ###
### 2015-01-09 ###


## longitudianl outcome can have missing values

# simulate survival data
sim_Ti = function(n=250, alpha1, alpha2, delta=c(1,1), gamma=c(1,1)){
	library(mvtnorm)
	# random effects
	Sigma=matrix(c(0.09, 0.09*0.16, 0.09*0.16, 0.09), nrow=2, byrow=T)
	u = rmvnorm(n,c(0,0),Sigma)
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

	list(Ti=Ti2, event=Delta, H=H, W=W, U=u)
}

# simulate quantile longitudinal data in long format
sim_longitudinal_data = function(survival_data=surdata, n=250, time=c(0, 0.25, 0.5, 0.75, 1, 3), tau, sigma=1, beta=c(1,1), delta=c(1,1)){
	# survival_data - data simulated from survival model
	# n - # of subjects
	# time - time points of observations
	# tau - quantile
	# sigma - scale parameter
	n_obs = length(time)
	time_long = rep(time, n) # at most # = length(time) observations per patient
	H = U = matrix(NA, n*n_obs, 2)
	Y = rep(NA, n*n_obs) # long format
	U[,1] = rep(survival_data$U[,1], each=n_obs) # random effects
	U[,2] = rep(survival_data$U[,2], each=n_obs)
	H[,1] = rep(survival_data$H[,1], each=n_obs)
	H[,2] = rep(survival_data$H[,2], each=n_obs)*time_long
	X = cbind(1, rep(rnorm(n), each=n_obs))
	Ti = survival_data$Ti
	count = sapply(Ti, function(x) sum(x > time)) # number of observations after drop-outs
	end = seq(6, n*n_obs, 6)

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
		Y[i] = X[i,] %*% beta + H[i,] %*% delta + U[i,] %*% c(1, time_long[i]) + rald(1, scale=sigma, p=tau)	

	}

	# create drop out missing
	for(i in seq_along(end)){
		if(count[i]!=6){
					Y[min(count[i]+6*(i-1)+1,end[i]):end[i]] = NA
		}

	}

	obs_id = rep(seq(1:6), n)
	data.frame(obs_id, Y, X, H=H, time_long)
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
		sur_data = sur_fun(n=n, alpha1=alpha[1], alpha2=alpha[2])
		longi_data = longi_fun(n=n, sur_data, tau=tau)
		outdata[[i]] = list(survival_data=sur_data, longitudinal_data=longi_data)		
	}
	outdata	
}

# surdata = sim_Ti(alpha1=1, alpha2=1)
# long_data = sim_longitudinal_data(surdata, tau=0.25)
multi_data = sim_multiple_data(10, n=500, alpha=c(0,1), tau=0.25)
# head(long_data, 24)


run_multi_inla = function(data){
	library(INLA)
	long_data = data$longitudinal_data
	surdata = data$survival_data

	## prepare data set
	N = 250
	ng =dim(long_data)[1]
	ns = N

	## prepare the response variable
	y_long = c(long_data$Y, rep(NA, ns))
	y_surv = inla.surv(time=c(rep(NA, ng), surdata$Ti), event=c(rep(NA, ng), surdata$event))
	Yjoint = list(y_long, y_surv)

	# prepare the fixed covariate
	linear.covariate = data.frame(mu = as.factor(c(rep(1, ng), rep(2, ns))),
                               b1.X1 = c(long_data$X1, rep(0, ns)),
                               b2.X2 = c(long_data$X2, rep(0, ns)),
                               d1.H1 = c(long_data$H.1, rep(0, ns)),
                               d2.H2 =  c(long_data$H.2, rep(0, ns)), 
                               time = c(long_data$time_long, rep(0, ns)),
                               g1.W1 = c(rep(0, ng), surdata$W[,1]),
                               g2.W2 = c(rep(0, ng), surdata$W[,2]),
                               Ti = c(rep(0, ng), surdata$Ti),
                               d12.H1 = c(rep(0, ng), surdata$H[,1]),
                               d22.H2 = c(rep(0, ng), surdata$H[,2]*surdata$Ti)
                               )


	# prepare random effects
	random.covariate <- list(U11 = c(rep(1:N, each=6),rep(NA, ns)),
                         U21 = c(rep(N+(1:N), each=6),rep(NA, ns)),
                         U12 = c(rep(NA,ng), 1:N),
                         U22 = c(rep(NA,ng), N+(1:N))
                         )


	idh1 = idh2 = c(rep(1, ng), rep(NA, ns))
	idh12 = idh22 = c(rep(NA, ng), rep(1, ns))

	joint.data <- c(linear.covariate,random.covariate)
	joint.data$Y <- Yjoint
	joint.data$idh1 = idh1
	joint.data$idh2 = idh2
	joint.data$idh12 = idh12
	joint.data$idh22 = idh22

	
	### MODEL
	formula = Y ~ b1.X1 + b2.X2 + g1.W1 + g2.W2 - 1 +
	f(idh1, d1.H1, hyper = list(prec = list(initial = -4, fixed=TRUE))) +
	f(idh2, d2.H2, hyper = list(prec = list(initial = -4, fixed=TRUE))) +
	f(idh12, d12.H1, copy="idh1", hyper = list(beta = list(fixed=FALSE))) +
	f(idh22, d22.H2, copy="idh2", hyper = list(beta = list(fixed=FALSE))) +
	f(U11, model="iid2d", n=2*N) +
	f(U21, time, copy="U11") +
	f(U12, copy="U11", fixed=FALSE) +
	f(U22, Ti, copy="U11", same.as="U12") 


	# enable.model.likelihood.laplace=TRUE
	mod = inla(formula, family=c("laplace", "exponential"), data = joint.data, verbose=TRUE, control.family=list(list(alpha=0.25), list()))

	out = list(summary(mod)$fixed, summary(mod)$hyperpar)

}

t1=Sys.time()
result = lapply(multi_data, run_multi_inla)
t2=Sys.time()
t2-t1

result = lapply(result, function(x)list(summary(x)$fixed_effects, summary(x)$hyperpar))


# another example
# data(Munich)
# g = system.file('demodata/munich.graph', package='INLA')

# formula = rent ~ f(location, model='besag', graph.file=g, param=c(1, .001)) +
# f(year, model='crw2', values=seq(1918,2001), param=c(1, .001)) + 
# f(floor.size, model='crw2', param=c(1,.001)) + Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh + Kein.Badkach + Besond.Bad + Gehobene.Kueche + zim1 + zim2 + zim3 + zim4 + zim5 + zim6 - 1

# mod = inla(formula, data = Munich, verbose = TRUE, family = 'laplace',
# control.family=list(alpha=.5, gamma=2, epsilon=.01), control.predictor = list(initial = 12), control.inla = list(h=1e-4))

sum = function(result){

		# fix_effects_mean = sapply(result, function(x) x[[1]][,1])
		# fix_effects_sd = sapply(result, function(x) x[[1]][,2])
		re = sapply(result, function(x) x[[2]][5,1:2])
		# mean_fix_mean = round(rowMeans(fix_effects_mean), 3)
		# mean_fix_sd = round(rowMeans(fix_effects_sd), 3)
		re_mean = round(rowMeans(re), 3)

		list(alpha2 = re_mean)


}

sum(result)


