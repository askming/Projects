### given a set of posterior samples of model parameters loop over it, in each iteration run a 5000-iteration MCMC to get the posterior sample of random effects

## input data
# longitudinal outcome: Y; survival outcome: Ti, Deltai
# posterior sample of the parameters: post_theta 

# load JAGS function to run QRJM
source(run_jags)

row_col = dim(post_theta)
# post_u = null 

# load data for specific subject
y = Y 
Ti = Ti 
death = Deltai

for (i in seq_along(post_theta)){

	parameter = post_theta[i,] # ith posterior sample of the parameter

	# predict new y
	new_y = 
	# sample random effect
	post_re = runJAGS(y, Ti, death, parameter, iter=5000, burnin=2500)

	# calculate pi

	# calculate sensitivity

}
