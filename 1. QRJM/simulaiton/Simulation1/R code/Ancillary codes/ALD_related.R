###############################################################
################# functions related with ALD  #################
###############################################################
# desity function of ALD(mu, sigma,p)
dald = function(x, p, mu, sigma){
	# density function of ALD(mu, sigma, tau) based on the definition
	rho = (abs(x-mu) + (2*p -1)*(x-mu)) / (2*sigma)
	den = p*(1-p)/sigma*exp(-rho)
}

r_ald = function(n, location=0, sigma, p){
	# random variable simulation based on normal+exponetial mixture
		u = rnorm(n)
		z = rexp(n)
		v = sigma*z
		theta = (1-2*p)/(p*(1-p))
		tau = sqrt(2/(p*(1-p)))
		sample = theta * z + tau * sqrt(z) * u
	}

# check the validity of above functions
x = seq(-20,20,0.01)
plot(x, dald(x, sigma=1, mu=0, p=0.25), ylim=c(0,0.45), main='', ylab="density",type="l", lwd=2)
lines(x, dald(x, p=0.75, mu=0, sigma=1), col="blue", lty=2, lwd=2)
lines(x, dald(x, p=0.5, mu=0, sigma=1), col="red",lty=3, lwd=2)
lines(x, dnorm(x, 0, 1),col="green",lty=4, lwd=2)
legend("topright", legend=c("N(0,1)", "LD(0, 1)", "ALD(0, 1, 0.75)", "ALD(0, 1, 0.25)"), lty=c(4, 3, 2, 1), col=c("green", "red", "blue", "black"), lwd=c(2,2,2,2), cex=0.75)