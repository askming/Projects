
rm(list=ls())
set.seed(1234)
library(lme4)
library(pscl)
library(HI)

expit <- function(x) exp(x)/(1+exp(x))

I <- 1000
J <- 5
beta0 <- -0.5
beta1 <- 1
p <- 2 # p: number of regression parameters
sigma.u2 <- 1.69
ncycles <- 5000  # number of MCMC cycles
myshape <- 0.5*I-1 # the alpha parameter for the full conditional of parameter sigma2 

y <- matrix(NA, nrow=I, ncol=J)
U <- rnorm(I, sd=sqrt(sigma.u2))
covariate <- rnorm(I)
subject <- as.vector(mapply(rep, 1:I, J)); str(subject)

for (i in 1:I) y[i,] <- rbinom(J, size=1, prob=expit(beta0+beta1*covariate[i]+U[i]))
str(y) # 1000 X 5
table(y)
head(y)
y <- as.vector(t(y)) # change to long format
covariate <- as.vector(mapply(rep, covariate, J)) # change to long format

design.matrix <- cbind(rep(1,length(covariate)), covariate)
nobs <- rep(J, 1000) # number of observations for each subject

output <- data.frame(ID=subject, y = y, x=covariate)
write.table(output, file="mydata.txt", sep="\t", col.names=T, row.names=F)

gm2 <- glmer(y ~ covariate + (1 | subject), family = binomial)
summary(gm2)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -0.47521    0.05402  -8.797   <2e-16 ***
#covariate    1.09507    0.05949  18.407   <2e-16 ***
# Groups  Name        Variance Std.Dev.
# subject (Intercept) 1.6562   1.2869  

fit1 <- glmmPQL(y ~ covariate, random = ~1 | subject, family = binomial)
summary(fit1)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -0.47521    0.05402  -8.797   <2e-16 ***
#covariate    1.09507    0.05949  18.407   <2e-16 ***
# Groups  Name        Variance Std.Dev.
# subject (Intercept) 1.6562   1.2869  

