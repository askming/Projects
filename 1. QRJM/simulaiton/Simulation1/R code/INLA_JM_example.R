rm(list=ls())
n = 1000
x1 = rnorm(n)
x2 = rnorm(n)
eta1 = 1 + x1 + 2*x2
eta2 = 2 + 3*(x1 + 2*x2)
y1 = eta1 + rnorm(n)
y2 = eta2 + rnorm(n)

## code this as
## eta1 = ... + beta*x
## eta2 = ... + beta*x*ff
## where 'ff' is a scaling factor

Y = matrix(NA, 2*n, 2)
Y[1:n, 1] = y1
Y[n + 1:n, 2] = y2
X1 = c(x1, x1)
X2 = c(x2, x2)
intercept = as.factor(c(rep(1, n), rep(2, n)))
idx1 = c(rep(1, n), rep(NA, n))
idx2 = c(rep(1, n), rep(NA, n))
idxx1 = c(rep(NA, n), rep(1, n))
idxx2 = c(rep(NA, n), rep(1, n))

formula = Y ~ -1 + intercept +
    f(idx1, X1,
      hyper = list(prec = list(initial = -4, fixed=TRUE))) +
    f(idx2, X2,
      hyper = list(prec = list(initial = -4, fixed=TRUE))) +
    f(idxx1, X1, copy="idx1",
      hyper = list(beta = list(fixed=FALSE)))+
    f(idxx2, X2, copy="idx2",
      hyper = list(beta = list(fixed=FALSE)))

r = inla(formula, data = list(Y=Y, intercept = intercept, idx1=idx1, idx2=idx2, X1=X1, X2=X2),
    family = rep("gaussian", 2))

summary(r)
r$summary.fixed



n = 1000
x = rnorm(n)
eta1 = 1 + 1.5*x
eta2 = 2 + 2*(1.5*x)
y1 = eta1 + rnorm(n)
y2 = eta2 + rnorm(n)

## code this as
## eta1 = ... + beta*x
## eta2 = ... + beta*x*ff
## where 'ff' is a scaling factor

Y = matrix(NA, 2*n, 2)
Y[1:n, 1] = y1
Y[n + 1:n, 2] = y2
X = c(x, x)
intercept = as.factor(c(rep(1, n), rep(2, n)))
idx = c(rep(1, n), rep(NA, n))
idxx = c(rep(NA, n), rep(1, n))
formula = Y ~ -1 + intercept + 
    f(idx, X,
      hyper = list(prec = list(initial = -4, fixed=TRUE))) +
    f(idxx, X, copy="idx",
      hyper = list(beta = list(fixed=FALSE)))

r = inla(formula, data = list(Y=Y, intercept = intercept, idx=idx, X=X),
    family = rep("gaussian", 2))

