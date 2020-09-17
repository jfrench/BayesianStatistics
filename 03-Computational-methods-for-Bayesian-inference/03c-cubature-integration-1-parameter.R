# An example of finding the posterior mean
# From an unnormalized density

# Data distribution: y_1, ..., y_n | \theta ~ i.i.d Poisson(theta)
# Prior distribution: theta ~ Exp(1)

# Generate some synthetic data
# poisson with mean of 0.5
# in practice, we just have some data
set.seed(3) # for reproducibility
y = rpois(10, lambda = 0.5)

# unnormalized posterior
qtheta = function(theta, y) {
  prod(dpois(y, lambda = theta)) * dexp(theta, rate = 1)
}

# vectorize qtheta
vqtheta = Vectorize(qtheta, vectorize.args = "theta")
# plot of unnormalized posterior
mytheta = seq(0, 10, len = 1000)
plot(mytheta, vqtheta(mytheta, y = y), type = "l",
     xlab = expression(theta), ylab = expression(q(theta*"|"*y)))

# determine normalizing constant deterministically
# plot above suggests 10 is adequate
# Note: From integrate documentation
# f must accept a vector of inputs and
# produce a vector of function evaluations at
# those points. The Vectorize function may be
# helpful to convert f to this form.
(const = integrate(f = vqtheta, lower = 0, upper = 10, y = y))

# posterior density function
dpost = function(theta, y) {
  qtheta(theta, y)/const$value
}

# vectorized version of dpost for integrate function
vdpost = Vectorize(dpost, vectorize.args = "theta")

# double-check that posterior is proper!
integrate(vdpost, lower = 0, upper = 10, y = y)

# determine mean of posterior
# define function to integrate over to find the posterior mean
e = function(theta) {
  theta * vdpost(theta, y = y)
}

# Mean of posterior
integrate(e, 0, 10)

# the true mean is
a = 1 + sum(y)
b = 1 + length(y)
a/b
