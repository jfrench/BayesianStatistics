# An approximating the posterior with a normal
# From an unnormalized density

# Data distribution: y_1, ..., y_n | theta ~ i.i.d Poisson(theta)
# Prior distribution: theta ~ Exp(1)

# Generate some synthetic data
# poisson with mean of 0.5
# in practice, we just have some data
set.seed(3) # for reproducibility
y = rpois(10, lambda = 0.5)

# unnormalized posterior
qtheta = function(theta, y) {
  prod(dpois(y, lambda = theta)) *
    dexp(theta, rate = 1)
}

# vectorize qtheta
vqtheta = Vectorize(qtheta, vectorize.args = "theta")

# determine normalizing constant deterministically
# plot above suggests 10 is adequate
# Note: From integrate documentation
# f must accept a vector of inputs and
# produce a vector of function evaluations at
# those points. The Vectorize function may be
# helpful to convert f to this form.
(nconst = integrate(f = vqtheta, lower = 0, upper = 10, y = y))

# posterior density function
dpost = function(theta, y) {
  qtheta(theta, y)/nconst$value
}

# vectorized version of dpost for integrate function
vdpost = Vectorize(dpost, vectorize.args = "theta")

# double-check that posterior is proper!
integrate(vdpost, lower = 0, upper = 10, y = y)

# determine map
(map = optimize(qtheta, interval = c(0, 10), y = y, maximum = TRUE))
map$maximum

# normal approximation posterior
dpapprox = function(theta, y, thetahat) {
  dnorm(theta, mean = thetahat, sd = sqrt(thetahat^2/sum(y)))
}

# vectorize over theta
vdapprox = Vectorize(dpapprox, vectorize.args = "theta")

# range of theta values
theta = seq(0, 2, length = 1000)
# plot true density
plot(theta, vdpost(theta, y), ylab = "density", type = "l", col = "orange")
# plot normal approximation
lines(theta, vdapprox(theta, y, thetahat = map$maximum), col = "blue")
legend("topright", legend = c("true", "approximation"),
       col = c("orange", "blue"), lty = 1)


# Generate some synthetic data
# poisson with mean of 0.5
# in practice, we just have some data
set.seed(3) # for reproducibility
y = rpois(100, lambda = 0.5)

# unnormalized posterior
qtheta = function(theta, y) {
  log_pdata = sum(dpois(y, lambda = theta, log = TRUE))
  log_pprior = dexp(theta, rate = 1, log = TRUE)
  exp(log_pdata + log_pprior)
}

# vectorize qtheta
vqtheta = Vectorize(qtheta, vectorize.args = "theta")

# determine normalizing constant
(nconst = integrate(f = vqtheta, lower = 0, upper = 10, y = y))

# there is an issue numerically
# the log_pdata is too small (exp(log_pdata) approx 0)
# we'll add a constant to deal with this
# how to choose?
# trial and error
# too big and everything blows up
# too small and numerical underflow remains
# unnormalized posterior
# add constant to deal with numerical underflow
qtheta_robust = function(theta, y) {
  log_pdata = sum(dpois(y, lambda = theta, log = TRUE))
  log_pprior = dexp(theta, rate = 1, log = TRUE)
  exp(log_pdata + log_pprior + 80)
}

# vectorize qtheta
vqtheta_robust = Vectorize(qtheta_robust, vectorize.args = "theta")

# determine normalizing constant
(nconst = integrate(f = vqtheta_robust, lower = 0, upper = 10, y = y))

# posterior density function
dpost = function(theta, y) {
  qtheta_robust(theta, y)/nconst$value
}

# vectorized version of dpost for integrate function
vdpost = Vectorize(dpost, vectorize.args = "theta")

# double-check that posterior is proper!
integrate(vdpost, lower = 0, upper = 10, y = y)

# determine map
(map = optimize(qtheta, interval = c(0, 10), y = y, maximum = TRUE))
map$maximum

# normal approximation posterior
dpapprox = function(theta, y, thetahat) {
  dnorm(theta, mean = thetahat, sd = sqrt(thetahat^2/sum(y)))
}

# vectorize over theta
vdapprox = Vectorize(dpapprox, vectorize.args = "theta")

# range of theta values
theta = seq(0, 2, length = 1000)
# plot true density
plot(theta, vdpost(theta, y), ylab = "density", type = "l", col = "orange")
# plot normal approximation
lines(theta, vdapprox(theta, y, thetahat = map$maximum), col = "blue")
legend("topright", legend = c("true", "approximation"),
       col = c("orange", "blue"), lty = 1)

