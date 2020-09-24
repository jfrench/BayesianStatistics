# library(mvtnorm)
# library(cubature)

# An example of finding the posterior mean
# for two-parameter posterior

# Data distribution: y_1, ..., y_n | mu, sigma^2
# ~ i.i.d. N(mu, sigma^2)
# Prior distribution: mu ~ U(10, 15)
# sigma^2 ~ N(0.5, 0.25) * I(sigma^2 > 0)
# Essentially a truncated normal with a
# different scaling constant

# Generate some synthetic data
# in practice, we just have some data
set.seed(7)
y = rnorm(100, mean = 11, sd = 0.47)

# create function for unnormalized posterior
# takes the vector theta = (mu, sigma^2) and y = data.
# add constant to handle numerical underflow
qtheta = function(theta, y) {
  mu = theta[1]
  sigma = sqrt(theta[2])
  log_pdata = sum(dnorm(y, mean = mu, sd = sigma, log = TRUE))
  log_pmu = dunif(mu, 10, 15, log = TRUE)
  log_psigmasq = dnorm(sigma^2, mean = 0.5, sd = 0.5, log = TRUE)
  log_qtheta = log_pdata + log_pmu + log_psigmasq
  exp(log_qtheta + 70)
}

# determine normalizing constant of posterior
(py = cubintegrate(f = qtheta,
                   lower = c(10, 0),
                   upper = c(15, 2),
                   y = y))
const = py$integral

# posterior distribution
dpost = function(theta, y, const) {
  qtheta(theta, y)/const
}

# determine map estimates
# log of unnormalized density
log_qtheta = function(theta, y) {
  mu = theta[1]
  sigma = sqrt(theta[2])
  log_pdata = sum(dnorm(y, mean = mu, sd = sigma, log = TRUE))
  log_pmu = dunif(mu, 10, 15, log = TRUE)
  log_psigmasq = dnorm(sigma^2, mean = 0.5, sd = 0.5, log = TRUE)
  log_pdata + log_pmu + log_psigmasq
}

# fnscale = -1 is to make optim perform maximization instead of minimization
map = optim(par = c(12.5, 1),
            f = log_qtheta,
            lower = c(10.0001, 0.05),
            upper = c(14.9999, 2),
            method = "L-BFGS-B",
            y = y,
            control = list(fnscale = -1))

# return the observed information matrix
obs_I = function(thetahat, n) {
  mu = thetahat[1]
  sigma = sqrt(thetahat[2])
  # 2nd derivative of log likelihood w/r to mu
  d2dmu2 = -n/sigma^2
  # derivative of log likelihood w/r to mu and sigmasq
  d2dmudsigmasq = -sum((y - mu))/sigma^4
  # 2nd derivative of log likelihood w/r to sigmasq
  d2dsigmasq2 = n/(2 * sigma^4) - 1/sigma^6 * sum((y-mu)^2)
  #
  H = cbind(c(d2dmu2, d2dmudsigmasq), c(d2dmudsigmasq, d2dsigmasq2))
  solve(-H)
}

# posterior approximation
dpapprox = function(theta, thetahat, Ihat) {
  mvtnorm::dmvnorm(x = theta, mean = thetahat, sigma = Ihat)
}

# compare results
# create sequence of values for mu and sigmasq
mymu = seq(10.95, 11.2, length = 200)
mysigmasq = seq(0.14, 0.29, length = 200)
# create grid
mytheta = expand.grid(mymu, mysigmasq)
# for each row of mytheta, plug it into dpost function
z = apply(mytheta, 1, function(theta) {
  dpost(theta, y = y, const = const)
})
# convert z to matrix for plotting
zmat = matrix(z, nrow = length(mymu))

# compute observed Information outside of loop
Ihat = obs_I(map$par, length(y))
zhat = mvtnorm::dmvnorm(x = mytheta, mean = map$par, sigma = Ihat)
# convert z to matrix for plotting
zhatmat = matrix(zhat, nrow = length(mymu))

# side-by-side results
par(mfrow = c(1, 2))

# create heat map of objective surface
image(mymu, mysigmasq, zmat,
      col = hcl.colors(64, "YlOrRd", rev = TRUE),
      xlab = expression(mu), ylab = expression(sigma^2),
      zlim = range(c(z, zhat)))
# add contours
contour(mymu, mysigmasq, zmat,
        add = TRUE)
title("true posterior")

# create heat map of approximation surface
image(mymu, mysigmasq, zhatmat,
      col = hcl.colors(64, "YlOrRd", rev = TRUE),
      xlab = expression(mu), ylab = expression(sigma^2),
      zlim = c(range(c(z, zhat))))
# add contours
contour(mymu, mysigmasq, zhatmat, add = TRUE)
title("normal approximation")
par(mfrow = c(1, 1))