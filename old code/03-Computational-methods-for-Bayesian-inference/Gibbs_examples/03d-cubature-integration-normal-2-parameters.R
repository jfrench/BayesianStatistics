# An example of finding the posterior mean
# for two-parameter posterior

# need to do multivariate integration (easily)
library(cubature)

# Data distribution: y_1, ..., y_n | mu, sigma^2
# ~ i.i.d. N(mu, sigma^2)
# Prior distribution: mu ~ U(10, 15)
# sigma^2 ~ N(0.5, 0.25) * I(sigma^2 > 0)
# Essentially a truncated normal with a
# different scaling constant

# Generate some synthetic data
# in practice, we just have some data
set.seed(7)
y = rnorm(10, mean = 11, sd = 0.47)

# create function for unnormalized posterior
# takes the vector theta = (mu, sigma^2)
# and y = data.
qtheta = function(theta, y) {
  mu = theta[1]
  sigma = sqrt(theta[2])
  log_pdata = sum(dnorm(y, mean = mu, sd = sigma, log = TRUE))
  log_pmu = dunif(mu, 10, 15, log = TRUE)
  log_psigmasq = dnorm(sigma^2, mean = 0.5, sd = 0.5, log = TRUE)
  log_qtheta = log_pdata + log_pmu + log_psigmasq
  exp(log_qtheta)
}

# determine normalizing constant of posterior
(py = cubintegrate(f = qtheta,
                   lower = c(10, 0),
                   upper = c(15, 10),
                   y = y))
const = py$integral

# write posterior function
dpost = function(theta, y, const) {
  qtheta(theta, y) / const
}

# double-check normalizing constant
cubintegrate(f = dpost,
             lower = c(10, 0),
             upper = c(15, 10),
             y = y, const = const)

# mean of marginal posteriors
e_mu = function(theta, y, const) {
  theta[1] * dpost(theta, y, const)
}
e_sigmasq = function(theta, y, const) {
  theta[2] * dpost(theta, y, const)
}

# marginal posterior mean for mu
(mu_mean = cubintegrate(f = e_mu,
             lower = c(10, 0),
             upper = c(15, 10),
             y = y, const = const))
# posterior mean
mu_mean$integral

# marginal posterior mean for sigmasq
(sigmasq_mean = cubintegrate(f = e_sigmasq,
                             lower = c(10, 0),
                             upper = c(15, 10),
                             y = y, const = const))
# posterior mean for sigmasq
sigmasq_mean$integral

# create sequence of values for mu and sigmasq
mymu = seq(10, 15, length = 200)
mysigmasq = seq(0, 2, length = 200)

# double-checking the results
# determine marginal density for mu
# write posterior as function of sigmasq,
# holding mu constant
dpost_sgm = function(sigmasq, mu, y, const) {
  dpost(c(mu, sigmasq), y = y, const = const)
}

# vectorize dpost_sgm
vdpost_sgm = Vectorize(dpost_sgm, vectorize.args = "sigmasq")

# create function that evaluates the marginal
# posterior density of mu after integrating
# out sigmasq
dpostmu = function(mu) {
  integrate(vdpost_sgm, lower = 0, upper = 10,
            mu = mu, y = y, const = const)$value
}

# vectorize dmu
vdpostmu = Vectorize(dpostmu)

# evaluate marginal posterior density of mu
vdpostmu_out = vdpostmu(mymu)

# plot results and compare to estimated
# posterior mean
plot(mymu, vdpostmu_out, type = "l",
     xlab = expression(mu), ylab = expression(p(mu*"|"*y)))
abline(v = 11.04896)
title("marginal posterior distribution of mu")

# double-checking marginal posterior for mu integrates to 1
integrate(vdpostmu, 10, 15)

# double-checking the results
# determine marginal density for sigmasq
# write posterior as function of mu given sigma
dpost_mgs = function(mu, sigmasq, y, const) {
  dpost(c(mu, sigmasq), y = y, const = const)
}

# vectorize dpost_mgs
vdpost_mgs = Vectorize(dpost_mgs, vectorize.args = "mu")

# create function that evaluates the marginal
# posterior density of mu after integrating
# out sigmasq
dpost_sigmasq = function(sigmasq) {
  integrate(vdpost_mgs, lower = 10, upper = 15,
            sigmasq = sigmasq, y = y, const = const)$value
}

# vectorize dsigmasq
vdpost_sigmasq = Vectorize(dpost_sigmasq)

# evaluate marginal posterior density of mu
vdpost_sigmasq_out = vdpost_sigmasq(mysigmasq)

# plot results and compare to estimated
# posterior mean
plot(mysigmasq, vdpost_sigmasq_out, type = "l",
     xlab = expression(sigma^2), ylab = expression(p(sigma^2*"|"*y)))
abline(v = 0.5154897)
title("marginal posterior distribution of mu")

# create sequence of values for mu and sigmasq
mymu = seq(10, 12, length = 200)
mysigmasq = seq(0, 1, length = 200)
# create grid
mytheta = expand.grid(mymu, mysigmasq)
# for each row of mytheta, plug it into lup function
z = apply(mytheta, 1, function(theta) {
  dpost(theta, y = y, const = const)
})
# convert z to matrix for plotting
zmat = matrix(z, nrow = length(mymu))

# create heat map of objective surface
image(mymu, mysigmasq, zmat,
      col = hcl.colors(64, "YlOrRd", rev = TRUE),
      xlab = expression(mu), ylab = expression(sigma^2))
# add contours
contour(mymu, mysigmasq, zmat,
        add = TRUE)

# determine map
optim(par = c(12.5, 1),
      f = dpost,
      lower = c(10.0001, 0.05),
      upper = c(14.9999, 2),
      method = "L-BFGS-B",
      y = y, const = const,
      control = list(fnscale = -1))
# place point for posterior mode
points(map$par[1], map$par[2], pch = 20)
# place point for posterior means
points(mu_mean$integral, sigmasq_mean$integral, pch = 4)
title("log posterior density (unnormalized)")
legend("topright", legend = c("map", "mean"),
       pch = c(20, 4))


