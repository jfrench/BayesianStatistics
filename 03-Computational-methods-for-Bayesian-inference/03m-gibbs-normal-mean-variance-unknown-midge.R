### Create Gibbs sampler for N(mu, sigma^2) sampling distribution
### with mu and sigma^2 unknown.
### Midge Example
### Taken from Chapter 5 of A First Course in Bayesian Statistical Methods by
### Peter Hoff

# Grogan and Wirth (1981) provide data on the wing length in
# millimeters of nine members of a species of midge (small, two-winged
# flies). From these nine measurements we wish to make inference on
# the population mean theta.

# Data distribution: y_1, ..., y_n| mu, sigma^2 ~ iid N(mu, sigma^2)
# Prior for mu: mu | sigma^2, mu0, kappa0 ~  N(mu0, sigma^2/kapp0)
# w/ mu0 = 1.9 and kappa0 = 1
# Prior for sigma^2: sigma^2 | nu0, sigma0^2 ~ Inv-Chisq(nu0, sigma0^2),
# w/ nu = 1 and sigma0 = 0.1.

# Studies of other populations suggest that the true mean should be
# around 1.9 mm with a standard deviation of 0.1.  However, this
# population may be different from the others, so we choose k0 = vu0 =
# 1 so that the prior distributions are only weakly centered around
# these estimates from other populations.

# load useful distributions
source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")

set.seed(90)

#Set prior parameters
mu0 = 1.9
k0 = 1
nu0 = 1
sigma0 = 0.1

#define number of simulations, data, sample size, sample st dev, sample mean
B = 500000
y = c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
n = length(y)
s = sd(y)
ybar = mean(y)

#Initial value for sigma
sigma = sd(y)
#Initial value for mu
mu = mean(y)

#parameters for posterior
kn = k0 + n
mun = (k0 * mu0 + n * ybar) / kn
nun = nu0 + n + 1

sigmasqpost = numeric(B)
mupost = numeric(B)
for (i in 1:B) {
	mu = rnorm(1, mun, sigma/sqrt(kn))

	#parameter for full conditional posterior of sigma^2
	ssqn = (nu0 * sigma0^2 + k0 * (mu - mu0)^2 + (n - 1) * s^2 + n * (ybar - mu)^2)/nun

	sigmasq = rinvchisq(1, df = nun, scale = sqrt(ssqn))
	sigma = sqrt(sigmasq)

	mupost[i] = mu
	sigmasqpost[i] = sigmasq
}

# Bayesian Data Analysis, 3e by Gelman et al. provides the exact marginal
# posteriors for mu and sigmasq
# the marginal posterior for mu is
# p(mu | y) ~ t_vn(mu_n, taunsq/kn)
# the marginal posterior for sigmasq is
# p(sigmasq | y) ~ Inv-Chisq(vn, taunsq) with
vn = nu0 + n
taunsq = (nu0*sigma0^2 + (n - 1) * s^2 + k0 * n/kn * (ybar - mu0)^2)/vn

# plot approximate posterior density for mu
plot(density(mupost), main = "", xlab = "mu", xlim = c(1.6, 2))
# plot true posterior density for mu
x = seq(-1, 3, len = 1001)
lines(x, dst(museq, df = nun, mean = mun, sd = sqrt(taunsq/kn)),
      col = "orange")
title("Posterior Density for mu")
legend("topright", legend = c("Gibbs", "true"), col = c("black", "orange"),
	lwd = 1)

# plot marginal posterior for sigmasq
plot(density(sigmasqpost), main = "", xlab = "sigmasq",
	xlim = c(0, 0.1))
x = seq(0, 0.1, len = 1001)
lines(x, dinvchisq(x, df = vn, scale = sqrt(taunsq)), col = "orange")
title("Posterior Density for sigmasq")
legend("topright", legend = c("Gibbs", "True"), col = c("black", "orange"),
	lwd = 1)

#posterior quantiles for mu
p = c(0.01, 0.10, 0.25, 0.5, 0.75, 0.90, 0.99)
quantile(mupost, prob = p)
qst(p, df = nun, mean = mun, sd = sqrt(taunsq/kn))

#posterior quantiles for sigmasq
quantile(sigmasqpost, prob = p)
qinvchisq(p, df = vn, scale = sqrt(taunsq))

#95% posterior interval for mu
quantile(mupost, c(.025, .975))
qst(c(0.025, 0.975), df = nun, mean = mun, sd = sqrt(taunsq/kn))

#95% posterior interval for sigmasq
quantile(sigmasqpost, c(.025, .975))
qinvchisq(c(0.025, 0.975), df = vn, scale = sqrt(taunsq))

#plot path of gibbs sampler for 100 iterations
plot.mcmc.path(cbind(mupost, sigmasqpost)[1:100,],
               xlab = expression(mu), ylab = expression(sigma^2))

#plot path of cycles for 100 iterations
plot(cbind(mupost, sigmasqpost)[1:100,],
     xlab = expression(mu), ylab = expression(sigma^2), type = "l")