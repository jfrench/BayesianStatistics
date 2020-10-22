# Create sampler for N(mu, sigma^2) sampling distribution
# with mu and sigma^2 unknown.

# Midge Example
# Chapter 5 of A First Course in Bayesian Statistical Methods by PD Hoff

# Grogan and Wirth (1981) provide data on the wing length in
# millimeters of nine members of a species of midge (small, two-winged
# flies). From these nine measurements we wish to make inference on
# the population mean theta.

# Data distribution: y1, ..., yn | mu, sigma^2 ~ iid N(mu, sigma^2)
# Prior for mu: mu | sigma^2 ~  N(mu0, sigma^2/k0)
# w/ mu0 = 1.9 and k0 = 1
# Prior for sigma^2: sigma^2 ~ Inv-Chisq(nu0, sigma0^2),
# w/ nu = 1 and sigma0 = 0.1.
# Note: this is equivalent to
# sigma^2 ~ Inv-Gamma(nu0/2, nu0/2 * sigma0^2)

# Studies of other populations suggest that the true mean
# should be around 1.9 mm
# with a standard deviation of 0.1.  However, this
# population may be different
# from the others, so we choose k0 and nu0 = 1 so
# that the prior distributions
# are only weakly centered around these estimates
# from other populations.

library(rstan)
library(bayesplot)

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // measurements
  real mu0;
  real<lower=0> k0;
  real<lower=0> nu0;
  real<lower=0> sigma0;
}
parameters {
  real<lower=0> sigmasq;
  real mu;
}
model {
  y ~ normal(mu, sqrt(sigmasq));  // data distribution
  mu ~ normal(mu0, sqrt(sigmasq/k0)); // prior for mu
  //prior for sigmasq
  sigmasq ~ inv_gamma(nu0/2, nu0/2*sigma0^2);
}
generated quantities {
  // create vector to store replicated values
  vector[n] yrep;
  // draw replicated values of each observation
  // from the posterior predictive distribution
  for (j in 1:n) {
    yrep[j] = normal_rng(mu, sqrt(sigmasq));
  }
}
"
# load data
y = c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
stan_dat = list(n = length(y), y = y, mu0 = 1.9, k0 = 1,
                 nu0 = 1, sigma0 = 0.1)

# # midge_fit model using stand with 4 chains
# midge_fit = stan(model_code = stanmod, data = stan_dat, iter = 1e5)
# # compile model
# midge_mod = stan_model(model_code = stanmod)
# # save model
# save(midge_mod, file = "midge_mod.rda", compress = "xz")
load(file = "midge_mod.rda")
# draw samples from the model
midge_fit = sampling(midge_mod, data = stan_dat, iter = 1000, chains = 4)

summary(midge_fit, par = c("mu", "sigmasq"), prob = c(0.025, 0.975))
# Approximate posterior from previous analysis
# > #95% central posterior interval for mu
# > quantile(mu, c(.025, .975))
#       2.5%    97.5%
#   1.727237 1.900704
# >
# > #95% posterior interval for sigmasq
# > quantile(sigmasq, c(.025, .975))
#          2.5%       97.5%
#   0.007488022 0.047175081

posterior = as.array(midge_fit)
mcmc_trace(posterior, pars = c("mu", "sigmasq"))

# density plots of results
mcmc_dens_overlay(posterior, pars = c("mu", "sigmasq"))

# plot of acf of chains
stan_ac(midge_fit, "mu")
stan_ac(midge_fit, "sigmasq")

# a comparison of the errors from y - yrep
yrep = extract(midge_fit, "yrep")$yrep

# histogram comparing y, yrep
ppc_hist(y, yrep[21:28,])

# scatterplot of y vs yrep
ppc_scatter(y, yrep[1:9,])

# scatter of y vs yhat, where yhat is the mean of the
# posterior predictive distribution
ppc_scatter_avg(y, yrep)

# scatter of y vs avg (y - yrep)
ppc_error_scatter_avg(y, yrep)
