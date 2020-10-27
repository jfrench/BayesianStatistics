library(rstan)
library(coda)
library(bayesplot)

# Example taken from Bayesian Data Analysis, 3rd edition
### Speed of light example
# We have 66 observations from an experiment to estimate the speed of light.
# We will make posterior inference about mu and sigma^2, the mean and
# variance of the light measurements.

# Data distribution:  y1, ..., yn | mu, sigma^2 ~ iid N(mu, sigma^2)
# Prior distributions: p(mu) propto 1
# p(log(sigma)) propto 1, implying p(sigma^2) propto 1/sigma^2

# Note: If the user doesn't specify a prior distribution in Stan,
# then Stan assumes a uniform prior over the valid range of the
# parameter.

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // measurements
}
parameters {
  real logsigma;
  real mu;
}
transformed parameters {
  real<lower=0> sigma;
  sigma = exp(logsigma);
}
model {
  // specify data distribution
  y ~ normal(mu, sigma);
  // no priors for mu and logsigma means that mu and
  // logsigma have unnormalized uniform priors over the real line
  // this corresponds to p(sigmasq) propto 1/sigmasq
}
generated quantities {
  //variables must be the first things declared in a block
  real<lower=0> sigmasq;
  vector[n] yrep;
  sigmasq = square(exp(logsigma));
  for (j in 1:n) {
    yrep[j] = normal_rng(mu, sigma);
  }
}
"

# load data
data(newcomb, package = "MASS") # load data
stan_dat = list(n = length(newcomb), y = newcomb)

# # fit model using stan with 4 chains
# sol_fit = stan(model_code = stanmod, data = stan_dat, iter = 100000, chains = 4)
# # compile model
# sol_mod = stan_model(model_code = stanmod)
# # save model
# save(sol_mod, file = "sol_mod.rda", compress = "xz")
load(file = "sol_mod.rda")
# draw samples from the model
sol_fit = sampling(sol_mod, data = stan_dat, iter = 1000, chains = 4)

# trace plots of results
posterior = as.array(sol_fit)
mcmc_trace(posterior, pars = c("mu", "sigmasq"))

# density plots of results
mcmc_dens_overlay(posterior, pars = c("mu", "sigmasq"))

# plot of acf of chains
stan_ac(sol_fit, "mu")
stan_ac(sol_fit, "sigmasq")

# summary of chains
summary(sol_fit, pars = c("mu", "sigmasq"),
        prob = c(0.025, 0.975))$summary
# exact 95% posterior intervals for mu and sigma^2, respectively
# > qst(c(.025, .975), df  = n - 1, mean = m, sd = s/sqrt(n))
# [1] 23.57059 28.85365
#  > ### 95% central posterior interval for sigma^2
#  > qinvchisq(c(.025, .975), df = n - 1, scale = s^2)
# [1]  84.15867 168.26293

# compare the empirical distribution of the data y
# to the distributions of simulated/replicated data
# yrep from the posterior predictive distribution
# important for model checking!

y = newcomb
yrep = extract(sol_fit, "yrep")$yrep
# histogram comparing y, yrep
ppc_hist(y, yrep[31:38,])

# scatterplot of y vs yrep
ppc_scatter(y, yrep[10:18,])

# scatter of y vs yhat, where yhat is the mean of the
# posterior predictive distribution
ppc_scatter_avg(y, yrep)

# scatter of y vs avg (y - yrep)
ppc_error_scatter_avg(y, yrep)
