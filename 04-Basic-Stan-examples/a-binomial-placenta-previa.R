### Example from Bayesian Data Analysis, 3rd edition
### 2.5 Placenta previa example
# Placenta previa is a condition in which the placenta of an unborn child is
# implanted very low in the uterus, obstructing the child from
# a normal vaginal delivery. An early study concerning
# the sex of placenta previa births in Germany
# found that of a total of 980 births, 437 were female.
# How strong is the evidence that the proportion of female
# births in the population of placenta previa births is less than 0.485
# (the proportion in the general population)?

# Data distribution: y ~ Bin(n, theta) with n = 980 and observed y = 437
# Prior: theta ~ Beta(1, 1)
# Posterior: Beta(y + 1, n - y + 1) = Beta(438, 544)

# Load libraries rstan (to use stan) and coda (to analyze results)
# and bayesplot for cool plots
library(rstan)
library(coda)
library(bayesplot)

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of births
  int<lower=0> y; // number of female placenta previa births
  real<lower=0> alpha; // prior parameters
  real<lower=0> beta; // prior parameters
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  y ~ binomial(n, theta);  //data distribution
  theta ~ beta(alpha, beta); //prior distribution
}
generated quantities{
  int<lower=0> ytilde;
  ytilde = binomial_rng(n, theta);
}
"

# Specify the data in R, using a list format compatible with STAN:
stan_dat = list(n = 980, y = 437, alpha = 1, beta = 1)

# compile the model, sample from model, returns object of class stanfit
fit = stan(model_code = stanmod, data = stan_dat, iter = 1000, chains = 4)
# # alternatively
# # describe model
# stan_mod = stan_model(model_code = stanmod)
# # draw samples from the model
# stan_samples = sampling(stan_mod, data = stan_dat, iter = 1000, chains = 4)

# summary of stanfit object
summary(fit, pars = c("theta", "ytilde"), probs = c(0.025, 0.975))

# rstan plotting functions
# posterior intervals and point estimates
# traceplot of chains
stan_trace(fit)
# density plot of theta
stan_dens(fit, "theta")
# histogram of ytilde
stan_hist(fit, "ytilde")
# ACF plot of chain
stan_ac(fit, "theta")

# coda plots
# convert samples to coda object
codasamples = As.mcmc.list(fit)

# trace plots of posterior samples
coda::traceplot(codasamples)

# density plots of samples
densplot(codasamples)

# summarize codasamples
summary(codasamples, quantiles = c(.025, .975))

# bayesplot plots
posterior = as.array(fit) # convert to array for plotting
# central posterior interval (median shown as point, by default)
# shows 50% and 95% posterior intervals, by default
mcmc_intervals(posterior, pars = "theta")
# area plot of each parameter (shows 50% credible interval by default)
# shows median as vertical line
mcmc_areas(posterior, pars = "theta")
mcmc_areas(posterior, pars = "ytilde")
# histogram
mcmc_hist(posterior, pars = c("theta", "ytilde"))
# histogram by chain
mcmc_hist_by_chain(posterior, pars = c("theta", "ytilde"))
# density plot
mcmc_dens(posterior, pars = c("theta", "ytilde"))
# overlay density for each chain
mcmc_dens_overlay(posterior, pars = c("theta", "ytilde"))
# violin plots
mcmc_violin(posterior, pars = c("theta", "ytilde"))
# trace plots
mcmc_trace(posterior, pars = "theta")
# highlight a single chain
mcmc_trace_highlight(posterior, pars = "theta", highlight = 1)
