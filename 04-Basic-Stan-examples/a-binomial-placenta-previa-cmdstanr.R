library(cmdstanr)
library(bayesplot)
library(posterior)

### Example from Bayesian Data Analysis, 3rd edition
### 2.5 Placenta previa example
# Placenta previa is a condition in which the placenta of an
# unborn child is implanted very low in the uterus,
# obstructing the child from a normal vaginal delivery. An
# early study concerning the sex of placenta previa births
# in Germany found that of a total of 980 births, 437 were
# female. How strong is the evidence that the proportion of
# female births in the population of placenta previa births
# is less than 0.485 (the proportion in the general
# population)?

# Data distribution: y ~ Bin(n, theta) with n = 980, y = 437
# Prior: theta ~ Beta(1, 1)
# Posterior: Beta(y + 1, n - y + 1) = Beta(438, 544)

# Create model.  Notice the quotes
birth_code <- "
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

# Specify the data in R in Stan format
stan_dat <- list(n = 980, y = 437, alpha = 1, beta = 1)

# compile and save model if it doesn't already exist
if (!file.exists("stan/birth_mod.rda")) {
  # write code to file
  birth_file <- write_stan_file(code = birth_code, dir = "stan")
  # compile model
  birth_mod <- cmdstan_model(birth_file)
  # save model
  save(birth_mod, file = "stan/birth_mod.rda", compress = "xz")
}
load(file = "stan/birth_mod.rda")

# draw samples from the model
birth_sample <-
  birth_mod$sample(data = stan_dat, seed = 18, chains = 4,
    parallel_chains = 2, iter_warmup = 1000,
    iter_sampling = 1000, refresh = 500)

# extract posterior draws in data.frame format
draws <- birth_sample$draws(var = c("theta", "ytilde"))

# summarize posterior draws
summarize_draws(draws)

# trace plot
mcmc_trace(draws, pars = c("theta"))
# central posterior interval (median shown as point, by default)
mcmc_intervals(draws, pars = "theta")
# area plot of each parameter (shows 50% credible interval by default)
mcmc_areas(draws, pars = c("theta"))
mcmc_areas(draws, pars = c("ytilde"))
# histogram
mcmc_hist(draws, pars = c("theta", "ytilde"))
# histogram by chain
mcmc_hist_by_chain(draws, pars = c("theta", "ytilde"))
# density plot
mcmc_dens(draws, pars = c("theta", "ytilde"))
# overlay density for each chain
mcmc_dens_overlay(draws, pars = c("theta", "ytilde"))
# violin plots
mcmc_violin(draws, pars = c("theta", "ytilde"))