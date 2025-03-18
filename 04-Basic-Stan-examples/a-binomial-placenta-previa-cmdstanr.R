# Load libraries rstan (to use stan) and coda (to analyze results)
# and bayesplot for cool plots
library(cmdstanr)
library(coda)
library(ggplot2)
library(bayesplot)
library(posterior)

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

# Create model.  Notice the quotes
stan_mod <- "
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
stan_dat <- list(n = 980, y = 437, alpha = 1, beta = 1)

# compile and save model if it doesn't already exist
if (!file.exists("pprevia_mod.rda")) {
  # write code to file
  stan_file <- write_stan_file(code = stan_mod)
  # compile model
  pprevia_mod <- cmdstan_model(stan_file)
  # save model
  save(pprevia_mod,
       file = "pprevia_mod.rda",
       compress = "xz")
}

load(file = "pprevia_mod.rda")
# draw samples from the model
pprevia_sample <-
  pprevia_mod$sample(
    data = stan_dat,
    seed = 314,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 500 # print update every 500 iters
)

# extract posterior draws in data.frame format
post_draws <- pprevia_sample$draws(
  variables = c("theta", "ytilde"),
  format = "df"
)

options(tibble.width = Inf)

# summarize posterior draws
summarize_draws(post_draws)

# subset of summary
summarize_draws(post_draws,
                "mean", "sd",
                ~quantile2(., probs = 0.05),
                ~quantile2(., probs = 0.95),
                "rhat", "ess_bulk")
# trace plot
mcmc_trace(post_draws, pars = c("theta"))

# central posterior interval (median shown as point, by default)
# shows 50% and 95% posterior intervals, by default
mcmc_intervals(post_draws, pars = "theta")
# area plot of each parameter (shows 50% credible interval by default)
# shows median as vertical line
mcmc_areas(post_draws, pars = "theta")
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
