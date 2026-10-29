# Soda regression with data-scaled weakly informative priors
# CmdStanR fits and CSV files are saved in stan/ for reuse between runs.
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

# Soda: observed data
soda <- matrix(c(
  16.68, 7, 560,
  11.5, 3, 220,
  12.03, 3, 340,
  14.88, 4, 80,
  13.75, 6, 150,
  18.11, 7, 330,
  8, 2, 110,
  17.83, 7, 210,
  79.24, 30, 1460,
  21.5, 5, 605,
  40.33, 16, 688,
  21, 10, 215,
  13.5, 4, 255,
  19.75, 6, 462,
  24, 9, 448,
  29, 10, 776,
  15.35, 6, 200,
  19, 7, 132,
  9.5, 3, 36,
  35.1, 17, 770,
  17.9, 10, 140,
  52.32, 26, 810,
  18.75, 9, 450,
  19.83, 8, 635,
  10.75, 4, 150),
  ncol = 3, byrow = TRUE)
soda <- as.data.frame(soda)
colnames(soda) <- c("Time", "Cases", "Distance")

# Soda: data-scaled default prior parameters
sy <- sd(soda$Time)
prior_mean <- mean(soda$Time)
prior_intercept_sd <- 2.5 * sy
prior_beta1_sd <- 2.5 * sy / sd(soda$Cases)
prior_beta2_sd <- 2.5 * sy / sd(soda$Distance)
prior_sigma_rate <- 1 / sy

# Soda: center predictors and prepare Stan data
cases_c <- soda$Cases - mean(soda$Cases)
dist_c <- soda$Distance - mean(soda$Distance)
stan_dat_default <- list(
  n = nrow(soda),
  y = soda$Time,
  cases = cases_c,
  dist = dist_c,
  prior_mean = prior_mean,
  prior_intercept_sd = prior_intercept_sd,
  prior_beta1_sd = prior_beta1_sd,
  prior_beta2_sd = prior_beta2_sd,
  prior_sigma_rate = prior_sigma_rate
)

# Soda: Stan program
soda_default_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] cases;
  vector[n] dist;
  real prior_mean;
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_beta1_sd;
  real<lower=0> prior_beta2_sd;
  real<lower=0> prior_sigma_rate;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}
model {
  beta0 ~ normal(prior_mean, prior_intercept_sd);
  beta1 ~ normal(0, prior_beta1_sd);
  beta2 ~ normal(0, prior_beta2_sd);
  sigma ~ exponential(prior_sigma_rate);
  y ~ normal(beta0 + beta1*cases + beta2*dist, sigma);
}
generated quantities {
  real sigmasq = square(sigma);
  vector[n] mu = beta0 + beta1*cases + beta2*dist;
  real<lower=0,upper=1> R2 = variance(mu) / (variance(mu) + sigmasq);
}
"

# Soda: compile and save model
if (!file.exists("stan/soda_default_mod.rda")) {
  soda_default_file <- write_stan_file(soda_default_code, dir = "stan")
  soda_default_mod <- cmdstan_model(soda_default_file)
  save(soda_default_mod, file = "stan/soda_default_mod.rda")
}
load("stan/soda_default_mod.rda")

# Soda: fit and save posterior samples
if (!file.exists("stan/soda_default_fit.rda")) {
  soda_default_fit <- soda_default_mod$sample(
    data = stan_dat_default, seed = 30, chains = 4,
    parallel_chains = 2, iter_warmup = 1000,
    iter_sampling = 2000, refresh = 500,
    output_dir = "stan")
  save(soda_default_fit, file = "stan/soda_default_fit.rda")
}
load("stan/soda_default_fit.rda")

# Posterior summaries
soda_default_fit$summary(
  variables = c("beta0", "beta1", "beta2", "sigmasq"),
  default_summary_measures()
)

# Posterior densities
mcmc_dens(soda_default_fit$draws(
  variables = c("beta0", "beta1", "beta2", "sigmasq")
))

# Convergence measures
soda_default_fit$summary(
  variables = c("beta0", "beta1", "beta2", "sigmasq"),
  default_convergence_measures())

# rank histograms
mcmc_rank_hist(soda_default_fit$draws(
  variables = c("beta0", "beta1", "beta2", "sigmasq")
))

# HMC/NUTS diagnostics
soda_default_fit$diagnostic_summary()
