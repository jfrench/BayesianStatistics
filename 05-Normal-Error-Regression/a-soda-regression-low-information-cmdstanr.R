# Updated from the accompanying Bayesian linear regression Quarto presentation.
# CmdStanR fits and CSV files are saved in stan/ for reuse between runs.
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

# Soda: observed data
soda <- matrix(c(
  16.68,	7,	560,
  11.5,	3,	220,
  12.03,	3,	340,
  14.88,	4,	80,
  13.75,	6,	150,
  18.11,	7,	330,
  8,	2,	110,
  17.83,	7,	210,
  79.24,	30,	1460,
  21.5,	5,	605,
  40.33,	16,	688,
  21,	10,	215,
  13.5,	4,	255,
  19.75,	6,	462,
  24,	9,	448,
  29,	10,	776,
  15.35,	6,	200,
  19,	7,	132,
  9.5,	3,	36,
  35.1,	17,	770,
  17.9,	10,	140,
  52.32,	26,	810,
  18.75,	9,	450,
  19.83,	8,	635,
  10.75,	4,	150),
  ncol = 3, byrow = TRUE)
soda <- as.data.frame(soda)
colnames(soda) <- c("Time", "Cases", "Distance")

X_soda <- cbind(1, soda$Cases, soda$Distance)
n_soda <- nrow(soda)
v_soda <- var(soda$Time)

# Soda: Stan program
soda_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] cases;
  vector[n] dist;
  real<lower=0> v;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> prec;
}
transformed parameters {
  real<lower=0> sigma = inv_sqrt(prec);
}
model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  prec ~ gamma(0.01, 0.01);
  y ~ normal(beta0 + beta1 * cases + beta2 * dist, sigma);
}
generated quantities {
  real<lower=0> sigmasq = 1 / prec;
  real Rbsq = 1 - sigmasq / v;
  vector[n] mu = beta0 + beta1 * cases + beta2 * dist;
  real<lower=0,upper=1> R2 = variance(mu) / (variance(mu) + sigmasq);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
}
"

# Soda: compile and save model
soda_data <- list(n = n_soda, y = soda$Time,
  cases = soda$Cases, dist = soda$Distance, v = v_soda)
if (!file.exists("stan/soda_mod.rda")) {
  stan_file <- write_stan_file(soda_code, dir = "stan")
  soda_mod <- cmdstan_model(stan_file)
  save(soda_mod, file = "stan/soda_mod.rda")
}
load("stan/soda_mod.rda")

# Soda: fit the regression
if (!file.exists("stan/soda_fit.rda")) {
  soda_fit <- soda_mod$sample(
    data = soda_data, seed = 30, chains = 4,
    parallel_chains = 2, iter_warmup = 2000,
    iter_sampling = 2000, refresh = 500,
    output_dir = "stan"
  )
  save(soda_fit, file = "stan/soda_fit.rda")
}
load("stan/soda_fit.rda")

# Soda: posterior summaries
soda_fit$summary(variables = c("beta0", "beta1", "beta2", "sigmasq"))

# Soda: 95% central posterior intervals
posterior::summarise_draws(
  soda_fit$draws(variables = c("beta0", "beta1", "beta2", "sigmasq")),
  mean,
  q2.5 = function(x) quantile(x, 0.025),
  q97.5 = function(x) quantile(x, 0.975)
)

# Soda: posterior densities
bayesplot::mcmc_dens_overlay(
  soda_fit$draws(variables = c("beta0", "beta1", "beta2", "sigmasq"))
)

# Soda: posterior probability of positive slopes
soda_draws <- posterior::as_draws_matrix(
  soda_fit$draws(variables = c("beta0", "beta1", "beta2", "R2")))
mean(soda_draws[, "beta1"] > 0)
mean(soda_draws[, "beta2"] > 0)

# Soda: posterior of Bayesian $R^2$
bayesplot::mcmc_dens(soda_fit$draws(variables = "R2"))

# MCMC convergence and HMC/NUTS diagnostics
soda_fit$summary(variables = c("beta0", "beta1", "beta2", "sigmasq"))
mcmc_trace(soda_fit$draws(variables = c("beta0", "beta1", "beta2", "sigmasq")))
mcmc_rank_hist(soda_fit$draws(variables = c("beta0", "beta1", "beta2", "sigmasq")))
soda_fit$diagnostic_summary()  # divergences, maximum treedepth, E-BFMI
soda_fit$cmdstan_diagnose()
