# Updated from the accompanying Bayesian linear regression Quarto presentation.
# CmdStanR fits and CSV files are saved in stan/ for reuse between runs.
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
dir.create("stan", showWarnings = FALSE)

library(loo)

# Oxygen uptake: observations
D2 <- c(0,0,0,0,0,0,1,1,1,1,1,1)
# The age of participant
x <- c(23,22,22,25,27,20,31,23,27,28,22,24)

# The change in maximal oxygen uptake
y <- c(-0.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,0.26,2.51)
# number of observations
n <- length(y)
# sample variance of responses
v <- var(y)

# D2 = 1 for running, 0 for aerobic.
oxygen_data <- list(n = length(y), y = y, x = x, d2 = as.numeric(D2),
  v = var(y))

# Oxygen: parallel-lines Stan program
oxygen_pl_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
  vector[n] d2;
  real<lower=0> v;
}
parameters {
  real beta0;
  real beta1;
  real alpha2;
  real<lower=0> sigmasq;
}
transformed parameters {
  vector[n] mu = beta0 + beta1*x + alpha2*d2;
}
model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  sigmasq ~ inv_gamma(0.01, 0.01);
  y ~ normal(mu, sqrt(sigmasq));
}
generated quantities {
  real Rbsq = 1 - sigmasq / v;
  real<lower=0,upper=1> R2 = variance(mu)/(variance(mu)+sigmasq);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | mu[i], sqrt(sigmasq));
}
"

# Oxygen: separate-lines Stan program
oxygen_sl_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
  vector[n] d2;
  real<lower=0> v;
}
parameters {
  real beta0;
  real beta1;
  real alpha2;
  real delta2;
  real<lower=0> sigmasq;
}
transformed parameters {
  vector[n] mu = beta0 + beta1*x + alpha2*d2 + delta2*(x .* d2);
}
model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  delta2 ~ normal(0, 100);
  sigmasq ~ inv_gamma(0.01, 0.01);
  y ~ normal(mu, sqrt(sigmasq));
}
generated quantities {
  real Rbsq = 1 - sigmasq / v;
  real<lower=0,upper=1> R2 = variance(mu)/(variance(mu)+sigmasq);
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | mu[i], sqrt(sigmasq));
}
"

# Oxygen: compile both models
if (!file.exists("stan/oxygen_pl_mod.rda")) {
  stan_file <- write_stan_file(oxygen_pl_code, dir = "stan")
  oxygen_pl_mod <- cmdstan_model(stan_file)
  save(oxygen_pl_mod, file = "stan/oxygen_pl_mod.rda")
}
load("stan/oxygen_pl_mod.rda")
if (!file.exists("stan/oxygen_sl_mod.rda")) {
  stan_file <- write_stan_file(oxygen_sl_code, dir = "stan")
  oxygen_sl_mod <- cmdstan_model(stan_file)
  save(oxygen_sl_mod, file = "stan/oxygen_sl_mod.rda")
}
load("stan/oxygen_sl_mod.rda")

# Oxygen: fit both models
if (!file.exists("stan/oxygen_pl_fit.rda")) {
  oxygen_pl_fit <- oxygen_pl_mod$sample(data = oxygen_data,
    seed = 43, chains = 4, parallel_chains = 2,
    iter_warmup = 2500, iter_sampling = 2500,
    adapt_delta = 0.99, output_dir = "stan")
  save(oxygen_pl_fit, file = "stan/oxygen_pl_fit.rda")
}
load("stan/oxygen_pl_fit.rda")
if (!file.exists("stan/oxygen_sl_fit.rda")) {
  oxygen_sl_fit <- oxygen_sl_mod$sample(data = oxygen_data,
    seed = 43, chains = 4, parallel_chains = 2,
    iter_warmup = 2500, iter_sampling = 2500,
    adapt_delta = 0.99, output_dir = "stan")
  save(oxygen_sl_fit, file = "stan/oxygen_sl_fit.rda")
}
load("stan/oxygen_sl_fit.rda")

# Oxygen: parallel-lines posterior summaries
oxygen_pl_fit$summary(variables = c("beta0", "beta1", "alpha2", "sigmasq"))

# Oxygen: parallel-lines posterior densities
bayesplot::mcmc_dens_overlay(
  oxygen_pl_fit$draws(variables = c("beta0", "beta1", "alpha2", "sigmasq")))

# Oxygen: fitted parallel lines
pl_means <- oxygen_pl_fit$summary(variables = c("beta0", "beta1", "alpha2"))
pl_b <- setNames(pl_means$mean, pl_means$variable)
plot(x, y, pch = D2 + 1, xlab = "Age", ylab = "Change in oxygen uptake")
legend("topleft", c("aerobic", "running"), pch = 1:2)
abline(pl_b["beta0"], pl_b["beta1"], lty = 1)
abline(pl_b["beta0"] + pl_b["alpha2"], pl_b["beta1"], lty = 2)

# Oxygen: separate-lines posterior summaries
oxygen_sl_fit$summary(variables = c("beta0", "beta1", "alpha2", "delta2", "sigmasq"))

# Oxygen: separate-lines posterior densities
bayesplot::mcmc_dens_overlay(
  oxygen_sl_fit$draws(variables = c("beta0", "beta1", "alpha2", "delta2", "sigmasq")))

# Oxygen: fitted separate lines
sl_means <- oxygen_sl_fit$summary(variables = c("beta0", "beta1", "alpha2", "delta2"))
sl_b <- setNames(sl_means$mean, sl_means$variable)
plot(x, y, pch = D2 + 1, xlab = "Age", ylab = "Change in oxygen uptake")
legend("topleft", c("aerobic", "running"), pch = 1:2)
abline(sl_b["beta0"], sl_b["beta1"], lty = 1)
abline(sl_b["beta0"] + sl_b["alpha2"], sl_b["beta1"] + sl_b["delta2"], lty = 2)

# Oxygen: posterior difference at age 25
od <- posterior::as_draws_matrix(
  oxygen_sl_fit$draws(variables = c("alpha2", "delta2")))
contrast25 <- od[, "alpha2"] + 25 * od[, "delta2"]
quantile(contrast25, c(0.025, 0.5, 0.975))
mean(contrast25 > 0)
# MCMC convergence and HMC/NUTS diagnostics
oxygen_pl_fit$summary(variables = c("beta0", "beta1", "alpha2", "sigmasq"))
bayesplot::mcmc_trace(oxygen_pl_fit$draws(variables = c("beta0", "beta1", "alpha2", "sigmasq")))
bayesplot::mcmc_rank_hist(oxygen_pl_fit$draws(variables = c("beta0", "beta1", "alpha2", "sigmasq")))
oxygen_pl_fit$diagnostic_summary()  # divergences, maximum treedepth, E-BFMI
oxygen_pl_fit$cmdstan_diagnose()

# MCMC convergence and HMC/NUTS diagnostics
oxygen_sl_fit$summary(variables = c("beta0", "beta1", "alpha2", "delta2", "sigmasq"))
bayesplot::mcmc_trace(oxygen_sl_fit$draws(variables = c("beta0", "beta1", "alpha2", "delta2", "sigmasq")))
bayesplot::mcmc_rank_hist(oxygen_sl_fit$draws(variables = c("beta0", "beta1", "alpha2", "delta2", "sigmasq")))
oxygen_sl_fit$diagnostic_summary()  # divergences, maximum treedepth, E-BFMI
oxygen_sl_fit$cmdstan_diagnose()

# WAIC and PSIS-LOO (requires pointwise log_lik in generated quantities)
ll_pl <- oxygen_pl_fit$draws(variables = "log_lik", format = "array")
ll_sl <- oxygen_sl_fit$draws(variables = "log_lik", format = "array")
waic_pl <- loo::waic(ll_pl)
waic_sl <- loo::waic(ll_sl)
loo_pl <- oxygen_pl_fit$loo()
loo_sl <- oxygen_sl_fit$loo()
waic_pl
waic_sl
loo_pl
loo_sl
loo::loo_compare(loo_pl, loo_sl)
loo::loo_compare(waic_pl, waic_sl)
loo::pareto_k_table(loo_pl)
loo::pareto_k_table(loo_sl)
