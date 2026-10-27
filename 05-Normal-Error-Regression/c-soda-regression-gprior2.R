# Updated from the accompanying Bayesian linear regression Quarto presentation.
# CmdStanR fits and CSV files are saved in stan/ for reuse between runs.
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
dir.create("stan", showWarnings = FALSE)

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

# Soda $g$-prior: Stan model
soda_g_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  matrix[n,3] X;
  vector[3] mu0;
  cov_matrix[3] V;
  real<lower=0> g;
  real<lower=0> v;
}
parameters {
  vector[3] beta;
  real<lower=0> prec;
}
transformed parameters {
  real<lower=0> sigmasq = 1 / prec;
}
model {
  beta ~ multi_normal(mu0, sigmasq * g * V);
  prec ~ gamma(0.01, 0.01);
  y ~ normal(X * beta, sqrt(sigmasq));
}
generated quantities {
  real Rbsq = 1 - sigmasq / v;
  vector[n] mu = X * beta;
  real<lower=0,upper=1> R2 = variance(mu) / (variance(mu) + sigmasq);
}
"

# Soda $g$-prior: compile once
if (!file.exists("stan/soda_g_mod.rda")) {
  stan_file <- write_stan_file(soda_g_code, dir = "stan")
  soda_g_mod <- cmdstan_model(stan_file)
  save(soda_g_mod, file = "stan/soda_g_mod.rda")
}
load("stan/soda_g_mod.rda")
soda_g_data <- list(n = n_soda, y = soda$Time,
  X = X_soda, mu0 = rep(0, 3),
  V = solve(crossprod(X_soda)), v = v_soda, g = n_soda)

soda_g_data$g <- 100^2
if (!file.exists("stan/soda_g2_fit.rda")) {
  soda_g2_fit <- soda_g_mod$sample(data = soda_g_data, seed = 31,
    chains = 4, parallel_chains = 2, iter_warmup = 2000,
    iter_sampling = 2000, output_dir = "stan")
  save(soda_g2_fit, file = "stan/soda_g2_fit.rda")
}
load("stan/soda_g2_fit.rda")

soda_g2_fit$summary(variables = c("beta", "sigmasq", "R2"))
bayesplot::mcmc_dens_overlay(soda_g2_fit$draws(variables = c("beta", "sigmasq")))

# MCMC convergence and HMC/NUTS diagnostics
soda_g2_fit$summary(variables = c("beta", "sigmasq"))
bayesplot::mcmc_trace(soda_g2_fit$draws(variables = c("beta", "sigmasq")))
bayesplot::mcmc_rank_hist(soda_g2_fit$draws(variables = c("beta", "sigmasq")))
soda_g2_fit$diagnostic_summary()  # divergences, maximum treedepth, E-BFMI
soda_g2_fit$cmdstan_diagnose()
