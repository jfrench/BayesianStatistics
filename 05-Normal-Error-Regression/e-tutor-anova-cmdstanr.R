# Updated from the accompanying Bayesian linear regression Quarto presentation.
# CmdStanR fits and CSV files are saved in stan/ for reuse between runs.
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

# Tutors: scores and group assignments
grades <- c(84, 58, 100, 51, 28, 89,
			97, 50, 76, 83, 45, 42, 83,
			64, 47, 83, 81, 83, 34, 61,
			77, 69, 94, 80, 55, 79)

tutor <- rep(1:4, times = c(6, 7, 7, 6))

tutor2 <- as.numeric(tutor == 2)
tutor3 <- as.numeric(tutor == 3)
tutor4 <- as.numeric(tutor == 4)

# Tutors: descriptive summaries
tapply(grades, tutor, mean)
boxplot(grades ~ tutor, xlab = "Tutor", ylab = "Test grade")

# Tutors: Stan model
tutor_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] tutor2;
  vector[n] tutor3;
  vector[n] tutor4;
}
parameters {
  real<lower=0> sigmasq;
  real mu0;
  real alpha2;
  real alpha3;
  real alpha4;
}
model {
  vector[n] mu = mu0 + alpha2*tutor2 + alpha3*tutor3 + alpha4*tutor4;
  mu0 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  alpha3 ~ normal(0, 100);
  alpha4 ~ normal(0, 100);
  sigmasq ~ inv_gamma(0.01, 0.01);
  y ~ normal(mu, sqrt(sigmasq));
}
generated quantities {
  vector[4] group_mean;
  group_mean[1] = mu0;
  group_mean[2] = mu0 + alpha2;
  group_mean[3] = mu0 + alpha3;
  group_mean[4] = mu0 + alpha4;
}
"

# Tutors: compile and sample
tutor_data <- list(n = length(grades), y = grades,
  tutor2 = tutor2, tutor3 = tutor3, tutor4 = tutor4)
if (!file.exists("stan/tutor_mod.rda")) {
  stan_file <- write_stan_file(tutor_code, dir = "stan")
  tutor_mod <- cmdstan_model(stan_file)
  save(tutor_mod, file = "stan/tutor_mod.rda")
}
load("stan/tutor_mod.rda")
if (!file.exists("stan/tutor_fit.rda")) {
  tutor_fit <- tutor_mod$sample(data = tutor_data,
    seed = 41, chains = 4, parallel_chains = 2,
    iter_warmup = 2000, iter_sampling = 2000, output_dir = "stan")
  save(tutor_fit, file = "stan/tutor_fit.rda")
}
load("stan/tutor_fit.rda")

# Tutors: posterior summaries
tutor_fit$summary(variables = c("mu0", "alpha2", "alpha3", "alpha4", "sigmasq"))

# Tutors: posterior distributions of group means
bayesplot::mcmc_dens_overlay(tutor_fit$draws(variables = "group_mean"))

# Tutors: calculate posterior probabilities
gm <- posterior::as_draws_matrix(tutor_fit$draws(variables = "group_mean"))
mean(gm[, "group_mean[4]"] > gm[, "group_mean[1]"])
mean(gm[, "group_mean[4]"] > gm[, "group_mean[2]"])
mean(gm[, "group_mean[4]"] > gm[, "group_mean[3]"])
mean(max.col(gm, ties.method = "first") == 4)
# MCMC convergence and HMC/NUTS diagnostics
tutor_fit$summary(variables = c("mu0", "alpha2", "alpha3", "alpha4", "sigmasq"))
bayesplot::mcmc_trace(tutor_fit$draws(variables = c("mu0", "alpha2", "alpha3", "alpha4", "sigmasq")))
bayesplot::mcmc_rank_hist(tutor_fit$draws(variables = c("mu0", "alpha2", "alpha3", "alpha4", "sigmasq")))
tutor_fit$diagnostic_summary()  # divergences, maximum treedepth, E-BFMI
tutor_fit$cmdstan_diagnose()
