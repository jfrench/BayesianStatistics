library(cmdstanr)
library(posterior)
library(bayesplot)

# Chapter 7 Example from A First Course in Bayesian
# Statistical Methods by P. Hoff.
# A sample of 22 children are given reading comprehension
# tests before and after receiving a particular
# instructional method.  Each student will then have two
# scores denoting the pre- and post-instructional scores
# Each row of the data set is a two-dimensional vector
# representing a single case, distributed N(mu, Sigma).

# Data distribution:
# y_1, ..., y_n | mu, Sigma ~ i.i.d. N_2(mu, Sigma)
# For priors:
# mu | Sigma ~ N_2(mu0, Sigma/k0)
# Sigma ~ Inv-Wishart(nu, L0^(-1))

# The exam was designed to give average score around 50 out
# of 100, so set mu0 = c(50, 50).  Since the true mean
# cannot be below 0 or above 100, we want the prior
# variances to keep mu inside this range with high
# probability. Assume k0 = 1 (a single prior observation on
# the Sigma scale). We want Sigma[1,1] = Sigma[2,2] to be
# approximately 625, so that the prior probability that the
# mean is outside [0, 100] is only 0.05.  Since the two
# exams measure similar things, we think the correlation
# between the means is around 0.5, so Sigma[1,2] =
# Sigma[2,1] = 0.5(625) = 312.5. Since E(Sigma) ~
# Inv-Wishart(nu0, L0^(-1)) = L0/(nu0 - d - 1), choosing nu0
# = 4 loosely centers samples of Sigma around L0. So we
# choose L0[1,1]=L0[2,2] = 625 and the diagonal elements
# equal to 312.5.

# load data
reading <- matrix(c(
  59, 77, 43, 39, 34, 46, 32, 26, 42, 38, 38, 43, 55, 68,
  67, 86, 64, 77, 45, 60, 49, 50, 72, 59, 34, 38, 70, 48,
  34, 55, 50, 58, 41, 54, 52, 60, 60, 75, 34, 47, 28, 48,
  35, 33), ncol = 2, byrow = TRUE)
# reformat
reading <- data.frame(pretest = reading[,1], posttest = reading[,2])

# Create model.  Notice the quotes
reading_code <- "
data {
  int<lower=1> n; // number of trials
  array[n] vector[2] y; // 2d array of size n
  vector[2] mu0; // prior mean for mu
  real<lower=0> k0; // prior number of observations for mu
  real<lower=0> nu0; // df for prior for Sigma
  cov_matrix[2] L0; // scale of prior for Sigma
}
parameters {
  vector[2] mu;
  cov_matrix[2] Sigma;
}
model {
  // specify my prior distributions for mu and Sigma
  mu ~ multi_normal(mu0, Sigma/k0);
  Sigma ~ inv_wishart(nu0, L0);
  // specify the data distribution for each observation
  for(i in 1:n){
    y[i] ~ multi_normal(mu, Sigma);
  }
}
generated quantities {
  vector[2] ytilde; //pre- and post-test measurements for new students
  ytilde = multi_normal_rng(mu, Sigma);
}
"

L0 <- matrix(c(625, 312.5, 312.5, 625), nrow = 2)

reading_data <- list(n = 22, y = as.matrix(reading),
                 mu0 = c(50, 50), k0 = 1, nu0 = 4, L0 = L0)

# compile and save model if it doesn't already exist
if (!file.exists("stan/reading_mod.rda")) {
  stan_file <- write_stan_file(code = reading_code, dir = "stan")
  reading_mod <- cmdstan_model(stan_file)
  save(reading_mod, file = "stan/reading_mod.rda", compress = "xz")
}
load(file = "stan/reading_mod.rda")

# Draw samples once and reuse the saved fit
if (!file.exists("stan/reading_fit.rda")) {
  reading_fit <- reading_mod$sample(data = reading_data,
    seed = 18, chains = 4, parallel_chains = 2,
    iter_warmup = 1000, iter_sampling = 1000, refresh = 500,
    output_dir = "stan")
  save(reading_fit, file = "stan/reading_fit.rda")
}
load("stan/reading_fit.rda")
reading_draws <- reading_fit$draws(var = c("mu", "Sigma"))

### quantiles of mu from original example
#           1%      25%      50%      75%      99%
# mu[1] 40.21531 45.37341 47.30597 49.26025 54.17412
# mu[2] 45.90558 51.56652 53.72375 55.86592 61.38130

# results should be similar
# posterior summaries for the two means
summarize_draws(
  reading_fit$draws(variables = "mu"),
  mean, median, sd,
  q1 = function(x) quantile(x, 0.01),
  q25 = function(x) quantile(x, 0.25),
  q50 = function(x) quantile(x, 0.50),
  q75 = function(x) quantile(x, 0.75),
  q99 = function(x) quantile(x, 0.99)
)

# extract posterior samples as a matrix
reading_draws_matrix <- as_draws_matrix(
  reading_fit$draws(variables = c("mu", "ytilde")))

# determine posterior probability that
# post-instruction test score mean greater than pre-instruction test score mean
mean(reading_draws_matrix[, "mu[2]"] > reading_draws_matrix[, "mu[1]"])

# determine probability that post-instruction test score
# greater than pre-instruction test score for new student
mean(reading_draws_matrix[, "ytilde[2]"] > reading_draws_matrix[, "ytilde[1]"])

# Reading comprehension: graphical diagnostics
mcmc_trace(reading_draws, pars = c("mu[1]", "mu[2]"))
mcmc_dens_overlay(reading_draws, pars = c("mu[1]", "mu[2]"))

# Joint posterior: plot a subset to avoid overplotting
reading_mu <- reading_fit$draws(variables = "mu", format = "matrix")
set.seed(318)
idx <- sample(seq_len(nrow(reading_mu)),
              size = min(1000, nrow(reading_mu)))
plot(reading_mu[idx, "mu[1]"], reading_mu[idx, "mu[2]"],
     pch = 20, col = adjustcolor("black", alpha.f = 0.35),
     xlab = expression(mu[1]), ylab = expression(mu[2]))

# HMC/NUTS diagnostics
reading_fit$diagnostic_summary()
reading_fit$cmdstan_diagnose()

np <- nuts_params(reading_fit)
lp <- log_posterior(reading_fit)

mcmc_nuts_divergence(np, lp)
reading_fit$diagnostic_summary()$num_max_treedepth
mcmc_nuts_treedepth(np, lp)
mcmc_nuts_energy(np)

mcmc_pairs(reading_fit$draws(variables = "mu", format = "array"),
  pars = c("mu[1]", "mu[2]"), np = np)

# Optional second fit with higher target acceptance probability
if (!file.exists("stan/reading_fit2.rda")) {
  reading_fit2 <- reading_mod$sample(
    data = reading_data,
    seed = 22,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = 5000,
    iter_sampling = 5000,
    adapt_delta = 0.95,
    refresh = 500,
    output_dir = "stan"
  )
  save(reading_fit2, file = "stan/reading_fit2.rda")
}
load("stan/reading_fit2.rda")
reading_fit2$diagnostic_summary()
