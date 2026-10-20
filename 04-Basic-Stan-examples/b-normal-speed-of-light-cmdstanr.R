library(cmdstanr)
library(bayesplot)
library(posterior)

# Example from Bayesian Data Analysis, 3rd edition Speed of
# light example We have 66 observations from an experiment
# to estimate the speed of light. We will make posterior
# inference about mu and sigma^2, the mean and variance of
# the light measurements.

# Data distribution:
# y1, ..., yn | mu, sigma^2 ~ iid N(mu, sigma^2)
# Prior distributions:
# p(mu) propto 1
# p(log(sigma)) propto 1, implying p(sigma^2) propto 1/sigma^2

# Stan parameters are mu and logsigma, each with an improper
# flat prior. The induced prior for sigma^2 is proportional
# to 1/sigma^2.

# Create model. Notice the quotes.
sol_code <- "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // measurements
}
parameters {
  real logsigma;
  real mu;
}
transformed parameters {
  real<lower=0> sigma = exp(logsigma);
}
model {
  y ~ normal(mu, sigma);
  // Flat improper priors for mu and logsigma
}
generated quantities {
  real<lower=0> sigmasq = square(sigma);
  vector[n] yrep;
  for (j in 1:n) {
    yrep[j] = normal_rng(mu, sigma);
  }
}
"

# Load data
data("newcomb", package = "MASS")
y <- newcomb
stan_dat <- list(n = length(y), y = y)

# Compile if the saved model is missing.
if (!file.exists("stan/sol_mod.rda")) {
  sol_file <- write_stan_file(code = sol_code, dir = "rda")
  sol_mod <- cmdstan_model(sol_file)
  save(sol_mod, file = "stan/sol_mod.rda", compress = "xz")
}
load("stan/sol_mod.rda")

# Draw samples from the model
sol_fit <- sol_mod$sample(data = stan_dat,
  seed = 814, chains = 4, parallel_chains = 2,
  iter_warmup = 1000, iter_sampling = 1000, refresh = 500)

# Extract posterior draws in data-frame format. Exclude the
# large yrep array from general parameter diagnostics.
draws <- sol_fit$draws(var = c("mu", "sigmasq"))

# Summarize posterior draws
summarize_draws(draws)

# Rank histograms
mcmc_rank_hist(draws, pars = c("mu", "sigmasq"))

# Density overlays
mcmc_dens_overlay(draws, pars = c("mu", "sigmasq"))

# Exact central 90% posterior intervals for comparison:
# mu | y ~ t_(n-1)(mean(y), sd(y)/sqrt(n))
# sigma^2 | y ~ Inv-Chisq(n-1, var(y))
n <- length(y)
m <- mean(y)
s <- sd(y)
qt(c(0.05, 0.95), df = n - 1) * s / sqrt(n) + m
# Quantiles of Inv-Chisq(n-1, s^2): (n-1)*s^2 / chi-square quantiles
(n - 1) * s^2 / qchisq(c(0.95, 0.05), df = n - 1)

# Posterior predictive checks
# Convert generated yrep variables into a draws-by-observations matrix.
yrep <- as_draws_matrix(sol_fit$draws(variables = "yrep"))

# Histogram comparing observed and replicated data
ppc_hist(y, yrep[51:58, , drop = FALSE])

# Scatterplot comparing observed and replicated data
ppc_scatter(y, yrep[10:18, , drop = FALSE])

# Observed values versus average replicated values
ppc_scatter_avg(y, yrep)