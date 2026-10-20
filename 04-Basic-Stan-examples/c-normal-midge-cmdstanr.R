library(cmdstanr)
library(bayesplot)
library(posterior)

# Midge example, Chapter 5 of A First Course in Bayesian Statistical Methods
# by P. D. Hoff.
# Grogan and Wirth (1981) measured wing lengths (mm) for nine midges.
# We want inference about the population mean and variance.
#
# Data distribution: y_i | mu, sigmasq ~ Normal(mu, sigmasq)
# Prior: mu | sigmasq ~ Normal(mu0, sigmasq/k0)
#        sigmasq ~ Inv-Chisq(nu0, sigma0^2)
# The inverse-chi-squared prior is equivalent to
# Inv-Gamma(nu0/2, nu0*sigma0^2/2).
# We use mu0 = 1.9, k0 = 1, nu0 = 1, sigma0 = 0.1.

# Stan model
midge_code <- "
data {
  int<lower=1> n;
  vector[n] y;
  real mu0;
  real<lower=0> k0;
  real<lower=0> nu0;
  real<lower=0> sigma0;
}
parameters {
  real<lower=0> sigmasq;
  real mu;
}
model {
  y ~ normal(mu, sqrt(sigmasq));
  mu ~ normal(mu0, sqrt(sigmasq / k0));
  sigmasq ~ inv_gamma(nu0 / 2, nu0 * square(sigma0) / 2);
}
generated quantities {
  vector[n] yrep;
  for (j in 1:n) {
    yrep[j] = normal_rng(mu, sqrt(sigmasq));
  }
}
"

# Observed data and prior parameters
y <- c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
stan_dat <- list(n = length(y), y = y, mu0 = 1.9,
                 k0 = 1, nu0 = 1, sigma0 = 0.1)

# compile if the saved model is missing.
if (!file.exists("stan/midge_mod.rda")) {
  midge_file <- write_stan_file(code = midge_code, dir = "stan")
  midge_mod <- cmdstan_model(midge_file)
  save(midge_mod, file = "stan/midge_mod.rda", compress = "xz")
}
load("stan/midge_mod.rda")

# Draw samples (1000 warmup + 1000 retained draws per chain)
midge_fit <- midge_mod$sample(
  data = stan_dat, seed = 4, chains = 4, parallel_chains = 2,
  iter_warmup = 1000, iter_sampling = 1000, refresh = 500)

# Posterior draws for model parameters
draws <- midge_fit$draws(var = c("mu", "sigmasq"))

# Diagnostics: R-hat, bulk/tail ESS
summarize_draws(draws)

# rank plots
mcmc_rank_hist(draws, pars = c("mu", "sigmasq"))

# autocorrelation
mcmc_acf(draws, pars = c("mu", "sigmasq"))

# central 95% posterior interval from previous analysis
# mu:      (1.727237, 1.900704)
# sigmasq: (0.007488022, 0.047175081)
summarize_draws(draws, lb = function(x) quantile(x, 0.025),
                ub = function (x) quantile(x, 0.975))

# Posterior predictive checks
# Rows = posterior draws; columns = replicated observations.
yrep <- as_draws_matrix(midge_fit$draws(variables = "yrep"))

# posterior predictive checks
ppc_hist(y, yrep[21:28, , drop = FALSE])
ppc_scatter(y, yrep[1:9, , drop = FALSE])
ppc_scatter_avg(y, yrep)
ppc_error_scatter_avg(y, yrep)