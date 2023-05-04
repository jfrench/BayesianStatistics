# Example 9.3, Bayesian Modeling Using WinBUGS by Ntzoufras
# We consider modeling the number of times each survey
# participant had intercourse within the previous month.
# The data are accumulated by gender by Agresti (2002, pp.
# 569-70).  The samples means are 5.9 and 4.3 for males and
# females, respectively, while the variances are much higher
# (54.8 and 34.4, respectively).  The data are
# overdispersed!  A simple Poisson log-linear regression
# model would be inadequate.
#
# Poisson-gamma model
# Data distribution: Y_i ~ Poisson(lambda_gender * u)
# Prior distribution:
# u ~ Gamma(r_gender], r_gender])
# r_gender ~ Gamma(0.001, 0.001), gender = 0, 1
# lambda_gender ~ Gamma(0.001, 0.001), gender = 0, 1

# observed data
gss = list(
  n = 550,
  y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12,
        12, 12, 12, 12, 13, 13, 13, 15, 15, 15, 16, 16, 16, 20, 20, 20, 20, 20, 20, 20, 24, 25, 30, 30, 30, 50, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 17, 18, 20, 20, 20, 20, 20, 20, 22, 23, 25, 25, 25, 27, 30),
  gender = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

# poisson-gamma model
pg_code = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // count of sexual intercourse in previous month
  int<lower=0, upper=1> gender[n];  // male = 0, female = 1
}
parameters {
  vector[n] u;
  vector[2] lambda;
  vector[2] r;
}
transformed parameters {
  vector[n] mu;
  for (i in 1:n) {
    mu[i] = lambda[gender[i] + 1] * u[i];
  }
}
model {
  // prior distributions
  // different parameter for each gender
  for (j in 1:2) {
    lambda[j] ~ gamma(0.001, 0.001);
    r[j] ~ gamma(0.001, 0.001);
  }
  // data/mixture distributions
  for(i in 1:n) {
    u[i] ~ gamma(r[gender[i] + 1], r[gender[i] + 1]);
    y[i] ~ poisson(mu[i]);
  }
}
generated quantities {
  vector[n] log_lik; //log likelihood of observations
  vector[n] yrep;  //replicated data
  for(i in 1:n) {
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
    yrep[i] = poisson_rng(mu[i]);
  }
}
"

# the init functions below generate appropriate starting
# values for each chain.  Otherwise, impossible values can
# be generated, causing Stan to fail

# poisson gamma starting values
init_fun_pg = function() {
  list(lambda = rexp(2, 0.1),
       r = runif(2, 0.2, 0.8),
       u = rexp(gss$n, 0.1))
}

library(rstan)
library(loo)
library(bayesplot)

if (!file.exists("example_9_3_pg_output.rda")) {
pg = stan(model_code = pg_code, data = gss,
          iter = 5e1, seed = 99, init = init_fun_pg,
          chains = 4,
          control = list(adapt_delta = 0.9))

# summary information
summary_pg = summary(pg,
                      pars = c("r", "beta", "mu", "v"))
# compute ic
ll_pg <- as.array(pg, pars = "log_lik")
waic_pg = waic(ll_pg)
looic_pg = rstan::loo(pg, pars = "log_lik")

# for plotting
params_pg = as.array(pg,
                      pars = c("r", "beta", "mu", "v"))

# extract small subset of yrep
samples_pg = extract(pg)
s = sample(seq_len(nrow(samples_pg$log_lik)), 1000)
yrep_pg = samples_pg$yrep[s, ]

save(summary_pg,
     waic_pg, looic_pg,
     params_pg,
     yrep_pg,
     file = "example_9_3_pg_output.rda")
}
load(file = "example_9_3_pg_output.rda")

# summary of posteriors
summary_pg$summary

# ic results
waic_pg
looic_pg

# posterior predictive checks
ppc_hist(gss$y, yrep_pg[1:8, ])
ppc_ecdf_overlay(gss$y, yrep_pg)
ppc_intervals(gss$y, yrep_pg)

# trace plots
mcmc_trace(params_pg, regex_pars = c("r", "beta"))
mcmc_trace(params_pg, regex_pars = c("mu", "v"))