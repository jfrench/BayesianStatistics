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
# Negative-Binomial model
# Data distribution: Y_i ~ NegBinomial2(lambda_i, r_i) with
# log(lambda_i) = beta_0 + beta_1 * gender_i,
# r_i = r_gender_i
# Prior distribution:
# beta_j ~ N(0, 1000), j = 0, 1
# r_gender ~ Gamma(0.001, 0.001), gender = 0, 1 

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

# negative binomial model
nb_glm_code = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // count of sexual intercourse in previous month
  int<lower=0, upper=1> gender[n];  // male = 0, female = 1
}
parameters {
  vector[2] beta; //vector of coefficients
  vector[2] r; //related to dispersion
}
transformed parameters {
  vector[n] lambda;
  for (i in 1:n) {
    lambda[i] = exp(beta[1] + beta[2]* gender[i]);
  }
}
model {
  // prior distributions
  // different parameter for each gender
  for (j in 1:2) {
    r[j] ~ gamma(0.001, 0.001);
    beta[j] ~ normal(0, sqrt(1000));
  }
  // data distribution
  for(i in 1:n) {
    y[i] ~ neg_binomial_2(lambda[i], r[gender[i] + 1]);
  }
}
generated quantities {
  vector[2] mu;  //mean for each gender
  vector[2] v;   //variance for each gender
  vector[n] log_lik; //log likelihood of observations
  vector[n] yrep;  //replicated data
  mu[1] = exp(beta[1]);
  mu[2] = exp(beta[1] + beta[2]);
  for (j in 1:2) {
    v[j] = lambda[j] + lambda[j]^2/r[j];
  }
  for(i in 1:n) {
    log_lik[i] = neg_binomial_2_lpmf(y[i] | lambda[i], r[gender[i] + 1]);
    yrep[i] = neg_binomial_2_rng(lambda[i], r[gender[i] + 1]);
  }
}
"

# the init functions below generate appropriate starting
# values for each chain.  Otherwise, impossible values can 
# be generated, causing Stan to fail
# negative binomial starting values
init_fun_nb = function() { 
  list(r = rexp(2), beta = rnorm(2, 0, 3)) 
} 

library(rstan)
# nb_glm = stan(model_code = nb_glm_code, data = gss,
#                iter = 5e4, seed = 7,
#                init = init_fun_nb,
#                control = list(adapt_delta = 0.9))
# 
# # summary information
# summary_nb_glm = summary(nb_glm,
#                          pars = c("r", "beta", "mu", "v"))
# save(summary_nb_glm, file = "example_9_3_nb_glm_summary.rda")
load(file = "example_9_3_nb_glm_summary.rda")
summary_nb_glm$summary

## compute model fit information
# extract log likelihood
library(loo)
# ll_nb_glm = extract_log_lik(nb_glm, merge_chains = FALSE)
# # compute waic
# waic_nb_glm = waic(ll_nb_glm)
# # compute effective sample size of liklihood
# r_eff_nb_glm = relative_eff(exp(ll_nb_glm))
# # compute looic (using effective sample size)
# looic_nb_glm = loo(ll_nb_glm, r_eff = r_eff_nb_glm)
# save(waic_nb_glm, looic_nb_glm,
#      file = "example_9_3_nb_glm_ic.rda")
load(file = "example_9_3_nb_glm_ic.rda")
waic_nb_glm
looic_nb_glm

# # extract small subset of yrep
# samples_nb_glm = extract(nb_glm)
# s = sample(seq_len(nrow(samples_nb_glm$log_lik)), 1000)
# yrep_nb_glm = samples_nb_glm$yrep[s, ]
# save(yrep_nb_glm, file = "example_9_3_nb_glm_yrep.rda")
load(file = "example_9_3_nb_glm_yrep.rda")

library(bayesplot)
ppc_hist(gss$y, yrep_nb_glm[1:8, ])
ppc_ecdf_overlay(gss$y, yrep_nb_glm)
ppc_intervals(gss$y, yrep_nb_glm)

# params_nb_glm = as.array(nb_glm,
#                          pars = c("r", "beta", "mu", "v"))
# save(params_nb_glm, file = "example_9_3_nb_glm_params.rda")
load(file = "example_9_3_nb_glm_params.rda")

# trace plots
mcmc_trace(params_nb_glm, regex_pars = c("r", "beta"))
mcmc_trace(params_nb_glm, regex_pars = c("mu", "v"))
