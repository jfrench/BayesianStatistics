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
# Poisson-lognormal model
# Data distribution: Poisson(lambda_i) with
# Y_i ~ log(lambda_i) = beta_0 + beta_1 * gender_i + epsilon_i
# Prior distribution:
# beta_j ~ N(0, 1000), j = 0, 1
# epsilon_i ~ N(0, sigmasq_gender)
# sigmasq_gender ~ IG(0.001, 0.001), gender = 0, 1 

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

# poisson log-normal model
pln_code = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // count of sexual intercourse in previous month
  int<lower=0, upper=1> gender[n];  // male = 0, female = 1
}
parameters {
  vector[n] e;  //epsilon vector
  vector[2] beta; //vector of coefficients
  vector[2] sigmasq_e; //variance of errors
}
transformed parameters {
  vector[n] log_mu;
  for (i in 1:n) {
    log_mu[i] = beta[1] + beta[2]* gender[i] + e[i];
  }
}
model {
  // prior distributions
  // different parameter for each gender
  for (j in 1:2) {
    sigmasq_e[j] ~ inv_gamma(0.001, 0.001);
    beta[j] ~ normal(0, sqrt(1000));
  }
  // data/mixture distributions
  for(i in 1:n) {
    e[i] ~ normal(0, sqrt(sigmasq_e[gender[i] + 1]));
    y[i] ~ poisson_log(log_mu[i]);
  }
}
generated quantities {
  vector[2] lambda_z; //exp of mean of log_mu
  vector[2] v_z; //variance of log_mu for each gender
  vector[2] mu;  //mean for each gender
  vector[2] v; //variance for each gener
  vector[n] log_lik; //log likelihood of observations
  vector[n] yrep;  //replicated data
  for(j in 1:2) {
    lambda_z[j] = exp(beta[1] + beta[2] * (j == 2));
    //v_z[j] = (exp(sigmasq_e[j]) - 1) * exp(2 * mu_z[j] + sigmasq_e[j]);
    mu[j] = lambda_z[j] * exp(sigmasq_e[j]/2);
    v_z[j] = mu[j] + lambda_z[j]^2 * (exp(2 * sigmasq_e[j]) - exp(sigmasq_e[j]));
    v[j] = v_z[j] + mu[j];
  }
  for(i in 1:n) {
    log_lik[i] = poisson_log_lpmf(y[i] | log_mu[i]);
    yrep[i] = poisson_log_rng(log_mu[i]);
  }
}
"
# poisson lognormal starting values
init_fun_pln = function() { 
  list(sigmasq_e = rexp(2), 
       beta = rnorm(2, 0, 1), 
       e = rnorm(gss$n, 0, 0.1)) 
} 

library(rstan)
# pln = stan(model_code = pln_code, data = gss,
#            iter = 5e4, seed = 1004,
#            init = init_fun_pln)
# 
# # summary information
# summary_pln = summary(pln, pars = c("beta", "sigmasq_e", "mu", "v"))
# save(summary_pln, file = "example_9_3_pln_summary.rda")
load(file = "example_9_3_pln_summary.rda")
summary_pln$summary

# ## compute model fit information
# # extract log likelihood
# library(loo)
# ll_pln = extract_log_lik(pln, merge_chains = FALSE)
# # compute waic
# waic_pln = waic(ll_pln)
# # compute effective sample size of liklihood
# r_eff_pln = relative_eff(exp(ll_pln))
# # compute looic (using effective sample size)
# looic_pln = loo(ll_pln, r_eff = r_eff_pln)
# save(waic_pln, looic_pln,
#      file = "example_9_3_pln_ic.rda")
load(file = "example_9_3_pln_ic.rda")
waic_pln
looic_pln

# # extract small subset of yrep
# samples_pln = extract(pln)
# s = sample(seq_len(nrow(samples_pln$log_lik)), 1000)
# yrep_pln = samples_pln$yrep[s, ]
# save(yrep_pln, file = "example_9_3_pln_yrep.rda")
load(file = "example_9_3_pln_yrep.rda")

library(bayesplot)
ppc_hist(gss$y, yrep_pln[1:8, ])
ppc_ecdf_overlay(gss$y, yrep_pln)
ppc_intervals(gss$y, yrep_pln)

# params_pln = as.array(pln, pars = c("beta", "sigmasq_e", "mu", "v"))
# save(params_pln, file = "example_9_3_pln_params.rda")
load(file = "example_9_3_pln_params.rda")

mcmc_trace(params_pln, regex_pars = c("beta", "sigmasq_e"))
mcmc_trace(params_pln, regex_pars = c("mu", "v"))
