# Example:  Example 1 of Chapter 9 of Bayesian Modeling
# Using Winbugs by Ntzoufras
#
# Consider a problem related to blood pressure measurements.  We consider
# the repeated measurements of blood pressure from 20 health individuals.
# Our aim is to estimate both the between-subject and within-subject variability.
#
# Random effects model
# We assume that is y_ij = mu + a_i + epsilon_ij
# We consider a hierarchical model
# a_i ~ N(0, sigmasq_a), epsilon_ij ~ N(0, sigmasq)
#
# Equivalently, Y_ij ~ N(mu_i, sigma^2) with mui = mu + a_i
# and a_i ~ N(0, sigmasq_a)
#
# We assume low information prior distributions for the regression coefficients,
# mu ~ N(0, 1000), sigmasq ~ IG(0.001, 0.001), sigmasq_a ~ IG(0.001, 0.001)
#
# Fixed effects model
# Data distribution:
# Y_ij ~ N(mu_i, sigma^2) with mui = mu + a_i
# Prior distributions: a_i ~ N(0, 1000),
# mu ~ N(0, 1000), sigmasq ~ IG(0.001, 0.001)
#
# Common mean model
# Data distribution:
# Y_ij ~ N(mu, sigma^2)
# Prior distributions:
# mu ~ N(0, 1000), sigmasq ~ IG(0.001, 0.001)

library(rstan)

# Enter data manually
n = 20
y = matrix(c(108, NA, 91, NA, NA, NA, 104, NA, 99, 97, 95,
             98, 93, 97, 99, 96, 90, 100, 92, 95, 101, 89,
             97, 97, 97, 100, 96, 95, 106, 100, 100, 98,
             90, 99, 88, 98, 92, 92, 100, 101),
           ncol = 2, byrow = TRUE)

# convert data to vector
y = c(y)
# determine observed and missing indices, other related
# information
obs = which(!is.na(y))
miss = which(is.na(y))
nobs = length(obs)
nmiss = length(miss)
subject_id = rep(1:20, times = 2)
subject_id_obs = subject_id[obs]
subject_id_miss = subject_id[miss]
measurement = rep(1:2, each = 20)
measurement_obs = measurement[obs]
measurement_miss = measurement[miss]
yobs = y[obs]
nsubjects = 20

# random effects model
recode = "
data {
  int<lower=1> nobs; // number of observed data values
  int<lower=1> nmiss;  // number of missing data values
  int<lower=1> nsubjects; // number of subjects measured
  int subject_id_obs[nobs]; //subject index for yobs.  Must be int!
  int subject_id_miss[nmiss]; //subject index for ymiss.  Must be int!
  vector[nobs] yobs;  // vector of observed blood pressure measurements
}
parameters {
  real mu;      // common mean
  real<lower=0> sigmasq_a; // variance of between-subject random effects
  real<lower=0> sigmasq;  // variance of within-subject random effects
  vector[nsubjects] a;  // random effect for mean
}
model {
  // prior distributions
  mu ~ normal(0, sqrt(1000));
  sigmasq_a ~ inv_gamma(0.001, 0.001);
  sigmasq ~ inv_gamma(0.001, 0.001);

  // distribution of random effects for mean
  for(i in 1:nsubjects) a[i] ~ normal(0, sqrt(sigmasq_a));

  // data distribution
  for(i in 1:nobs) {
    yobs[i] ~ normal(mu + a[subject_id_obs[i]], sqrt(sigmasq));
  }
}
generated quantities {
  vector[nmiss] ymiss;
  for (i in 1:nmiss) {
    ymiss[i] = normal_rng(mu + a[subject_id_miss[i]], sqrt(sigmasq));
  }
}
"

# fixed effects model
fecode = "
data {
  int<lower=1> nobs; // number of observed data values
  int<lower=1> nmiss;  // number of missing data values
  int<lower=1> nsubjects; // number of subjects measured
  int subject_id_obs[nobs]; //subject index for yobs.  Must be int!
  int subject_id_miss[nmiss]; //subject index for ymiss.  Must be int!
  vector[nobs] yobs;  // vector of observed blood pressure measurements
}
parameters {
  real mu;      // common mean
  real<lower=0> sigmasq;  // variance of within-subject random effects
  vector[nsubjects] a;  // random effect for mean
}
model {
  // prior distributions
  mu ~ normal(0, sqrt(1000));
  sigmasq ~ inv_gamma(0.001, 0.001);

  // distribution of random effects for mean
  for(i in 1:nsubjects) a[i] ~ normal(0, sqrt(1000));

  // data distribution
  for(i in 1:nobs) {
      yobs[i] ~ normal(mu + a[subject_id_obs[i]], sqrt(sigmasq));
    }
  }
  generated quantities {
  vector[nmiss] ymiss;
  for (i in 1:nmiss) {
    ymiss[i] = normal_rng(mu + a[subject_id_miss[i]], sqrt(sigmasq));
  }
}
"

# common mean model
cmcode = "
data {
  int<lower=1> nobs; // number of observed data values
  int<lower=1> nmiss;  // number of missing data values
  int<lower=1> nsubjects; // number of subjects measured
  int subject_id_obs[nobs]; //subject index for yobs.  Must be int!
  int subject_id_miss[nmiss]; //subject index for ymiss.  Must be int!
  vector[nobs] yobs;  // vector of observed blood pressure measurements
}
parameters {
  real mu;      // common mean
  real<lower=0> sigmasq;  // variance of within-subject random effects
}
model {
  // prior distributions
  mu ~ normal(0, sqrt(1000));
  sigmasq ~ inv_gamma(0.001, 0.001);

  // data distribution
  for(i in 1:nobs) {
    yobs[i] ~ normal(mu, sqrt(sigmasq));
  }
}
generated quantities {
  real ytilde;
  ytilde = normal_rng(mu, sqrt(sigmasq));
}
"

# create the data list
blood_data = list(nobs = nobs, nmiss = nmiss,
                  nsubjects = nsubjects, yobs = yobs,
                  subject_id_obs = subject_id_obs,
                  subject_id_miss = subject_id_miss)

if (!file.exists("example_9_1_re_mod.rda")) {
  # draw samples from the models
  re_fit = stan(model_code = recode, data = blood_data,
                iter = 10, seed = 23)
  re_mod = stan_model(model_code = recode)
  save(re_mod, file = "example_9_1_re_mod.rda")
}
# draw samples from the models
re_fit = sampling(re_mod, data = blood_data,
                  iter = 5e4, seed = 23)


if (!file.exists("example_9_1_fe_mod.rda")) {
  # draw samples from the models
  fe_fit = stan(model_code = fecode, data = blood_data,
                iter = 10, seed = 24)
  fe_mod = stan_model(model_code = fecode)
  save(fe_mod, file = "example_9_1_fe_mod.rda")
}
# draw samples from the models
fe_fit = sampling(fe_mod, data = blood_data,
                  iter = 5e4, seed = 24)


if (!file.exists("example_9_1_cm_mod.rda")) {
  # draw samples from the models
  cm_fit = stan(model_code = cmcode, data = blood_data,
                iter = 10, seed = 25)
  cm_mod = stan_model(model_code = cmcode)
  save(cm_mod, file = "example_9_1_cm_mod.rda")
}
# draw samples from the models
cm_fit = sampling(cm_mod, data = blood_data,
                  iter = 5e4, seed = 25)

# summarize information from each model
re_sum = summary(re_fit, prob = c(0.025, 0.975))$summary
fe_sum = summary(fe_fit, prob = c(0.025, 0.975))$summary
cm_sum = summary(cm_fit, prob = c(0.025, 0.975))$summary

# order to match book
o = c(2, 3, 1, 4, 5)
# posterior information about ymiss
cbind(sub = subject_id_miss, meas = measurement_miss,
      round(re_sum[24:28, c("mean", "sd", "2.5%", "97.5%")], 1))[o, ]
cbind(sub = subject_id_miss, meas = measurement_miss,
      round(fe_sum[23:27, c("mean", "sd", "2.5%", "97.5%")], 1))[o, ]
round(cm_sum[3, c("mean", "sd", "2.5%", "97.5%"), drop = FALSE], 1)

library(bayesplot)

re_df = as.data.frame(re_fit)
fe_df = as.data.frame(fe_fit)
mcmc_intervals(re_df[,4:23])
mcmc_intervals(fe_df[,3:22])

# compute the posterior predictive sd
# of ytilde_31
# from p. 41 of BDA3
# var(ytilde | y) = E(sigmasq | y) + var(mu_i | y)
# standard deviation of ytilde | y using sample
# = sqrt(var(mu) + var(a[3]) + mean(sigmasq)
# compare SD(y_31 | y) for FE model and RE model
# FE model
sqrt(var(fe_df$mu) + var(fe_df$`a[3]`) + mean(fe_df$sigmasq))
# RE model
sqrt(var(re_df$mu) + var(re_df$`a[3]`) + mean(re_df$sigmasq))
sqrt(var(re_df$mu) + mean(re_df$sigmasq_a) + mean(re_df$sigmasq))

round(re_sum[c(1:3, 6), c("mean", "sd")], 1)
round(fe_sum[c(1:2, 5), c("mean", "sd")], 1)
round(cm_sum[1:2, c("mean", "sd"), drop = FALSE], 1)
