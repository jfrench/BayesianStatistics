# Example:  Example 1 of Chapter 9 of Bayesian Modeling
# Using Winbugs by Ntzoufras
#
# Consider a problem related to blood pressure measurements.  We consider
# the repeated measurements of blood pressure from 20 health individuals.
# Our aim is to estimate both the between-subject and within-subject variability.
#
# We assume that is y_ij = mu + a_i + epsilon_ij
# We consider a hierarchical model
# a_i ~ N(0, sigmasq_a), epsilon_ij ~ N(0, sigmasq)
# 
# Equivalently, Y_ij ~ N(mu_i, sigma^2) with mui = mu + a_i
# and a_i ~ N(0, sigmasq_a)
#
# We assume low information prior distributions for the regression coefficients, 
# mu ~ N(0, 1000), sigmasq ~ IG(0.001, 0.001), sigmasq_a ~ IG(0.001, 0.001)

library(rstan)         

# Enter data manually
n = 20
y = matrix(c(108, 98, 91, 94, 93, 96, 104, 99, 99, 97, 95, 
             98, 93, 97, 99, 96, 90, 100, 92, 95, 101, 89, 
             97, 97, 97, 100, 96, 95, 106, 100, 100, 98, 
             90, 99, 88, 98, 92, 92, 100, 101), 
           ncol = 2, byrow = TRUE)

bloodmod = "
data {
  int<lower=1> n;  // number of observations
  matrix[n, 2] y;  // matrix of blood pressure measurements
}
parameters {
  real mu;      // common mean
  real<lower=0> sigmasq_a; // variance of between-subject random effects
  real<lower=0> sigmasq;  // variance of within-subject random effects
  vector[n] a;  // random effect for mean
}
model {
  // prior distributions
  mu ~ normal(0, sqrt(1000));
  sigmasq_a ~ inv_gamma(0.001, 0.001);
  sigmasq ~ inv_gamma(0.001, 0.001);
  
  // distribution of random effects for mean
  for(i in 1:n) a[i] ~ normal(0, sqrt(sigmasq_a));

  // data distribution
  for(i in 1:n) {
    for(j in 1:2) {
      y[i,j] ~ normal(mu + a[i], sqrt(sigmasq));
    }
  }
}
generated quantities {
  real<lower=0> tvar; // total variance
  real<lower=-1, upper=1> cor; //correlation
  real<lower=0> sigma; 
  real<lower=0> sigma_a;
  tvar = sigmasq + sigmasq_a;
  cor = sigmasq_a/tvar;
  sigma = sqrt(sigmasq);
  sigma_a = sqrt(sigmasq_a);
}
"

# create the data list
blood_data = list(n = n, y = y)

# # draw samples from the model
# blood_samples = stan(model_code = bloodmod, data = blood_data, 
#                       iter = 5e4, seed = 23, refresh = 0)
# save(blood_samples, file = "example_9_1.rda")
load(file = "example_9_1.rda")

# check convergence
summary(blood_samples)$summary[,"Rhat"]

# posterior summary
post_sum = summary(blood_samples, probs = c(0.025, 0.975))$summary
round(post_sum[c("mu", "sigmasq", "sigmasq_a", "tvar", "cor", 
           "sigma", "sigma_a"), c("mean", "sd", "2.5%", "97.5%")],
      2)

# Convert coda-object coda_samples to matrix object for easier handling.
chains = as.data.frame(blood_samples)

# display median, 50% posterior interval, 
# 95% posterior interval for random effects
library(bayesplot)
mcmc_intervals(chains[, 4:23], prob_outer = 0.95)
