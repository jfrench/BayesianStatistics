# Example 9.4, Bayesian Modeling Using WinBUGS by Ntzoufras
# We consider modeling the odds ratio of cancer for
# smokers relative to nonsmokers.
#
# We have estimate or OR and CI information from 7 studies,
# with actual case counts from 3 other studies.
#
# Data distributions:
# log(or_k)|theta_k, sigmasq_k ~ N(theta_k, sigmasq_k), k = 1, ..., 7
# Y_ik ~ Binom(N_ik, invlogit(theta_ik)), k = 8, 9, 10, i = 1, 2
# i = 1 - smokers
# i = 2 - nonsmokers
# theta_ik = a_k + theta_k * I(i = 1)
# Prior distributions:
# theta_k | mu_theta, sigmasq_theta ~ N(mu_theta, sigmasq_theta)
# k = 1, ..., 10
# a_k ~ N(0, 1000), k = 8, 9, 10
# mu_theta ~ N(0, 1000)
# sigmasq_theta ~ IG(0.001, 0.001)

# observed odds ratios
or = c(3.89, 3.97, 3.88, 17.47, 5.35, 9.1, 3.41, 3.48, 33.10, 3.43)
# lower confidence limits
L = c(0.92, 2.2, 2.47, 14.24, 2.44, 5.57, 2.94)
# upper confidence limits
U = c(16.3, 7.16, 6.08, 21.43, 11.74, 14.86, 3.96)
# observed counts
# column 1 is cancers cases, 2 is non-cancer
Y = array(c(49, 33, 29958, 70186, 12, 89, 0, 118, 29, 171, 4, 81),
          dim = c(2, 2, 3))
N = matrix(0, 2, 3)
for (k in 1:3) {
  for (i in 1:2) {
    N[i,k] = Y[i,1,k] + Y[i,2,k];
  }
}
# cancer cases
# row 1 is smokers, row 2 is non-smokers
# column 1 is study 8, column 2 is study 9, etc.
cases = Y[,1,]

ma_data = list(K1 = 7, K2 = 3,
               logor = log(or),
               selogor = log(U/L)/(2 * 1.96),
               cases = cases, N = N)

ma_code = "
data {
  int<lower=1> K1;  // number of observations for CIs
  int<lower=1> K2;  // number of complete data
  real logor[K1 + K2];   // log of odds ratios
  real selogor[K1]; // estimated standard error of logor
  int cases[2, K2]; // number of cases in studies 8-10
  int N[2, K2];     // number of trial in studies 8-10
}
parameters {
  real mu_theta;         //mean of theta
  real sigmasq_theta;    //variance of theta
  vector[K1 + K2] theta; // mean of log odds ratio for each study
  vector[3] a;           // a random effect on the log odds of cancer for non-smokers
}
transformed parameters {
  real logit_p[2, K2];
  real p[2, K2];
  for (k in 1:K2){
    for (i in 1:2){
      logit_p[i,k] = a[k] + theta[K1 + k] * (i < 2);
      p[i,k] = inv_logit(a[k] + theta[K1 + k] * (i < 2));
    }
  }
}
model{
  // prior distribution
  for( k in 1:3) {
    a[k] ~ normal(0, sqrt(1000));
  }
  mu_theta  ~ normal(0, sqrt(1000));
  sigmasq_theta ~ inv_gamma(0.001, 0.001);

  for (k in 1:K1){
    logor[k] ~ normal(theta[k],  selogor[k] );
    theta[k] ~ normal(mu_theta, sqrt(sigmasq_theta));
  }

  for (k in 1:K2){
    for (i in 1:2){
      cases[i, k] ~ binomial_logit(N[i,k], logit_p[i,k]);
    }
    theta[K1+k] ~ normal(mu_theta, sqrt(sigmasq_theta));
  }
}
generated quantities {
  vector[K1 + K2] OR;
  real OR_mu;
  OR = exp(theta);
  OR_mu = exp(mu_theta);
}
"

# poisson lognormal starting values
init_fun_ma = function() {
  list(theta = rnorm(10, 0, 0.1),
       a = rnorm(3, 0, 0.1),
       mu_theta = rnorm(1, 0, 0.1),
       sigmasq_theta = runif(1, 0.04, 0.10)
  )
}

library(rstan)
library(bayesplot)

if (!file.exists("example_9_4.rda")) {
ma_mod = stan(model_code = ma_code, data = ma_data,
           iter = 5e5, seed = 11,
           init = init_fun_ma,
           warmup = floor(5e5/20*19),
           control = list(adapt_delta = 0.99))
save(ma_mod,
     file = "example_9_4.rda", compress = "bzip2")
}
load(file = "example_9_4.rda")

sum_ma = summary(ma_mod, pars = c("OR", "OR_mu"))$summary[,c("mean", "sd", "2.5%", "97.5%", "Rhat")]
round(sum_ma, 3)

ma_array = as.array(ma_mod)
# trace plot of various parameters
mcmc_trace(ma_array, regex_pars = c("mu_theta", "sigmasq_theta"))
mcmc_trace(ma_array[,,28:37]) # odds ratio

# posterior intervals for odds ratio of cancer for smokers
# and nonsmokers
mcmc_intervals(ma_array, regex_pars = "OR")

# posterior predictive check for odds ratios
# y = c(or, 3.48, 33.10, 3.43)
samples_ma = extract(ma_mod)
# yrep = samples_ma$OR
ppc_intervals(or, samples_ma$OR) +
  xlab("study") + ylab("odds ratio")

# posterior predictive check for cases
# for studies 8-10
sim_cases = function(idx) {
  p = samples_ma$p[idx,,]
  Ytemp = matrix(0, 2, 3)
  for (k in 1:3) {
    for (i in 1:2) {
      Ytemp[i, k] = rbinom(1, size = N[i,k],
                           prob = p[i,k])
    }
  }
  c(Ytemp)
}

set.seed(483)
# observed cases
cases = c(Y[,1,])
# replicated cases based on p
cases_rep = t(sapply(1:nrow(samples_ma$logit_p), sim_cases))
# posterior predictive intervals
ppc_intervals(cases, cases_rep) + ylab("cancer cases")
