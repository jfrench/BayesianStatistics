library(rstan)
library(loo) # compute waic and looic

# Example:  Oxygen uptake (Kuehl 2000)

# Twelve healthy men who did not exercise regularly were recruited
# to take part in a study of the effects of two different exercise
# regimens on oxygen uptake.  Six of the twelve men were randomly
# assigned to a 12-week flat-terrain running program, and the
# remaining six were assigned to a 12-week step aerobics program.
# The maximum oxygen uptake of each subject was measured
# (in liters per minutes) while running on an inclined treadmill,

# Let the aerobic group be group1 and the running group be group 2.

# We will fit the parallel lines model to this data

# Data distribution: y_i | beta, alpha, D, X ~ indep. N(mui, sigma^2)
# where mui = beta0 + beta1 * x + alpha2 * D2, where D2 is an
# indicator variable indicating whether an observation
# comes from the running group.

# Prior distributions:
# beta0, beta1, alpha2 ~ N(0, 10000)
# sigma^2 ~ Inv-Gamma(.01, .01)

# The indicator variable is 1 if a person is in the running group, 0 otherwise
D2 = c(0,0,0,0,0,0,1,1,1,1,1,1)
# The age of participant
x = c(23,22,22,25,27,20,31,23,27,28,22,24)

# The change in maximal oxygen uptake
y = c(-0.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,0.26,2.51)
# number of observations
n = length(y)
# sample variance of responses
v = var(y)

# Create model.  Notice the quotes
pl_mod = "
data {
  int<lower=1> n;  // number of observations
  vector[n] y;     // change in maximal oxygen uptake measurements
  vector[n] x;     // vector of covariates
  vector[n] d2;     // vector of indicators
  real<lower=0> v; // sample variance of y
}
parameters {
  real<lower=0> sigmasq;
  real beta0;
  real beta1;
  real alpha2;
}
transformed parameters {
  vector[n] mu;           // mean of observations
  mu = beta0 + beta1*x + alpha2*d2;
  //for(i in 1:n) mu[i] = beta0 + beta1*x[i] + alpha2*d2[i];
}
model {
  // prior distributions
  sigmasq ~ inv_gamma(.01, .01);
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  // data distribution
  for(i in 1:n) y[i] ~ normal(mu[i], sqrt(sigmasq));
}
generated quantities {
  real Rbsq;              // goodness-of-fit
  real log_lik[n];      // log likelihood of each observation
  Rbsq = 1 - sigmasq/v;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | mu[i], sqrt(sigmasq));
}
"

# Specify the data in R, using a list
# format compatible with Stan:
dat = list(n = n, y = y, x = x, d2 = D2, v = v)

# draw samples from the model
# oxygen_pl_mod = stan_model(model_code = pl_mod)
# save(oxygen_pl_mod, file = "oxygen_pl_mod.rda", compress = "xz")
load("oxygen_pl_mod.rda")
fit_oxygen_pl = stan(model_code = pl_mod, data = dat, iter = 5000,
            control = list(adapt_delta = 0.99), seed = 43)

# check convergence with gelman-rubin statistics
summary(fit_oxygen_pl)$summary[,"Rhat"]

# check convergence with trace plots
stan_trace(fit_oxygen_pl, c("beta0", "beta1", "alpha2", "sigmasq"))

# summary of fitted values
summary(fit_oxygen_pl)$summary[c("beta0", "beta1", "alpha2", "sigmasq"),]

# posterior means
summary(fit_oxygen_pl)$summary[c("beta0", "beta1", "alpha2", "sigmasq"),"mean"]

# 95% central posterior intervals
summary(fit_oxygen_pl)$summary[c("beta0", "beta1", "alpha2", "sigmasq"), c("2.5%", "97.5%")]

# plot of densities
stan_dens(fit_oxygen_pl, par = c("beta0", "beta1", "alpha2", "sigmasq"),
          separate_chains = TRUE)

# distribution of Rb^2
stan_dens(fit_oxygen_pl, "Rbsq") + xlim(c(0.5, 1))

# Plot with estimated regression lines (using means)
plot(x, y, xlab = "Age", ylab = "Maximum Oxygen Uptake", pch = D2 + 1,
	col = D2 + 3)
legend("topleft", pch = c(1, 2), col = c("green", "blue"), legend =
	c("aerobic", "running"))

# determine means of posteriors
pl_coef = summary(fit_oxygen_pl)$summary[c("beta0", "beta1", "alpha2"),"mean"]
abline(pl_coef[1], pl_coef[2], col = "green")
abline(pl_coef[1] + pl_coef[3], pl_coef[2], col = "blue")

# Create model.  Notice the quotes
sl_mod = "
data {
  int<lower=1> n;  // number of observations
  vector[n] y;     // change in maximal oxygen uptake measurements
  vector[n] x;     // vector of covariates
  vector[n] d2;    // vector of indicators
  real v;          // sample variance of y
}
parameters {
  real<lower=0> sigmasq;
  real beta0;
  real beta1;
  real alpha2;
  real delta2;
}
transformed parameters{
  vector[n] mu;   // mean of data distribution
  for (i in 1:n) {
    mu[i] = beta0 + beta1*x[i] + alpha2*d2[i] + delta2*x[i]*d2[i];
  }
}
model {
  // prior distributions
  sigmasq ~ inv_gamma(.01, .01);
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  alpha2 ~ normal(0, 100);
  delta2 ~ normal(0, 100);

  //data distribution
  for(i in 1:n) y[i] ~ normal(mu[i], sqrt(sigmasq));
}
generated quantities {
  real Rbsq;          // goodness-of-fit
  vector[n] log_lik;  // log likelihood of data
  Rbsq = 1 - sigmasq/v;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i]|mu[i], sqrt(sigmasq));
}
"

# draw samples from the model
# oxygen_sl_mod = stan_model(model_code = sl_mod)
# save(oxygen_sl_mod, file = "oxygen_sl_mod.rda", compress = "xz")
load("oxygen_sl_mod.rda")
fit_oxygen_sl = stan(model_code = sl_mod, data = dat, iter = 5000,
                     control = list(adapt_delta = 0.99), seed = 43)

# check convergence with gelman-rubin statistics
summary(fit_oxygen_sl)$summary[,"Rhat"]

# check convergence with trace plots
stan_trace(fit_oxygen_sl, c("beta0", "beta1", "alpha2", "delta2", "sigmasq"))

# summary of fitted values
summary(fit_oxygen_sl)$summary[c("beta0", "beta1", "alpha2", "delta2", "sigmasq"),]

# posterior means
pl_means = summary(fit_oxygen_sl)$summary[c("beta0", "beta1", "alpha2", "delta2", "sigmasq"),"mean"]
print(pl_means, digits = 2)

# 95% central posterior intervals
pl_cpi = summary(fit_oxygen_sl)$summary[c("beta0", "beta1", "alpha2", "delta2", "sigmasq"), c("2.5%", "97.5%")]
print(pl_cpi, digits = 2)

# plot of densities
stan_dens(fit_oxygen_sl, par = c("beta0", "beta1", "alpha2", "delta2", "sigmasq"),
          separate_chains = TRUE)

# distribution of Rb^2
stan_dens(fit_oxygen_sl, "Rbsq") + xlim(c(0.25, 1))

# Plot with estimated regression lines (using means)
plot(x, y, xlab = "Age", ylab = "Maximum Oxygen Uptake", pch = D2 + 1,
     col = D2 + 3)
legend("topleft", pch = c(1, 2), col = c("green", "blue"), legend =
         c("aerobic", "running"))

sl_coef = summary(fit_oxygen_sl)$summary[c("beta0", "beta1", "alpha2", "delta2"),"mean"]
abline(sl_coef[1], sl_coef[2], col = "green")
abline(sl_coef[1] + sl_coef[3], sl_coef[2] + sl_coef[4], col = "blue")

# extract log likelihoods
ll_pl = extract_log_lik(fit_oxygen_pl, merge_chains = FALSE)
ll_sl = extract_log_lik(fit_oxygen_sl, merge_chains = FALSE)

# compute waic on both fitted models
(waic_pl = waic(ll_pl))
(waic_sl = waic(ll_sl))

# compute relative efficiency of log likelihoods
r_eff_pl = exp(relative_eff(ll_pl))
r_eff_sl = exp(relative_eff(ll_sl))

# compute looic for each model
(looic_pl = loo(ll_pl,
                 r_eff = r_eff_pl))
(looic_sl = loo(ll_sl,
                 r_eff = r_eff_sl))

loo_compare(waic_pl, waic_sl)
loo_compare(looic_pl, looic_sl)
