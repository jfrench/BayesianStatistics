# Example:  Example 1 of Chapter 7 of Bayesian Modeling Using Winbugs
# by Ntzoufras
#
# Montgomery et al. (2006) examined the number of aircraft damages in
# 30 strike missions during the Vietnam war.  The data we will analyze
# consist of 30 observations with the variables: damage: the number of
# damaged locations on the aircraft type: binary variable which
# indicates the type of plane (0 for A4; 1 for A6) bombload: the
# aircraft bomb load in tons airexp: the total months of aircrew
# experience

# Data distribution: damage_i ∼ Poisson(lambda_i),
# with log(lambda_i) = beta0 + beta_1 * type + beta2 * bombload + beta3 * airexp
# Prior distribution: beta_j ∼ N(0, 10^3) (sigmasq = 10^3).

library(rstan)

# Enter data manually
damage = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 1, 1, 1, 2,
	3, 1, 1, 1, 2, 0, 1, 1, 2, 5, 1, 1, 5, 5, 7)
type = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
	 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
bombload = c(4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8,
	 7, 7, 7, 10, 10, 10, 12, 12, 12, 8, 8, 8, 14, 14, 14)
airexp = c(91.5, 84, 76.5, 69, 61.5, 80, 72.5, 65, 57.5, 50,
	 103, 95.5, 88, 80.5, 73, 116.1, 100.6, 85, 69.4, 53.9,
	 112.3, 96.7, 81.1, 65.6, 50, 120, 104.4, 88.9, 73.7, 57.8)

n = length(damage)
aircraft_code = "
data {
  int<lower=1> n;      // number of observations
  int<lower=0> y[n];   // number of damaged locations.
                       // very important to declare as int!
  int<lower=0> type[n];// type of aircraft
  real bombload[n];    // bombload of aircraft
  real airexp[n];      // aircrew experience
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
}
model {
  real logmu[n];
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta1 ~ normal(0.0, sqrt(10^3));
  beta2 ~ normal(0.0, sqrt(10^3));
  beta3 ~ normal(0.0, sqrt(10^3));
  // data distribution
  for(i in 1:n) logmu[i] = beta0 + beta1 * type[i] + beta2 * bombload[i] + beta3 * airexp[i];
  y ~ poisson_log(logmu);
}
generated quantities {
  real exp_beta0;
  real exp_beta1;
  real exp_beta2;
  real exp_beta3;
  exp_beta0 = exp(beta0);
  exp_beta1 = exp(beta1);
  exp_beta2 = exp(beta2);
  exp_beta3 = exp(beta3);
}
"

# Specify the data in R, using a list format compatible with STAN:
aircraft_data = list(n = n, y = damage, type = type, bombload = bombload,
            airexp = airexp)

# draw samples from the model
# aircraft_mod = stan_model(model_code = aircraft_code)
# save(aircraft_mod, file = "aircraft_mod.rda", compress = "xz")
load("aircraft_mod.rda")
# draw samples from the stan model
fit_aircraft= sampling(aircraft_mod, data = aircraft_data,
                       iter = 10000, seed = 101, chains = 2)

# check convergence with gelman-rubin statistics
summary(fit_aircraft)$summary[,"Rhat"]

# check convergence with trace plots
stan_trace(fit_aircraft, c("beta0", "beta1", "beta2", "beta3"))

# summary of fit_aircraft output
summary(fit_aircraft)$summary[c("beta0", "beta1", "beta2", "beta3"),]

# posterior means and 95% central posterior intervals
summary(fit_aircraft)$summary[c("beta0", "beta1", "beta2", "beta3"), c("mean", "2.5%", "97.5%")]

# plot of densities
stan_dens(fit_aircraft, par = c("beta0", "beta1", "beta2", "beta3"),
          separate_chains = TRUE)

# posterior means and 95% central posterior intervals
# of exponentiated parameters
summary(fit_aircraft)$summary[c("exp_beta0", "exp_beta1", "exp_beta2", "exp_beta3"),
                     c("mean", "2.5%", "97.5%")]

#mean, median, low, high profile predictor values
mean.x = apply(cbind(bombload, airexp), 2, mean)
median.x = apply(cbind(bombload, airexp), 2, median)
low.x = c(min(bombload), max(airexp))
high.x = c(max(bombload), min(airexp))

#mean, median, low, high profile posterior values.  Separate calculations by type
chains = as.data.frame(fit_aircraft)

low.a4 = with(chains, exp(beta0 + low.x[1] * beta2 + low.x[2] * beta3))
med.a4 = with(chains, exp(beta0 + median.x[1] * beta2 + median.x[2] * beta3))
mean.a4 = with(chains, exp(beta0 + mean.x[1] * beta2 + mean.x[2] * beta3))
high.a4 = with(chains, exp(beta0 + high.x[1] * beta2 + high.x[2] * beta3))
low.a6 = low.a4 * chains$exp_beta1
med.a6 = med.a4 * chains$exp_beta1
mean.a6 = mean.a4 * chains$exp_beta1
high.a6 = high.a4 * chains$exp_beta1

# determine mean of posterior of different pofiles
mean(low.a4)
mean(med.a4)
mean(mean.a4)
mean(high.a4)
mean(low.a6)
mean(med.a6)
mean(mean.a6)
mean(high.a6)

# plot profiles
plot(density(low.a4))
plot(density(med.a4))
plot(density(mean.a4))
plot(density(high.a4))

plot(density(low.a6))
plot(density(med.a6))
plot(density(mean.a6))
plot(density(high.a6))
