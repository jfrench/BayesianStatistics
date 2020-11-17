library(rstan)
library(loo)

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

# Data distribution: damage_i  ~ indep. Poisson(lambda_i),
# with log(lambda_i) = beta0 + beta_1 * type + beta2 * bombload + beta3 * airexp
# Prior distribution: beta_j ~ N(0, 10^3) (sigmasq = 10^3).

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

mod_t = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta1;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta1 * type[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta1 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"

mod_b = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta2;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta2 * bombload[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta2 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"

mod_a = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta3;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta3 * airexp[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta3 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"

mod_tb = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta1;
  real beta2;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta1 * type[i] + beta2 * bombload[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta1 ~ normal(0.0, sqrt(10^3));
  beta2 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"

mod_ta = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta1;
  real beta3;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta1 * type[i] + beta3 * airexp[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta1 ~ normal(0.0, sqrt(10^3));
  beta3 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"

mod_ba = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta2;
  real beta3;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta2 * bombload[i] + beta3 * airexp[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta2 ~ normal(0.0, sqrt(10^3));
  beta3 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"

mod_tba = "
data {
  int<lower=1> n;  // number of observations
  int y[n];        // number of damaged locations.
  // very important to declare as int!
  real type[n];    // type of aircraft
  real bombload[n];// bombload of aircraft
  real airexp[n];  // aircrew experience
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
}
transformed parameters{
  real mu[n];
  for(i in 1:n) mu[i] = exp(beta0 + beta1 * type[i] + beta2 * bombload[i] + beta3 * airexp[i]);
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta1 ~ normal(0.0, sqrt(10^3));
  beta2 ~ normal(0.0, sqrt(10^3));
  beta3 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson(mu);
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | mu[i]);
}
"
# model data
aircraft_data = list(n = n, y = damage, type = type, bombload = bombload,
                     airexp = airexp)

# compile, save the models
# aircraft_mod_t = stan_model(model_code = mod_t)
# save(aircraft_mod_t, file = "aircraft_mod_t.rda", compress = "xz")
# aircraft_mod_b = stan_model(model_code = mod_b)
# save(aircraft_mod_b, file = "aircraft_mod_b.rda", compress = "xz")
# aircraft_mod_a = stan_model(model_code = mod_a)
# save(aircraft_mod_a, file = "aircraft_mod_a.rda", compress = "xz")
# aircraft_mod_tb = stan_model(model_code = mod_tb)
# save(aircraft_mod_tb, file = "aircraft_mod_tb.rda", compress = "xz")
# aircraft_mod_ta = stan_model(model_code = mod_ta)
# save(aircraft_mod_ta, file = "aircraft_mod_ta.rda", compress = "xz")
# aircraft_mod_ba = stan_model(model_code = mod_ba)
# save(aircraft_mod_ba, file = "aircraft_mod_ba.rda", compress = "xz")
# aircraft_mod_tba = stan_model(model_code = mod_tba)
# save(aircraft_mod_tba, file = "aircraft_mod_tba.rda", compress = "xz")

# load all compiled models
load("aircraft_mod_t.rda")
load("aircraft_mod_b.rda")
load("aircraft_mod_a.rda")
load("aircraft_mod_tb.rda")
load("aircraft_mod_ta.rda")
load("aircraft_mod_ba.rda")
load("aircraft_mod_tba.rda")

# sampling from compiled models
fit_aircraft_mod_t = sampling(aircraft_mod_t, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)
fit_aircraft_mod_b = sampling(aircraft_mod_b, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)
fit_aircraft_mod_a = sampling(aircraft_mod_a, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)
fit_aircraft_mod_tb = sampling(aircraft_mod_tb, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)
fit_aircraft_mod_ta = sampling(aircraft_mod_ta, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)
fit_aircraft_mod_ba = sampling(aircraft_mod_ba, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)
fit_aircraft_mod_tba = sampling(aircraft_mod_tba, data = aircraft_data, iter = 10000, refresh = 0, chain = 2)

# compute and waic and looic
waic_t = waic(extract_log_lik(fit_aircraft_mod_t))
waic_b = waic(extract_log_lik(fit_aircraft_mod_b))
waic_a = waic(extract_log_lik(fit_aircraft_mod_a))
waic_tb = waic(extract_log_lik(fit_aircraft_mod_tb))
waic_ta = waic(extract_log_lik(fit_aircraft_mod_ta))
waic_ba = waic(extract_log_lik(fit_aircraft_mod_ba))
waic_tba = waic(extract_log_lik(fit_aircraft_mod_tba))

loo_t = loo(extract_log_lik(fit_aircraft_mod_t))
loo_b = loo(extract_log_lik(fit_aircraft_mod_b))
loo_a = loo(extract_log_lik(fit_aircraft_mod_a))
loo_tb = loo(extract_log_lik(fit_aircraft_mod_tb))
loo_ta = loo(extract_log_lik(fit_aircraft_mod_ta))
loo_ba = loo(extract_log_lik(fit_aircraft_mod_ba))
loo_tba = loo(extract_log_lik(fit_aircraft_mod_tba))

# compare results
loo_compare(waic_t, waic_b, waic_a, waic_tb,
            waic_ta, waic_ba, waic_tba)

loo_compare(loo_t, loo_b, loo_a, loo_tb,
            loo_ta, loo_ba, loo_tba)
