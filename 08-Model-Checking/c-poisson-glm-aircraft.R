library(rstan)
library(lattice)
library(plyr)
library(dplyr)
library(ggplot2)
library(lattice)

# Example:  Example 1 of Chapter 7 of Bayesian Modeling
# Using Winbugs by Ntzoufras
#
# Montgomery et al. (2006) examined the number of aircraft
# damages in 30 strike missions during the Vietnam war.  The
# data we will analyze consist of 30 observations with the
# variables: damage: the number of damaged locations on the
# aircraft type: binary variable which indicates the type of
# plane (0 for A4; 1 for A6) bombload: the aircraft bomb
# load in tons airexp: the total months of aircrew
# experience
#
# Data distribution: damage_i ~ Poisson(lambda_i), with
# log(lambda_i) = beta0 + beta_1 * type + beta2 * bombload +
# beta3 * airexp Prior distribution: beta_j ~ N(0, 10^3)
# (sigmasq = 10^3).

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
tba_code_pp_check = "
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
transformed parameters {
  real loglambda[n];
  for(i in 1:n) {
    loglambda[i] = beta0 + beta1 * type[i] + beta2 * bombload[i] + beta3 * airexp[i];
  }
}
model {
  // prior distributions
  beta0 ~ normal(0.0, sqrt(10^3));
  beta1 ~ normal(0.0, sqrt(10^3));
  beta2 ~ normal(0.0, sqrt(10^3));
  beta3 ~ normal(0.0, sqrt(10^3));
  // data distribution
  y ~ poisson_log(loglambda);
}
generated quantities {
    real log_lik[n];
    real lik[n];
    for(i in 1:n) {
      log_lik[i] = poisson_log_lpmf(y[i] | loglambda[i]);
      lik[i] = exp(log_lik[i]);
    }
}
"

# specify data
aircraft_data = list(n = n, y = damage, type = type,
                     bombload = bombload, airexp = airexp)
if (!file.exists("aircraf_mod_pp_check.rda")) {
  aircraft_mod_fit =
    stan(model_code = tba_code_pp_check,
         data = aircraft_data,
         iter = 10, chains = 2)
  aircraft_mod_pp_check =
    stan_model(model_code = tba_code_pp_check)
  save(aircraft_mod_pp_check,
       file = "aircraft_mod_pp_check.rda", compress = "xz")
}
# load all compiled models
load("aircraft_mod_pp_check.rda")
fit_aircraft_check = sampling(aircraft_mod_pp_check,
                              data = aircraft_data,
                              iter = 5000, chains = 2)

# check convergence with gelman-rubin statistics
summary(fit_aircraft_check)$summary[,"Rhat"]

# plot of densities
stan_dens(fit_aircraft_check,
          par = c("beta0", "beta1", "beta2", "beta3"),
          separate_chains = TRUE)

chains = as.data.frame(fit_aircraft_check)

# get lambda values for each observation
lambda = exp(chains[,5:34])
# 30 x nsim yrep
yrep = apply(lambda, 1, function(x) rpois(n = 30, x))

# count and organize counts
yrep_freq = apply(yrep, 2, plyr::count)
yrep_freqdf = dplyr::bind_rows(yrep_freq)
y_freqdf = plyr::count(data.frame(damage))

# compare boxplots of yrep counts to observed counts
# too many 0s, not enough 1s
ggplot(yrep_freqdf, aes(x = x, y = freq)) +
  geom_boxplot(aes(group = cut_width(x, 1))) +
  geom_path(data = y_freqdf,
            mapping = aes(x = damage, y = freq),
            col = "blue")

# plot ppo vs y
ppo = colMeans(chains[,65:94])
# ppo for each value of damage
plot(ppo ~ damage)
# relative frequency of each value of damage
points(x = y_freqdf$damage,
       y = y_freqdf$freq/sum(y_freqdf$freq),
       type = "h")
