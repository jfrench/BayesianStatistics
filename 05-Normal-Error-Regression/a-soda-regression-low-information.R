library(rstan)
library(bayesplot)

# Example:  Soft drink delivery times (Chapter 5 of Bayesian
# Modeling Using Winbugs by Ntzoufras
#
# We are interested in estimation of the required time
# needed by each employee in a delivery system network to
# refill an automatic vending machine.  For this reason, a
# small quality assurance study was set up by an industrial
# engineer of the company.  As the response variable, the
# engineer considered the total service time (measured in
# minutes) of each machine, including its stocking with
# beverage products and any required maintenance or
# housekeeping. After examining the problem, the industrial
# engineer recommended two important variables that affect
# delivery time:  the number of cases of stocked products
# and the distance walked by the employee (measured in
# feet).
#
# Data distribution: y_i ~ N(x_i^T*beta, 1/tau) for i = 1, 2, ..., n
# Prior distributions:
# beta_j ~ N(0, 10^4) (sigmasq = 10^4)
# tau ~ Gamma(0.01,0.01).

# load and format data set
soda = matrix(c(
  16.68,	7,	560,
  11.5,	3,	220,
  12.03,	3,	340,
  14.88,	4,	80,
  13.75,	6,	150,
  18.11,	7,	330,
  8,	2,	110,
  17.83,	7,	210,
  79.24,	30,	1460,
  21.5,	5,	605,
  40.33,	16,	688,
  21,	10,	215,
  13.5,	4,	255,
  19.75,	6,	462,
  24,	9,	448,
  29,	10,	776,
  15.35,	6,	200,
  19,	7,	132,
  9.5,	3,	36,
  35.1,	17,	770,
  17.9,	10,	140,
  52.32,	26,	810,
  18.75,	9,	450,
  19.83,	8,	635,
  10.75,	4,	150),
  ncol = 3, byrow = TRUE)
soda = as.data.frame(soda)
colnames(soda) = c("Time", "Cases", "Distance")

# obtain the number of observations
n = length(soda$Cases)
# obtain the sample variance of the data
v = var(soda$Time)

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  vector[n] cases; //covariate
  vector[n] dist; //covariate
  real<lower=0> v; // sample variance of y
}
parameters {
  real<lower=0> prec; // tau
  real beta0;
  real beta1;
  real beta2;
}
transformed parameters{
  real<lower=0> sigma; //get sigma from the precision
  sigma = sqrt(1/prec);
}
model {
  //specify priors
  beta0 ~ normal(0.0, 100);
  beta1 ~ normal(0.0, 100);
  beta2 ~ normal(0.0, 100);
  prec ~ gamma(0.01, 0.01);

  // data distribution
  for(i in 1:n){
    y[i] ~ normal(beta0 + beta1*cases[i] + beta2*dist[i], sigma);
  }
}
generated quantities {
  real<lower=0> sigmasq; //get sigmasq from the precision
  real Rbsq;
  sigmasq = 1/prec;
  Rbsq = 1 - sigmasq / v;
}
"

stan_dat = list(n = n, y = soda$Time, cases = soda$Cases,
                 dist = soda$Distance, v = v)

# if compiled model doesn't already exist,
# compile the model, sample from model,
# returns object of class stan
# save model
if (!file.exists("soda_mod.rda")) {
  # soda_fit = stan(model_code = stanmod, data = stan_dat,
  #                 iter = 10000, chains = 4)
  # save(soda_fit, file = "soda_fit.rda", compress = "xz")
  # load(file = "soda_fit.rda")
  soda_mod = stan_model(model_code = stanmod)
  # save model
  save(soda_mod, file = "soda_mod.rda", compress = "xz")
}
load(file = "soda_mod.rda")
# # draw samples from the model
soda_fit = sampling(soda_mod, data = stan_dat,
                    iter = 10000, chains = 4)

# summary of soda_fitted values
summary(soda_fit)$summary

# plot of densities
stan_dens(soda_fit, par = c("beta0", "beta1", "beta2", "sigmasq"),
          separate_chains = TRUE)

# check convergence with gelman-rubin statistics
summary(soda_fit)$summary[,"Rhat"]

# 95% central posterior intervals
summary(soda_fit)$summary[,c("2.5%", "97.5%")]

# posterior means
summary(soda_fit)$summary[,"mean"]

# distribution of Rb^2
stan_dens(soda_fit, "Rbsq")

# how to draw from the predictive distribution of each
# observation

# extracting the MCMC samples from the soda_fit object
samples = extract(soda_fit)
ncycles = length(samples[[1]])
# each row of yrep is a sample from the pp distribution
yrep = matrix(0, ncol = nrow(soda), nrow = ncycles)
for (i in seq_len(nrow(soda))) {
  mui = samples$beta0 + samples$beta1 * soda$Cases[i] + samples$beta2 * soda$Distance[i]
  yrep[, i] = rnorm(ncycles, mean = mui, sd = samples$sigma)
}

# approximate posterior predictive density for an observation
# with the same covariates as observation 1
plot(density(yrep[,1]))

# posterior predictive check
y = soda$Time
ppc_hist(y, yrep[sample(1:ncycles, 8),])

# scatterplot of y vs yrep
ppc_scatter(y, yrep[10:18,])

# scatterplot of y vs average yrep
ppc_scatter_avg(y, yrep)

