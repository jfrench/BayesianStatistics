library(rstan)

# Example:  Soft drink delivery times (Chapter 5 of Bayesian
# Modeling Using Winbugs by Ntzoufras
#
# We are interested in estimation of the required time needed
# by each employee in a delivery system network to refill an
# automatic vending machine.  For this reason, a small
# quality assurance study was set up by an industrial
# engineer of the company.  As the response variable, the
# engineer considered the total service time (measured in
# minutes) of each machine, including its stocking with
# beverage products and any required maintenance or
# housekeeping. After examining the problem, the industrial
# engineer recommended two important variables that affect
# delivery time:  the number of cases of stocked products and
# the distance walked by the employee (measured in feet).
#
# Data distribution: y_i ~ N(x_i * beta, 1/tau).
# Prior distribution:
# beta ~ N(0, sigma^2 * c^2 * (X'X)^(-1)) with c^2 = 100^2
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

# Create X matrix using explanatory variables
X = cbind(1, soda$Cases, soda$Distance)
# Determine number of observations
n = length(soda$Cases)
# Sample variance of responses
v = var(soda$Time)

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  matrix[n, 3] X; //covariates
  vector[3] mu0; //prior mean for beta
  cov_matrix[3] V; //V part of Zellner's g-prior
  cov_matrix[n] I; //nxn identity matrix
  real<lower=0> csq; //constant for Zellner's g-prior
  real<lower=0> v; //sample variance
}
parameters {
  real<lower=0> prec;
  vector[3] beta;
}
transformed parameters{
  real<lower=0> sigmasq; //get sigmasq from the precision
  sigmasq = 1/prec;
}
model {
  //prior distributions
  beta ~ multi_normal(mu0, sigmasq*csq*V);
  prec ~ gamma(0.01, 0.01);
  // data distribution
  y ~ multi_normal(X*beta, sigmasq*I);
}
generated quantities {
  real Rbsq;
  Rbsq = 1 - sigmasq / v;
}
"

stan_dat = list(n = n, y = soda$Time,
                X = X, mu0 = c(0, 0, 0), V = solve(crossprod(X)),
                I = diag(n), csq = 100^2, v = v)

# fit model using stan with 4 chains
soda_g2_fit = stan(model_code = stanmod, data = stan_dat,
                   iter = 10000, chains = 4)

soda_g2_mod = stan_model(model_code = stanmod)
# # save model
save(soda_g2_mod, file = "soda_g2_mod.rda", compress = "xz")
load(file = "soda_g2_mod.rda")
# draw samples from the model
soda_g2_fit = sampling(soda_g2_mod, data = stan_dat, iter = 10000, chains = 4)

# plot of densities
stan_dens(soda_g2_fit, par = c("beta", "sigmasq"), separate_chains = TRUE)

# check convergence with gelman-rubin statistics
summary(soda_g2_fit)$summary[,"Rhat"]

# 95% central posterior intervals
summary(soda_g2_fit)$summary[,c("2.5%", "97.5%")]

# posterior means
summary(soda_g2_fit)$summary[,"mean"]

# distribution of Rb^2
stan_dens(soda_g2_fit, "Rbsq")