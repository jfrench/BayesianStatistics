library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

#Example:  Soft drink delivery times (Chapter 5 of Bayesian Modeling
#Using Winbugs by Ntzoufras
#
#We are interested in estimation of the required time needed by
#each employee in a delivery system network to refill an automatic
#vending machine.  For this reason, a small quality assurance study
#was set up by an industrial engineer of the company.  As the
#response variable, the engineer considered the total service time
#(measured in minutes) of each machine, including its stocking with
#beverage products and any required maintenance or housekeeping.
#After examining the problem, the industrial engineer recommended
#two important variables that affect delivery time:  the number
#of cases of stocked products and the distance walked by the
#employee (measured in feet).
#
# Data distribution: y_i ~ N(x_i *beta, 1/tau).
# Prior distributions:
# beta ~ N(0, sigma^2 * c^2 * (X'X)^(-1))
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
mod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  matrix[n, 3] X; //covariates
  vector[3] mu0; //prior mean for beta
  cov_matrix[3] V; //V part of Zellner's g-prior
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
  vector[n] mu;  // mean of responses.  Temporary
  mu = X * beta;

  // prior distributions
  beta ~ multi_normal(mu0, sigmasq*csq*V);
  prec ~ gamma(0.01, 0.01);
  // data distribution
  for(i in 1:n) y[i] ~ normal(mu[i], sqrt(sigmasq));
}
generated quantities {
  real Rbsq;
  Rbsq = 1 - sigmasq / v;
}
"

dat1 = list(n = n, y = soda$Time, X = X, mu0 = c(0, 0, 0),
             V = solve(crossprod(X)), v = v, csq = n)

dat2 = list(n = n, y = soda$Time, X = X, mu0 = c(0, 0, 0),
            V = solve(crossprod(X)), v = v, csq = 100^2)

if (!file.exists("soda_gprior_mod.rda")) {
  # fit models using stan with for both g-priors
  soda_fit1 = stan(model_code = mod, data = dat1,
              iter = 10000, seed = 30)
  soda_fit2 = stan(model_code = mod, data = dat2,
               iter = 10000, seed = 31)
  # compile model
  soda_gprior_mod = stan_model(model_code = mod)
  # save model
  save(soda_gprior_mod, file = "soda_gprior_mod.rda",
       compress = "xz")
}
load(file = "soda_gprior_mod.rda")
# draw samples from the model
soda_gprior_fit1 = sampling(soda_gprior_mod, data = dat1,
                            iter = 10000, seed = 30)
soda_gprior_fit2 = sampling(soda_gprior_mod, data = dat2,
                            iter = 10000, seed = 31)

# combine results into a data frame
fits = rbind(cbind(as.data.frame(soda_gprior_fit1),
                   model = "model1"),
             cbind(as.data.frame(soda_gprior_fit2),
                   model = "model2"))

# restructure data frame for plotting
df = tidyr::pivot_longer(
  data = fits,
  cols = c("beta[1]", "beta[2]", "beta[3]"),
  names_to = "parameter"
)

# check convergence with gelman-rubin statistics
summary(soda_gprior_fit1)$summary[,"Rhat"]
summary(soda_gprior_fit2)$summary[,"Rhat"]

# density plots of betas
ggplot(df, aes(x = value)) + geom_density(aes(fill = model), alpha = 0.4) +
  facet_wrap(~ parameter, scales = "free") + theme_bw()

# 95% central posterior intervals and means
summarize(group_by(df, parameter, model),
          mean = mean(value),
          l = quantile(value, prob = 0.025),
          u = quantile(value, prob = 0.975)
          )

# comparison of Rbsq
ggplot(fits, aes(x = Rbsq)) + theme_bw() +
  geom_density(aes(fill = model), alpha = 0.4)
