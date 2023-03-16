library(rstan)
library(ggplot2)

# Example:  Evaluation of Tutors
# Chapter 5 - Bayesian Modeling Using Winbugs by Ntzoufras
#
# The director of a private school wishes to employ a new
# mathematics tutor.  For this reason, the ability of four
# candidates is examined using a small study.  A group of 25
# students was randomly divided into four classes (1 for each
# tutor.  In all classes, the same mathematical topic
# was taught for 2 hours per day for 1 week.
# After completing the short course, all
# students had to take the same test.  Their grades were
# recorded and compared.  The administrator wishes to employ
# the tutor whose students attained the higher performance
# for the given test.
#
# Data distribution: y_jk ~ N(mu0 + alpha_j, sigmasq)
# Prior distribution: mu0, alpha_j ~ N(0, 100^2)
# sigmasq ~ Inv-Gamma(0.01,0.01).
#
# We will use indicator variables when tutor_fitting this model and
# a corner constraint (alpha1 = 0) to tutor_fit this model.

grades = c(84, 58, 100, 51, 28, 89,
			97, 50, 76, 83, 45, 42, 83,
			64, 47, 83, 81, 83, 34, 61,
			77, 69, 94, 80, 55, 79)

tutor = rep(1:4, times = c(6, 7, 7, 6))

# create indicator variables
tutor2 = (tutor == 2) + 0
tutor3 = (tutor == 3) + 0
tutor4 = (tutor == 4) + 0

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  vector[n] tutor2; //indicator for tutor 2
  vector[n] tutor3; //indicator for tutor 3
  vector[n] tutor4; //indicator for tutor 4
}
parameters {
  real<lower=0> sigmasq;
  real mu0;
  real alpha2;
  real alpha3;
  real alpha4;
}
model {
  // create mean vector for use in data distribution
  vector[n] mu;

  //specify priors
  mu0 ~ normal(0.0, 100);
  alpha2 ~ normal(0.0, 100);
  alpha3 ~ normal(0.0, 100);
  alpha4 ~ normal(0.0, 100);
  sigmasq ~ inv_gamma(0.01, 0.01);

  // compute mu.  This is not stored since it's declared
  // in the model block
  mu = mu0 + alpha2*tutor2 + alpha3*tutor3 + alpha4*tutor4;

  // data distribution
  for(i in 1:n){
    y[i] ~ normal(mu[i], sqrt(sigmasq));
  }
}
"

stan_dat = list(n = length(grades), y = grades,
                tutor2 = tutor2, tutor3 = tutor3,
                tutor4 = tutor4)

if (!file.exists("tutor_mod.rda")) {
  # # tutor_fit = stan(model_code = stanmod, data = stan_dat, iter = 10000)
  tutor_mod = stan_model(model_code = stanmod)
  # save model
  save(tutor_mod, file = "tutor_mod.rda", compress = "xz")
}
load(file = "tutor_mod.rda")
# draw samples from the model
tutor_fit = sampling(tutor_mod, data = stan_dat,
                     iter = 10000)

# alternative example with specified starting values
# init_list = list(list(sigmasq = 1, mu0 = 9, alpha2 = 0, alpha3 = 9, alpha4 = 4),
#                  list(sigmasq = 2, mu0 = 25, alpha2 = 1, alpha3 = -5, alpha4 = 1),
#                  list(sigmasq = 3, mu0 = 50, alpha2 = 2, alpha3 = -19, alpha4 = 10),
#                  list(sigmasq = 4, mu0 = 75, alpha2 = 3, alpha3 = -10, alpha4 = 7))
#
# tutor_fit = sampling(tutor_mod, data = stan_dat, iter = 10000, chains = 4,
#                      init = init_list)

# summary of tutor_fitted values
summary(tutor_fit)$summary

# check convergence with gelman-rubin statistics
summary(tutor_fit)$summary[,"Rhat"]

# plot of densities
stan_dens(tutor_fit,
          par = c("mu0", "alpha2", "alpha3", "alpha4"),
          separate_chains = TRUE)

# 95% central posterior intervals
summary(tutor_fit)$summary[,c("2.5%", "97.5%")]

# posterior means
summary(tutor_fit)$summary[,"mean"]

# get samples
samples = rstan::extract(tutor_fit)
# create data frame of sample of means for
# each group
df = data.frame(mean = c(samples$mu0,
                samples$mu0 + samples$alpha2,
                samples$mu0 + samples$alpha3,
                samples$mu0 + samples$alpha4),
                tutor = rep(c("tutor 1", "tutor 2",
                              "tutor 3", "tutor 4"),
                            each = length(samples$mu)))

# parallel boxplots
ggplot(df) + geom_boxplot(aes(y = mean, x = tutor))
# parallel violin plots
ggplot(df) + geom_violin(aes(y = mean, x = tutor))
# overlaid density plots
ggplot(df) + geom_density(aes(x = mean, fill = tutor),
                          alpha = 0.4)

