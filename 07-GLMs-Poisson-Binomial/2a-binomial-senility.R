library(rstan)
library(tidyr)
library(ggplot2)

# Example:  Example 3 of Chapter 7 of Bayesian Modeling
# Using Winbugs by Ntzoufras
#
# In Agresti (1990), a subset of the Wechsler Adult Intelligence Scale (WAIS)
# was given to 54 elderly people to help identify people with senility symptoms.
# The WAIS score is a discrete value between 0 and 20 and senility is a binary
# variable indicating whether a person is suffering senility symptoms.  Interest
# also lies in determining which WAIS scores correspond to an increased
# probability of senility symptoms.

# Data distribution: senility_i∼Binomial(1, pi_i),
# with logit(pi_i) = beta0 + beta_1 * wais
#
# Prior distributions: beta_j∼N(0, 10^3) (sigmasq = 10^3).

# Enter data manually
n = 54
wais = c(9, 13, 6, 8, 10, 4, 14, 8, 11, 7, 9, 7, 5, 14,  13, 16,
	10, 12, 11, 14, 15, 18, 7, 16, 9, 9, 11, 13, 15, 13, 10, 11, 6,
	17, 14, 19, 9, 11, 14, 10, 16, 10, 16, 14, 13, 13, 9, 15, 10, 11,
	12, 4, 14, 20)
senility  =  c(1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0)

wais_code = "
data {
  int<lower=1> n; // number of observations
  int y[n];       // indicator of senility
  vector[n] x;    // vector of covariates
}
parameters {
  real beta0;
  real beta1;
}
model {
  // prior distributions
  beta0 ~ normal(0, sqrt(1000));
  beta1 ~ normal(0, sqrt(1000));
  // data distribution
  y ~ bernoulli_logit(beta0 + x*beta1);
}
generated quantities {
  real exp_beta0;
  real exp_beta1;
  real half_prob;
  exp_beta0 = exp(beta0);
  exp_beta1 = exp(beta1);
  half_prob = -beta0/beta1;
}
"
# create the data list
senility_data = list(n = n, y = senility, x = wais)

# compile from model
if (!file.exists("wais_mod.rda")) {
  wais_fit = stan(model_code = wais_code,
                  data = senility_data,
                  iter = 50000, seed = 90)
  wais_mod = stan_model(model_code = wais_code)
  save(wais_mod, file = "wais_mod.rda",
       compress = "xz")
}

# draw samples from model
load("wais_mod.rda")
wais_fit = sampling(wais_mod, data = senility_data,
                    iter = 50000, seed = 90)

# check convergence with gelman-rubin statistics
summary(wais_fit)$summary[,"Rhat"]

# check convergence with trace plots
stan_trace(wais_fit, c("beta0", "beta1"))

# summary of fitted values
summary(wais_fit)$summary[c("beta0", "beta1"),]

# posterior means, medians, and 95% central posterior intervals
summary(wais_fit)$summary[c("beta0", "beta1"), c("mean", "50%", "2.5%", "97.5%")]

# plot of densities
stan_dens(wais_fit, par = c("beta0", "beta1"),
          separate_chains = TRUE)

# posterior means, medians, and 95% central posterior intervals
# of exponentiated parameters
summary(wais_fit)$summary[c("exp_beta0", "exp_beta1"), c("mean", "50%", "2.5%", "97.5%")]

# posterior means, medians, and 95% central posterior intervals
# of half probability
summary(wais_fit)$summary["half_prob", c("mean", "50%", "2.5%", "97.5%")]

# get posterior samples
wais_chain = as.data.frame(wais_fit)


# create posterior values of pi for different values of wais
logit.pi = wais_chain$beta0 +
	wais_chain$beta1 * matrix(0:20, byrow = TRUE, ncol = 21, nrow = nrow(wais_chain))
pi = exp(logit.pi)/(1 + exp(logit.pi))
# when wais is less than or equal to 7, the probability of senility is more than 50%
# for both the mean and median model
# when wais is less than or equal to 3, the probability of senility is more than 50%
# for the lower 0.025 model
# when wais is less than or equal to 9, the probability of senility is more than 50%
# for the higher 0.975 model

# create graph of mean posterior probability for senility
# versus wais score
# also plot median and credible intervals
# place horizontal bar of y = 0.5
plot(0:20, apply(pi, 2, mean), xlab = "wais score",
	ylab = "probability of senility", type = "l", ylim = c(0, 1))
# median lines
lines(0:20, apply(pi, 2, median), lty = 2)
# lower 0.025 quantile
lines(0:20, apply(pi, 2, quantile, prob = 0.025), col = "blue", lty = 3)
# upper 0.975 quantile
lines(0:20, apply(pi, 2, quantile, prob = 0.975), col = "blue", lty = 3)
abline(h = 0.5)
legend("topright", legend = c("mean", "median", "95% posterior interval"),
	col = c("black", "black", "blue"), lty = c(1, 2, 3, 3))
title("Probability of senility vs Wais score")

df = data.frame(wais = 0:20,
                mean = apply(pi, 2, mean),
                median = apply(pi, 2, median),
                q_0.025 = apply(pi, 2, quantile, prob = 0.025),
                q_0.975 = apply(pi, 2, quantile, prob = 0.975))
ggplot(df) + geom_line(aes(x = wais, y = mean)) + geom_line(aes(x = wais, y = median), lty = 2) + theme_bw()

dfa = data.frame(mean = apply(pi, 2, mean),
                 median = apply(pi, 2, median))
dfa = gather(dfa)
dfa$wais = rep(0:20, 2)
dfb = data.frame(q_025 = apply(pi, 2, quantile, prob = 0.025),
                 q_975 = apply(pi, 2, quantile, prob = 0.975))
dfb$wais = 0:20

ggplot() +
  geom_ribbon(data = dfb, aes(x = wais, ymin = q_025, ymax = q_975), fill = "grey90") +
  geom_line(data = dfa, aes(x = wais, y = value, lty = key))  +
  theme_bw() + ggtitle("Probability of senility vs Wais score") +
  ylab("probability of senility")

