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

library(rstan)         

# Enter data manually
n <- 54 
wais <- c(9, 13, 6, 8, 10, 4, 14, 8, 11, 7, 9, 7, 5, 14,  13, 16, 
	10, 12, 11, 14, 15, 18, 7, 16, 9, 9, 11, 13, 15, 13, 10, 11, 6,
	17, 14, 19, 9, 11, 14, 10, 16, 10, 16, 14, 13, 13, 9, 15, 10, 11, 
	12, 4, 14, 20)
senility  <-  c(1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0)

df = data.frame(wais, senility)
glmod_logit = stan_glm(senility ~ wais, data = df,
                       family = binomial(link = "logit"),
                       prior = normal(0, sqrt(1000)),
                       prior_intercept = normal(0, sqrt(1000)),
                       iter = 10000)

glmod_probit = stan_glm(senility ~ wais, data = df,
                         family = binomial(link = "probit"),
                         prior = normal(0, sqrt(1000)),
                         prior_intercept = normal(0, sqrt(1000)),
                         iter = 10000)

glmod_cloglog = stan_glm(senility ~ wais, data = df,
                       family = binomial(link = "cloglog"),
                       prior = normal(0, sqrt(1000)),
                       prior_intercept = normal(0, sqrt(1000)),
                       iter = 10000)

summary(glmod_logit)
summary(glmod_probit)
summary(glmod_cloglog)

# compare model fits
library(loo)
waic_l <- waic(glmod_logit)
waic_p <- waic(glmod_probit)
waic_c <- waic(glmod_cloglog)

loo_l <- loo(glmod_logit)
loo_p <- loo(glmod_probit)
loo_c <- loo(glmod_cloglog)

compare(waic_l, waic_p, waic_c)
compare(loo_l, loo_p, loo_c)




