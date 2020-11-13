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

library(rstanarm)         

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

df = data.frame(damage, type, bombload, airexp)

# compile and sample each model
mod_t = stan_glm(damage ~ type, family = poisson(),
                 data = df,
                 prior = normal(0, sqrt(1000)),
                 prior_intercept = normal(0, sqrt(1000)),
                 iter = 10000)

mod_b = stan_glm(damage ~ bombload, family = poisson(),
                 data = df,
                 prior = normal(0, sqrt(1000)),
                 prior_intercept = normal(0, sqrt(1000)),
                 iter = 10000)

mod_a = stan_glm(damage ~ airexp, family = poisson(),
                 data = df,
                 prior = normal(0, sqrt(1000)),
                 prior_intercept = normal(0, sqrt(1000)),
                 iter = 10000)

mod_tb = stan_glm(damage ~ type + bombload, 
                  family = poisson(),
                  data = df,
                  prior = normal(0, sqrt(1000)),
                  prior_intercept = normal(0, sqrt(1000)),
                  iter = 10000)

mod_ta = stan_glm(damage ~ type + airexp, 
                  family = poisson(),
                  data = df,
                  prior = normal(0, sqrt(1000)),
                  prior_intercept = normal(0, sqrt(1000)),
                  iter = 10000)

mod_ba = stan_glm(damage ~ bombload + airexp, 
                  family = poisson(),
                  data = df,
                  prior = normal(0, sqrt(1000)),
                  prior_intercept = normal(0, sqrt(1000)),
                  iter = 10000)

mod_tba = stan_glm(damage ~ type + bombload + airexp, 
                  family = poisson(),
                  data = df,
                  prior = normal(0, sqrt(1000)),
                  prior_intercept = normal(0, sqrt(1000)),
                  iter = 10000)

# compute and waic and looic
library(loo)
waic_t <- waic(mod_t)
waic_b <- waic(mod_b)
waic_a <- waic(mod_a)
waic_tb <- waic(mod_tb)
waic_ta <- waic(mod_ta)
waic_ba <- waic(mod_ba)
waic_tba <- waic(mod_tba)

loo_t <- loo(mod_t)
loo_b <- loo(mod_b)
loo_a <- loo(mod_a)
loo_tb <- loo(mod_tb)
loo_ta <- loo(mod_ta)
loo_ba <- loo(mod_ba)
loo_tba <- loo(mod_tba)

compare(waic_t, waic_b, waic_a, waic_tb, 
        waic_ta, waic_ba, waic_tba)

compare(loo_t, loo_b, loo_a, loo_tb, 
        loo_ta, loo_ba, loo_tba)


