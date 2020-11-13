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

df = data.frame(damage = damage, type = factor(type), bombload = bombload, airexp = airexp)
library(rstanarm)
glmod = stan_glm(damage ~ type * airexp + type * bombload, 
                data = df, 
                family = poisson(link = "log"),
                prior = normal(0, 10^1.5), 
                prior_intercept = normal(0, 10^1.5))
plot(glmod)

library(bayesplot)
# density of y vs yrep
pp_check(glmod, plotfun = "dens_overlay") + xlim(c(0, 20))
# histograms of y and yrep
pp_check(glmod, plotfun = "hist")
# plot of test statistic
pp_check(glmod, plotfun = "stat", stat = "sd") + xlim(c(0, 6))
# boxplots of y vs yrep
pp_check(glmod, plotfun = "boxplot")
# scatterplot of y vs mean(y_rep)
pp_check(glmod, plotfun = "scatter_avg")
# scatterplot of y vs realization of yrep
pp_check(glmod, plotfun = "scatter", nreps = 9)
# posterior predictive intervals
pp_check(glmod, plotfun = "intervals")
# average errors versus an x variable
pp_check(glmod, plotfun = "error_scatter_avg_vs_x", x = "bombload")
