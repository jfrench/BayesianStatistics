# Model checking for Newcomb light data
# Data distribution: y_1, y_2, ..., y_n | mu, sigma^2 ~ iid N(mu, sigma^2)
# Prior distribution: p(mu, log(sigma)) propto 1

# load additional Bayesian distributions
library(bayesutils)
# convenient posterior predictive checks
library(bayesplot)
# compute leave-one-out predictive checks
library(loo)


### Example from Section 6.3 Estimating the Speed of Light
data(newcomb, package = "MASS") # load data

B = 100000 # number of simulations

# length, mean, and sd of data
n = length(newcomb)
m = mean(newcomb)
s = sd(newcomb)

### posterior interval for mu (simulation)
# sample from posterior of sigmasq
sigmasq = rinvchisq(B, df = n - 1, scale = s^2)
# sample from conditional posterior of mu
mu = rnorm(B, mean = m, sd = sqrt(sigmasq/n))

# store nyrep samples of yrep
# from from posterior predictive distribution
nyrep = 200
yrep = matrix(0, nrow = nyrep, ncol = n)

# rename for convenience
y = newcomb
# sample 66 observations from posterior predictive distribution
for (i in seq_len(nyrep)) {
  yrep[i, ] = rnorm(66, mean = mu[i], sd = sqrt(sigmasq[i]))
}

# compare minimum of observed data to minimum from
# replicated samples

# minimum of replicated sample
mins = apply(yrep, 1, min)
# estimated p-value
(sum(mins <= min(y)) + 1)/(length(mins) + 1)

# histogram comparing T(y) and T(yrep)
ppc_stat(y, yrep, stat = "min")

# look at asymmetry of distribution
# by comparing order statistics to
# samples of posterior mean
d_sim = d_obs = numeric(nrow(yrep))
sort_y = sort(newcomb)

for (i in 1:nrow(yrep)) {
  thetai = mu[i]
  sort_yrep = sort(yrep[i, ])
  d_sim[i] = abs(sort_yrep[61] - mu[i]) -
         abs(sort_yrep[6] - mu[i])
  d_obs[i] = abs(sort_y[61] - mu[i]) -
    abs(sort_y[6] - mu[i])
}
# estimated posterior predictive p-value
(sum(d_sim >= d_obs) + 1)/(length(d_sim) + 1)

# compare observed and simulated discrepancy measures
plot(d_sim ~ d_obs, xlab = "observed discrepancy", ylab = "simulated discrepancy")
abline(0, 1)

# compare histograms of y and yrep
ppc_hist(y, yrep[sample(nyrep, 8), ])
# compare boxplots of y and yrep
ppc_boxplot(y, yrep[sample(nyrep, 8), ])
# compare densities of y and yrep
ppc_dens_overlay(y, yrep[sample(nyrep, 20), ])
# compare ecdfs of y and yrep
ppc_ecdf_overlay(y, yrep[sample(nyrep, 20) , ])
# compare histograms of y - yrep
ppc_error_hist(y, yrep[sample(nyrep, 9) , ])
# compare scatterplots of y vs yrep
ppc_scatter(y, yrep[sample(nyrep, 9) , ])
# compare scatterplots of y vs y - yrep
ppc_error_scatter(y, yrep[sample(nyrep, 9) , ])

# marginal predictive checks
# comparison of observed data and 90% predictive cpis
ppc_intervals(y, yrep)
ppc_ribbon(y, yrep)

nyrep = 10000
yrep = loglik_yrep = matrix(0, nrow = nyrep, ncol = n)
for (i in seq_len(nyrep)) {
  yrep[i, ] = rnorm(n, mean = mu[i], sd = sqrt(sigmasq[i]))
  loglik_yrep[i,] = dnorm(yrep[i, ],
                          mean = mu[i],
                          sd = sqrt(sigmasq[i]),
                          log = TRUE)
}

# compute relative effective MCMC sample size
# divided by the total sample size for likelihood
r_eff = relative_eff(exp(loglik_yrep),
                     chain_id = rep(1, nyrep))
# compute leave-one-out information
loo_info = loo(loglik_yrep, r_eff = r_eff, save_psis = TRUE)

# construct leave-one-out prediction intervals
ppc_loo_intervals(y, yrep,
                  psis_object = loo_info$psis_object)

# construct leave-one-out quantiles
ppc_loo_pit_qq(y, yrep,
               lw = weights(loo_info$psis_object))

# ppo vs y
PPO = colMeans(exp(loglik_yrep))
plot(PPO ~ y, ylab = "PPO")
dy = density(y)
# scale density to match scale of PPO
dy$y = dy$y/max(dy$y)*(max(PPO) - min(PPO)) + min(PPO)
lines(dy)

