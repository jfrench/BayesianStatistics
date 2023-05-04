# Example 9.3, Bayesian Modeling Using WinBUGS by Ntzoufras
# We consider modeling the number of times each survey
# participant had intercourse within the previous month.
# The data are accumulated by gender by Agresti (2002, pp.
# 569-70).  The samples means are 5.9 and 4.3 for males and
# females, respectively, while the variances are much higher
# (54.8 and 34.4, respectively).  The data are
# overdispersed!  A simple Poisson log-linear regression
# model would be inadequate.
#
# Poisson-gamma model
# Data distribution: Y_i ~ Poisson(lambda_gender * u)
# Prior distribution:
# u ~ Gamma(r_gender], r_gender])
# r_gender ~ Gamma(0.001, 0.001), gender = 0, 1
# lambda_gender ~ Gamma(0.001, 0.001), gender = 0, 1 
#
# observed data
nimble_data = list(
  y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12, 
        12, 12, 12, 12, 13, 13, 13, 15, 15, 15, 16, 16, 16, 20, 20, 20, 20, 20, 20, 20, 24, 25, 30, 30, 30, 50, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 17, 18, 20, 20, 20, 20, 20, 20, 22, 23, 25, 25, 25, 27, 30), 
  gender = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

# fit the poisson-gamma model using nimble
library(nimble)
pgcode = nimbleCode({
  for (i in 1:550) {
    # Poisson part
    y[i] ~ dpois(mu[i])
    # defining the mean of the Poisson
    mu[i] <- lambda[ gender[i] + 1 ] * u[i]
    # mixing distribution
    u[i] ~ dgamma( r[ gender[i] + 1 ] , r[ gender[i]+1 ]  )
  }
  for (j in 1:2) {
    # prior distributions
    lambda[j] ~ dgamma(0.001, 0.001)
    r[j] ~ dgamma(0.001, 0.001)
    # dispersion index
    di[j] <- (1 + lambda[j]/r[j])
    # assumed variance
    var[j] <- lambda[j]*di[j]
    # negative binomial probability
    p[j] <- r[j]/(r[j] + lambda[j])
  }
})

# Rmodel = nimbleModel(pgcode,
#                      data = nimble_data,
#                      inits = list(lambda = c(4, 5),
#                                   r = c(1, 2),
#                                   u = runif(550, 1, 2),
#                                   mu = 5 * runif(550, 2, 2)))
# mcmcspec = configureMCMC(Rmodel,
#                          monitors = c('r', 'lambda', 'mu', 'var'),
#                          enableWAIC = TRUE)
# Rmcmc = buildMCMC(mcmcspec)
# Cmodel = compileNimble(Rmodel)
# Cmcmc = compileNimble(Rmcmc, project = Rmodel)
# pg = runMCMC(Cmcmc, niter = 5e4, 
#              nburnin = 5e4/2,
#              nchains = 4,
#              samplesAsCodaMCMC = TRUE,
#              WAIC = TRUE,
#              summary = TRUE)
# 
# # create array from samples
# params_pg = as.array(pg$samples)
# params_pg = aperm(params_pg, c(1, 3, 2))
# 
# # see summary results
# summary_pg = pg$summary$all.chains[c(1:2, 553:556), ]
# save(summary_pg, file = "example_9_3_pg_summary.rda")
load(file = "example_9_3_pg_summary.rda")
summary_pg

## compute model fit information
# extract log likelihood
library(loo)
# # extract relevant mean-related parameters for each observation
# mu_pg = params_pg[,,3:552]
# mu_pg = matrix(c(mu_pg), ncol = dim(mu_pg)[3])
# 
# # compute ll for pg
# ll_pg = apply(params_pg[,,3:552], 1:2,
#               function(x) dpois(nimble_data$y, x, log = TRUE))
# ll_pg = aperm(ll_pg, c(2, 3, 1))
# 
# # compute waic
# waic_pg = waic(ll_pg)
# # compute effective sample size of liklihood
# r_eff_pg = relative_eff(exp(ll_pg))
# # compute looic (using effective sample size)
# looic_pg = loo(ll_pg, r_eff = r_eff_pg)
# save(waic_pg, looic_pg,
#      file = "example_9_3_pg_ic.rda")
load(file = "example_9_3_pg_ic.rda")
waic_pg
looic_pg

# # generate yrep for pg model
# s = sample(nrow(mu_pg), 1000)
# yrep_pg = t(apply(mu_pg[s,], 1, function(x) rpois(550, lambda = x)))
# save(yrep_pg, file = "example_9_3_yrep_pg.rda")
load(file = "example_9_3_yrep_pg.rda")

library(bayesplot)
ppc_hist(nimble_data$y, yrep_pg[1:8, ])
ppc_ecdf_overlay(nimble_data$y, yrep_pg)
ppc_intervals(nimble_data$y, yrep_pg)

# saveRDS(params_pg[,,c(1:2, 553:554)], file = "example_9_3_params_pg.RDS")
params_pg_s = readRDS("example_9_3_params_pg.RDS")
mcmc_trace(params_pg_s, regex_pars = c("r", "lambda"))
