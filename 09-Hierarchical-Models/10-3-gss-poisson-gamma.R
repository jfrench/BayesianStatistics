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

# observed data
gss = list(
  n = 550, 
  y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12, 
        12, 12, 12, 12, 13, 13, 13, 15, 15, 15, 16, 16, 16, 20, 20, 20, 20, 20, 20, 20, 24, 25, 30, 30, 30, 50, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 17, 18, 20, 20, 20, 20, 20, 20, 22, 23, 25, 25, 25, 27, 30), 
  gender = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

# poisson-gamma model
pg_code = "
data {
  int<lower=1> n;  // number of observations
  int y[n];  // count of sexual intercourse in previous month
  int<lower=0, upper=1> gender[n];  // male = 0, female = 1
}
parameters {
  vector[n] u;
  vector[2] lambda;
  vector[2] r;
}
transformed parameters {
  vector[n] mu;
  for (i in 1:n) {
    mu[i] = lambda[gender[i] + 1] * u[i];
  }
}
model {
  // prior distributions
  // different parameter for each gender
  for (j in 1:2) {
    lambda[j] ~ gamma(0.001, 0.001);
    r[j] ~ gamma(0.001, 0.001);
  }
  // data/mixture distributions
  for(i in 1:n) {
    u[i] ~ gamma(r[gender[i] + 1], r[gender[i] + 1]);
    y[i] ~ poisson(mu[i]);
  }
}
generated quantities {
  vector[n] log_lik; //log likelihood of observations
  vector[n] yrep;  //replicated data
  for(i in 1:n) {
    log_lik[i] = poisson_lpmf(y[i] | mu[i]);
    yrep[i] = poisson_rng(mu[i]);
  }
}
"

# the init functions below generate appropriate starting
# values for each chain.  Otherwise, impossible values can 
# be generated, causing Stan to fail

# poisson gamma starting values
init_fun_pg = function() { 
  list(lambda = rexp(2, 0.1), 
       r = runif(2, 0.2, 0.8), 
       u = rexp(gss$n, 0.1)) 
} 

# poisson lognormal starting values
library(rstan)
pg_mod = stan(model_code = pg_code, data = gss,
              iter = 5e1, seed = 99, init = init_fun_pg,
              chains = 4,
              control = list(adapt_delta = 0.9))

# pln_mod <- stan(model_code = pln, data = gss,
#                iter = 5e4, seed = 1004, 
#                init = init_fun2)
# 
# nb_mod <- stan(model_code = nb, data = gss,
#                iter = 5e4, seed = 7,
#                init = init_fun3, 
#                control = list(adapt_delta = 0.9))
# 
# # fit the poisson-gamma model using nimble
# nimble_data = gss
# nimble_data$n = NULL
# pgcode = nimbleCode({
#   for (i in 1:550) { 
#     # Poisson part
#     y[i] ~ dpois(mu[i])
#     # defining the mean of the Poisson
#     mu[i] <- lambda[ gender[i] + 1 ] * u[i]
#     # mixing distribution 
#     u[i] ~ dgamma( r[ gender[i] + 1 ] , r[ gender[i]+1 ]  )
#   }	
#   for (j in 1:2) {
#     # prior distributions
#     lambda[j] ~ dgamma(0.001, 0.001)
#     r[j] ~ dgamma(0.001, 0.001)
#     # dispersion index
#     di[j] <- (1 + lambda[j]/r[j])
#     # assumed variance
#     var[j] <- lambda[j]*di[j]
#     # negative binomial probability
#     p[j] <- r[j]/(r[j] + lambda[j])
#   }
# })
# 
# Rmodel = nimbleModel(pgcode, 
#                      data = nimble_data,
#                      inits = list(lambda = c(4, 5),
#                                   r = c(1, 2),
#                                   u = runif(550, 1, 2),
#                                   mu = 5 * runif(550, 2, 2)))
# mcmcspec <- configureMCMC(Rmodel, 
#                           monitors = c('r', 'lambda', 'u', 'mu'))
# Rmcmc <- buildMCMC(mcmcspec)
# Cmodel <- compileNimble(Rmodel)
# Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# pg_samples <- runMCMC(Cmcmc, niter = 5e4, nchains = 4)

# library(abind)
# pg_fit = abind(pg_samples$chain1,
#                pg_samples$chain2,
#                pg_samples$chain3,
#                pg_samples$chain4, along = 3)
# pg_fit = aperm(pg_fit, c(1, 3, 2))[25001:50000,,]
# pln_fit = as.array(pln_mod)
# nb_fit = as.array(nb_mod)
# save(pg_samples, pln_mod, nb_mod, 
#      file = "example_9_3.rda", compress = "bzip2")
# 
# # save posterior samples for simple parameters
# pgs = pg_fit[,,c(1:2, 553:554)]
# plns = pln_fit[,,551:554]
# nbs = nb_fit[,,1:4]
# 
# # save these for plotting/other purposes
# save(pgs, plns, nbs, 
#      file = "example_9_3_param_samples.rda", 
#      compress = "bzip2")
load(file = "example_9_3_param_samples.rda")
#
# s = sample(100000, 50)
# 
# extract relevant mean-related parameters for each observation
# mu_pg = pg_fit[,,3:552]
# mu_pg = matrix(c(mu_pg), ncol = dim(mu_pg)[3])[s, ]
# 
# logmu_pln = pln_fit[,,555:(555 + 549)]
# logmu_pln = matrix(c(logmu_pln),
#                    ncol = dim(logmu_pln)[3])[s, ]
# r_nb = nb_fit[,,3:4]
# r_nb = matrix(c(r_nb), ncol = dim(r_nb)[3])[s, ]
# r_nb = r_nb[, gss$gender + 1]
# lambda_nb = nb_fit[,, 5:554]
# lambda_nb = matrix(c(lambda_nb),
#                    ncol = dim(lambda_nb)[3])[s, ]
# p_nb = r_nb/(r_nb + lambda_nb)
#
# # generate yrep for each model 
# yrep_pg = t(apply(mu_pg, 1, function(x) rpois(550, lambda = mu_pg)))
# yrep_pln = t(apply(logmu_pln, 1, function(x) {
#   rpois(550, lambda = exp(logmu_pln))
#   }))
# yrep_nb = t(sapply(seq_len(nrow(p_nb)), function(i) {
#   rnbinom(550, size = r_nb[i, ], prob = p_nb[i,])
# }))
# save(yrep_pg, yrep_pln, yrep_nb, file = "example_9_3_yrep.rda")
load(file = "example_9_3_yrep.rda")

ppc_hist(gss$y, yrep_pg[1:8, ])
ppc_hist(gss$y, yrep_pln[1:8, ])
ppc_hist(gss$y, yrep_nb[1:8, ])

library(bayesplot)
mcmc_trace(pgs, regex_pars = c("r", "lambda"))
mcmc_trace(plns, regex_pars = c("beta", "sigmasq_e"))
mcmc_trace(nb_fit, regex_pars = c("r", "beta"))

# library(loo)
# ll_pg = apply(pg_fit[,,3:552], 1:2, 
#               function(x) dpois(gss$y, x, log = TRUE))
# ll_pg = aperm(ll_pg, c(2, 3, 1))
# ll_pln = apply(pln_fit[,,555:(555 + 549)], 1:2, 
#                function(x) dpois(gss$y, exp(x), log = TRUE))
# ll_pln = aperm(ll_pln, c(2, 3, 1))
# ll_nb = apply(nb_fit[,,3:554], 1:2, 
#               function(x) {
#                 r = x[1:2]
#                 rg = r[gss$gender + 1]
#                 mu = x[-(1:2)]
#                 p = rg/(rg + mu)
#                 dnbinom(gss$y, size = rg, prob = p, log = TRUE)
#               })
# ll_nb = aperm(ll_nb, c(2, 3, 1))
# 
# # compute waic for each model
# waic_pg = waic(ll_pg)
# waic_pln = waic(ll_pln)
# waic_nb = waic(ll_nb)
# 
# # relative efficiencies of exp(log_lik) for each model
# r_eff_pg = relative_eff(exp(ll_pg))
# r_eff_pln = relative_eff(exp(ll_pln))
# r_eff_nb = relative_eff(exp(ll_nb))
# 
# # compute looic for each model
# looic_pg = loo(ll_pg, r_eff = r_eff_pg)
# looic_pln = loo(ll_pln, r_eff = r_eff_pln)
# looic_nb = loo(ll_nb, r_eff = r_eff_nb)
# save(waic_pg, waic_pln, waic_nb,
#      looic_pg, looic_pln, looic_nb,
#      file = "ic_example_9_3.rda")
load(file = "ic_example_9_3.rda")

compare(waic_pg, waic_pln, waic_nb)
compare(looic_pg, looic_pln, looic_nb)


