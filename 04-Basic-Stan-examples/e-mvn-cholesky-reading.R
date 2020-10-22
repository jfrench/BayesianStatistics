library(rstan)
library(bayesplot)
library(mvtnorm)

### Chapter 7 Example from A First Course in Bayesian Statistical Methods
### by P. Hoff.

#A sample of 22 children are given reading comprehension tests before
#and after receiving a particular instructional method.  Each student
#will then have two scores denoting the pre- and post-instructional scores

# Each row of the data set is a two-dimensional vector representing a single
# case, distributed N(mu, Sigma).

# For priors:
# mu | Sigma ~ N(mu0, Sigma/k0)
# Sigma ~ Inv-Wishart(nu, L0^(-1))

# The exam was designed to give average score around 50 out of 100, so
# set mu0 = c(50, 50).  Since the true mean cannot be below 0 or above
# 100, we want the prior variances to keep mu inside this range with
# high probability. Assume k0 = 1 (a single prior observation on the
# Sigma scale). We want Sigma[1,1] = Sigma[2,2] to be approximately 625, so
# that the prior probability that the mean is outside [0, 100] is only
# 0.05.  Since the two exams measure similar things, we think the
# correlation between the means is around 0.5, so
# Sigma[1,2] = Sigma[2,1] = 0.5(625) = 312.5.
# Since E(Sigma) ~ Inv-Wishart(nu0, L0^(-1)) = L0/(nu0 - d - 1),
# choosing nu0 = 4 loosely centers samples of Sigma around L0.
# So we choose L0[1,1]=L0[2,2] = 625 and the diagonal elements equal to
# 312.5.

# load data
y = matrix(c(
  59, 77, 43, 39, 34, 46, 32, 26, 42, 38, 38, 43, 55, 68,
  67, 86, 64, 77, 45, 60, 49, 50, 72, 59, 34, 38, 70, 48,
  34, 55, 50, 58, 41, 54, 52, 60, 60, 75, 34, 47, 28, 48,
  35, 33), ncol = 2, byrow = TRUE)
# reformat
y = data.frame(pretest = y[,1], posttest = y[,2])

# Create model.  Notice the quotes
stanmod = "
data {
  int<lower=1> n; // number of trials
  vector[2] y[n]; // 2d array of pre- and post-test measurements
  vector[2] mu0; //prior mean for mu
  real<lower=0> k0; //prior number of observations for mu
  real<lower=0> nu0; // df for prior for Sigma
  cov_matrix[2] L0; // scale of prior for Sigma
}
parameters {
  vector[2] mu;
  cov_matrix[2] Sigma;
}
transformed parameters{
  matrix[2, 2] chol_Sigma;
  chol_Sigma = cholesky_decompose(Sigma);
}
model {
  // specify prior distributions for mu and Sigma
  mu ~ multi_normal(mu0, Sigma/k0);
  Sigma ~ inv_wishart(nu0, L0);

  // specify data distribution in terms of Cholesky
  // parameterization of the multivariate normal
  for(i in 1:n){
    y[i] ~ multi_normal_cholesky(mu, chol_Sigma);
  }
}
generated quantities {
  vector[2] ytilde; //pre- and post-test measurements for new students
  // vector[2] yrep[n];
  ytilde = multi_normal_cholesky_rng(mu, chol_Sigma);
  // much faster to generate yrep manually
  //for(i in 1:n){
  //  yrep[i] = multi_normal_cholesky_rng(mu, chol_Sigma);
  //}
}
"

L0 = matrix(c(625, 312.5, 312.5, 625), nrow = 2)

stan_dat = list(n = 22, y = y,
                 mu0 = c(50, 50), k0 = 1, nu0 = 4, L0 = L0)

# # reading_chol_fit model
# reading_chol_fit = stan(model_code = stanmod, data = stan_dat,
#             iter = 1000)
# compile model
# reading_chol_mod = stan_model(model_code = stanmod)
# # save model
# save(reading_chol_mod, file = "reading_chol_mod.rda", compress = "xz")
# load compiled model
load(file = "reading_chol_mod.rda")
# draw samples from the model
reading_chol_fit = sampling(reading_chol_mod, data = stan_dat, iter = 1000, chains = 4)

### quantiles of mu from original example
#           1%      25%      50%      75%      99%
# mu1 40.21531 45.37341 47.30597 49.26025 54.17412
# mu2 45.90558 51.56652 53.72375 55.86592 61.38130

# results should be similar
summary(reading_chol_fit, par = "mu", probs = c(0.01, 0.25, 0.5, 0.75, 0.99))$summary[,4:8]

# extract samples list from reading_chol_fit
samples = extract(reading_chol_fit)

# determine posterior probability that
# post-test mean greater than pre-test mean
mean(samples$mu[,2] > samples$mu[,1])

# determine probability that post-test score
# greater than pre-test score for new student
mean(samples$ytilde[,2] > samples$ytilde[,1])

# trace plots of results
posterior = as.array(reading_chol_fit)
mcmc_trace(posterior, pars = c("mu[1]", "mu[2]"))
mcmc_trace(posterior, regex_pars = "mu")

# density plots of results
mcmc_dens_overlay(posterior, pars = c("mu[1]", "mu[2]"))
mcmc_dens_overlay(posterior, regex_pars = "mu")

# scatter plot of mu samples
mcmc_scatter(posterior, pars = c("mu[1]", "mu[2]"))
mcmc_scatter(posterior, regex_pars = "mu")

# plot of acf of chains
stan_ac(reading_chol_fit, "mu[1]")
stan_ac(reading_chol_fit, "mu[2]")
stan_ac(reading_chol_fit, "Sigma[1,1]")
stan_ac(reading_chol_fit, "Sigma[1,2]")
stan_ac(reading_chol_fit, "Sigma[2,2]")

# perform posterior predictive checks
nrep = 50
# generate yrep
yrep = matrix(0, nrow = nrep, ncol = prod(dim(y)))
for (i in seq_len(nrep)) {
  # generate multivariate normals from posterior predictive distribution
  ysim = mvtnorm::rmvnorm(nrow(y), mean = samples$mu[i,], sigma = samples$Sigma[i,,])
  yrep[i, ] = c(ysim)
}

# ppc yrep pretest
ppc_dens_overlay(y[,1], yrep[1:50,1:22])
# ppc yrep posttest
ppc_dens_overlay(y[,2], yrep[1:50,23:44])
