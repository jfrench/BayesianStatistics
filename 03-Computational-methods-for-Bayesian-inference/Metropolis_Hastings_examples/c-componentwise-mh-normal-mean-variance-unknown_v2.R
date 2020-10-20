# load necessary packages
library(bayesutils)
library(coda)

### Create componentwise metropolis hastings
### sampler for normal(mu, sigma^2) sampling distribution
### with mu and sigma^2 unknown.

# Data distribution: y1, y2, ..., yn |mu, sigma^2 ~ iid N(mu, sigma^2).
# Prior distribution:
# mu | sigma^2 ~ N(mu0, sigma^2/k0) w/ mu0 = 1.9 and k0 = 1
# sigma^2 ~ Inv-Chisq(nu0, sigmasq0) with nu0 = 1 and sigmasq0 = 0.01.

# Proposal distributions:
# mu_star | mu^(t-1) ~ N(mu^(t-1), mustar_sigma^2) w/ mustar_sigma^2 = 0.2^2
# sigmasq_star ~ Inv-Chisq(df_sigmasqstar, scale_sigmasqstar^2) with
# df_sigmasqstar = 2 and scale_sigmasqstar = 0.1.

### Midge Example
### Chapter 5 of A First Course in Bayesian Statistical Methods by PD Hoff

# Grogan and Wirth (1981) provide data on the wing length in millimeters
# of nine members of a species of midge (small, two-winged flies).
# From these nine measurements we wish to make inference on the population mean
# theta.

# Studies of other populations suggest that the true mean should be around 1.9 mm
# with a standard deviation of 0.1.  However, this population may be different
# from the others, so we choose k0 and vu0 = 1 so that the prior distributions
# are only weakly centered around these estimates from other populations.

#Set prior parameters
mu0 = 1.9
k0 = 1
nu0 = 1
sigmasq0 = 0.01

#define number of simulations, data, sample size, sample st dev, sample mean
B = 100000
y = c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
n = length(y)

#Initial value for sigma
s = sd(y)
#Initial value for mu
m = mean(y)

cmh = function(B, start, jump_parm1, jump_parm2) {
  mustar_sigma = jump_parm1
  df_sigmasqstar = jump_parm2[1]
  scale_sigmasqstar = jump_parm2[2]

  # first column for mu, second column for sigmasq
  theta = matrix(0, nrow = B + 1, ncol = length(start))
  theta[1, ] = start

  # get current values of mu and sigmasq
  mu_current = start[1]
  sigmasq_current = start[2]

  for (i in 2:(B + 1)) {
    # component for mu

    # proposal new value for mu (based on current value of mu)
    mu_star = rnorm(1, mean = mu_current, sd = mustar_sigma)
    # print(paste("mu", mu_star)) used for debugging
    # what are the current values during a crash?

    # notice where we use the updated mu and where we do not
    # the sigmasq value just uses the current value of sigmasq
    num_logr = sum(dnorm(y, mean = mu_star, sd = sqrt(sigmasq_current), log = TRUE)) +
      dnorm(mu_star, mean = mu0, sd = sqrt(sigmasq_current/k0), log = TRUE) +
      bayesutils::dinvchisq(sigmasq_current, df = nu0, scale = sigmasq0, log = TRUE) -
      dnorm(mu_star, mean = mu_current, sd = mustar_sigma, log = TRUE)
    den_logr = sum(dnorm(y, mean = mu_current, sd = sqrt(sigmasq_current), log = TRUE)) +
      dnorm(mu_current, mean = mu0, sd = sqrt(sigmasq_current/k0), log = TRUE) +
      bayesutils::dinvchisq(sigmasq_current, df = nu0, scale = sigmasq0, log = TRUE) -
      dnorm(mu_current, mean = mu_star, sd = mustar_sigma, log = TRUE)

    logr = num_logr - den_logr

    # if we accept the new value of mu (mu_star), then update mu_current,
    # otherwise do nothing
    if (log(runif(1)) <= min(logr, 0)) {
      mu_current = mu_star
    } # otherwise do nothing

    # component for sigmasq
    # use mu for current iteration
    # change value of sigmasq in numerator and denominator

    sigmasq_star = bayesutils::rinvchisq(1, df = df_sigmasqstar, scale = scale_sigmasqstar)
    # print(sigmasq_star) used for debugging
    # what are the current values during a crash?

    # notice where we use the updated mu and where we do not
    # the sigmasq value just uses the previous iteration sigmasq
    num_logr = sum(dnorm(y, mean = mu_current, sd = sqrt(sigmasq_star), log = TRUE)) +
      dnorm(mu_current, mean = mu0, sd = sqrt(sigmasq_star/k0), log = TRUE) +
      bayesutils::dinvchisq(sigmasq_star, df = nu0, scale = sigmasq0, log = TRUE) -
      bayesutils::dinvchisq(sigmasq_star, df = df_sigmasqstar, scale = scale_sigmasqstar, log = TRUE)
    den_logr = sum(dnorm(y, mean = mu_current, sd = sqrt(sigmasq_current), log = TRUE)) +
      dnorm(mu_current, mean = mu0, sd = sqrt(sigmasq_current/k0), log = TRUE) +
      bayesutils::dinvchisq(sigmasq_current, df = nu0, scale = sigmasq0, log = TRUE) -
      bayesutils::dinvchisq(sigmasq_current, df = df_sigmasqstar, scale = scale_sigmasqstar, log = TRUE)
    logr = num_logr - den_logr

    # if we accept the new value of sigmasq (sigmasq_star), then update
    # sigmasq_current, otherwise do nothing
    if (log(runif(1)) <= min(logr, 0)) {
      sigmasq_current = sigmasq_star
    }

    # after updating all parameters, we must store them
    theta[i, ] = c(mu_current, sigmasq_current)
  }
  return(theta)
}

# run 3 chains with different starting values
# the sigmasq parameter doesn't mix very well.
B = 100000 # run with 10000 and 100000
# mu_sigma = .2
# df_sigmasq = 2
scale_sigmasq = .01

chain1 = cmh(B, start = c(m, s^2), jump_parm1 = .2, jump_parm2 = c(2, scale_sigmasq))
chain2 = cmh(B, start = c(1, 0.5^2), jump_parm1 = .2, jump_parm2 = c(2, scale_sigmasq))
chain3 = cmh(B, start = c(3, 0.05^2), jump_parm1 = .2, jump_parm2 = c(2, scale_sigmasq))

# values to keep
keep = seq(B/2 + 1, B + 1, by = 1)
# combine chains into mcmc list.  Discard half as warmup.
mc = mcmc.list(mcmc(chain1[keep, ]),
               mcmc(chain2[keep, ]),
               mcmc(chain3[keep, ]))

# summarize results of M-H algorithm
summary(mc)

# trace plot for the 3 chains
traceplot(mc)

# acf plot
autocorr.plot(mc)

# effective sample size
effectiveSize(mc)

# Calculate scale reduction factor
gelman.diag(mc, autoburnin = FALSE)

# Plot Gelman statistic
gelman.plot(mc, autoburnin = FALSE)

# Heidelberg & Welch
heidel.diag(mc) # if the halfwidth test fails, extend the chain(s)

# Raftery-Lewis
# like to see dependence factor less than 6.  If not the case,
# we'll need more samples
raftery.diag(mc)

# Geweke
# z scores for a test of equality (should be between -2 and 2) if converged
geweke.diag(mc)
